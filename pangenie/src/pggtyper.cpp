#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include <sys/stat.h>
#include <mutex>
#include <thread>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include "kmercounter.hpp"
#include "jellyfishreader.hpp"
#include "jellyfishcounter.hpp"
#include "emissionprobabilitycomputer.hpp"
#include "copynumber.hpp"
#include "variantreader.hpp"
#include "uniquekmercomputer.hpp"
#include "hmm.hpp"
#include "commandlineparser.hpp"
#include "timer.hpp"
#include "threadpool.hpp"
#include "pathsampler.hpp"
#include "cli/load/load_graph.cpp"
#include "cli/stats.hpp"
#include "cli/query.hpp"
#include "cli/config/config.hpp"
#include "cli/load/load_annotated_graph.hpp"
#include "graph/annotated_dbg.hpp"

using namespace std;


bool ends_with (string const &filename, string const &ending) {
    if (filename.length() >= ending.length()) {
        return (0 == filename.compare (filename.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

void check_input_file(string &filename) {
	// check if file exists and can be opened
	ifstream file(filename);
	if (!file.good()) {
		stringstream ss;
		ss << "File " << filename << " cannot be opened." << endl;
		throw runtime_error(ss.str());
	}
	// make sure file is not compressed
	if (ends_with(filename, ".gz")) {
		stringstream ss;
		ss << "File " << filename << " seems to be gzip-compressed. PanGenie requires an uncompressed file." << endl;
		throw runtime_error(ss.str());
	}
}

struct UniqueKmersMap {
	mutex kmers_mutex;
	map<string,
            map<string, vector<UniqueKmers*> >
            > unique_kmers;
	map<string, double> runtimes;
};

struct Results {
	mutex result_mutex;
	map<string,
            map<string, vector<GenotypingResult>>
            > result;
	map<string, double> runtimes;
};

void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, KmerCounter* read_kmer_counts, VariantReader* variant_reader, map<string,ProbabilityTable>* probs, UniqueKmersMap* unique_kmers_map, map<string,size_t> kmer_coverage) {
	Timer timer;
	UniqueKmerComputer kmer_computer(genomic_kmer_counts, read_kmer_counts, variant_reader, chromosome, kmer_coverage);
	map<string,std::vector<UniqueKmers*> > unique_kmers;
        // this needs to be map or vector of vectors for each sample
	kmer_computer.compute_unique_kmers(&unique_kmers, probs);
	// store the results
	{
		lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
                for(auto uniqPair:unique_kmers) {
                  string sampleName=uniqPair.first;
                  if (unique_kmers_map->unique_kmers.find(sampleName) ==
                      unique_kmers_map->unique_kmers.end()) {
                    unique_kmers_map->unique_kmers[sampleName] =
                        map<string, vector<UniqueKmers *>>();
                  }

                  auto tmp = pair<string, vector<UniqueKmers *>>(
                      chromosome, move(unique_kmers[sampleName]));
                  unique_kmers_map->unique_kmers[sampleName].insert(tmp);
                }
                unique_kmers_map->runtimes[chromosome]+=timer.get_total_time();
	}
}

void run_genotyping(string chromosome, vector<UniqueKmers*>* unique_kmers, ProbabilityTable* probs, bool only_genotyping, bool only_phasing, long double effective_N, vector<unsigned short>* only_paths, Results* results) {
	Timer timer;
	// construct HMM and run genotyping/phasing
	HMM hmm(unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, false, effective_N, only_paths, false);
	// store the results
	{
		lock_guard<mutex> lock_result (results->result_mutex);
                cout<<"Here"<<endl;
                if (results->result.find("sample") == results->result.end()) {
                  results->result["sample"] =  map<string, vector<GenotypingResult>>();
                }
		// combine the new results to the already existing ones (if present)
		if (results->result["sample"].find(chromosome) == results->result["sample"].end()) {
			results->result["sample"].insert(pair<string, vector<GenotypingResult>> (chromosome, hmm.move_genotyping_result()));
		} else {
			// combine newly computed likelihoods with already exisiting ones
			size_t index = 0;
			vector<GenotypingResult> genotypes = hmm.move_genotyping_result();
			for (auto likelihoods : genotypes) {
				results->result["sample"].at(chromosome).at(index).combine(likelihoods);
				index += 1;
			}
		}
		// normalize the likelihoods after they have been combined
		for (size_t i = 0; i < results->result["sample"].at(chromosome).size(); ++i) {
			results->result["sample"].at(chromosome).at(i).normalize();
		}

                if (results->runtimes.find(chromosome) == results->runtimes.end()) {
                  results->runtimes.insert(pair<string,double>(chromosome, timer.get_total_time()));
                } else {
                  results->runtimes[chromosome] += timer.get_total_time();
                }
                cout<<"Out"<<endl;
	}

}

bool ends_with (string const &full_string, string const ending) {
	if (full_string.size() >= ending.size()) {
		return (0 == full_string.compare(full_string.size() - ending.size(), ending.size(), ending));
	} else {
		return false;
	}
}

int main (int argc, char* argv[])
{
//
//  cout<<"Hello World"<<endl;
//  string graph_path=argv[1];
//  string annotation_path=argv[2];
//
//  std::shared_ptr<mtg::graph::DeBruijnGraph> graph=mtg::cli::load_critical_dbg(graph_path);
//  mtg::cli::print_stats(*graph);
//  auto config = std::make_unique<mtg::cli::Config>();
//  config->query_counts=true;
//  config->infbase=graph_path;
//  config->infbase_annotators.push_back(annotation_path);
//  std::unique_ptr<mtg::graph::AnnotatedDBG> anno_graph=mtg::cli::initialize_annotated_dbg(graph, *config);;
//
//  mtg::cli::QuerySequence query;
//  query.id=0;
//  query.name="ERCC-00165";
//  query.sequence="GATATGCGTTACGTGAGTCTGATAGCAGTTCACTACCTGGATATCTGATCCACTAGCTCGATCATGCTCACCCATAGTTTATCTGCATCACTCGTACTGAAATGCTCACATCGCAGGTAGAGCAGCATCGTAGAGCGTCAAGCTGCATCCTAGCGTCATGAGTCATAGTACCTCATGCTCACGTGATCTACCCTAGCTGACCGCTAATGACGGCAGTGCAACCTGAGATACCGACGGCATACTGTCGTCAACGTCAGGCAATGTGTCCGAACGGCGAGCTACGTCGCCTCACGGAGTAATCGCGTCCCTCTAGGTATAGTGCCGTCGGTTCAGGTCATATGTCGCGGGTTCTGCACATATCACGGACGTATCGCTATCAGACGGACGCTCTCGGACCTAAACCGTAGCTCTCGGCAAGATCGTCCTCGTCTCGAATATAGCGCCCTAGTGCTGCAAATGTCACCGCTATCTCGTAAGGGGTCCGTCTGTTGAGTTAGGCCTCCTCTCGTTGGATGTGAGCTCGGTTGCTTGGATGGTGCAGCTTACTTCGCGTACCTGCTGTTTGCATCAGTCCTCTGCATCTATAATCGCGTATCTCTCTCTAGTAGACCATATAGCCATCTAAGCGCTCGATATTCCACCTAAGTGGCGCCTATTGAACTAAGTGGCAGCCGAATGGACTATCGCTCCTCGATATGTACGGATAGGCCACGGCATGTACGAGCATAAGCCGAACTGCACGAGCATACCCGACACTGATCTGAGAGTCGCTTAAATCATCTGCGTGTCTTAGAGCTTATCGCCATGTCTGTCAACTGTACTGTCATCCTGTAACTGTAGCGTATGTGAAAAAAAAAAAAAAAAAAAAAAAA";
//
//  cout<<"Query Size = "<<query.sequence.size()<<endl;
//  float discovery_fraction=0.7;
//  float  presence_fraction=0.0;
//
//  typedef std::vector<std::tuple<string, size_t, std::vector<size_t>>> LabelCountAbundancesVec;
//
//  LabelCountAbundancesVec result= anno_graph->get_kmer_counts(query.sequence,
//                                                    100,
//                                                    discovery_fraction,
//                                                    presence_fraction);
//  for(auto t: result)
//  {
//    cout<<"Label = "<<std::get<0>(t)<<endl;
//    cout<<"size = "<<std::get<1>(t)<<endl;
//    auto vec=std::get<2>(t);
//    cout<<"Vec Size = "<<vec.size();
//    unsigned nonZero=0;
//    for(unsigned i=0;i<vec.size();i++)
//    {
//      if(i%60 == 0)
//        cout<<endl;
//      cout<<vec[i]<<" ";
//      if(vec[i]!=0)
//        nonZero++;
//    }
//    cout<<endl;
//    cout<<"Non Zero = "<<nonZero<<endl;
//    cout<<endl;
//  }
//
//          return 0;
	Timer timer;
	double time_preprocessing;
	double time_kmer_counting;
	double time_unique_kmers;
	double time_path_sampling;
	double time_writing;
	double time_total;

	cerr << endl;
	cerr << "program: PanGenie - genotyping and phasing based on kmer-counting and known haplotype sequences." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	string readfile = "";
	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	string outname = "result";
	string sample_name = "sample";
	size_t nr_jellyfish_threads = 1;
	size_t nr_core_threads = 1;
	bool only_genotyping = true;
	bool only_phasing = false;
	long double effective_N = 0.00001L;
	long double regularization = 0.001L;
	bool count_only_graph = true;
	bool ignore_imputed = false;
	bool add_reference = true;
	size_t sampling_size = 0;
	uint64_t hash_size = 3000000000;

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("PanGenie [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format or Jellyfish database in jf format. NOTE: INPUT FASTA/Q FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders.");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
	argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
	argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
//	argument_parser.add_optional_argument('n', "0.00001", "effective population size");
	argument_parser.add_flag_argument('g', "run genotyping (Forward backward algorithm, default behaviour).");
	argument_parser.add_flag_argument('p', "run phasing (Viterbi algorithm). Experimental feature.");
//	argument_parser.add_optional_argument('m', "0.001", "regularization constant for copynumber probabilities");
	argument_parser.add_flag_argument('c', "count all read kmers instead of only those located in graph.");
	argument_parser.add_flag_argument('u', "output genotype ./. for variants not covered by any unique kmers.");
	argument_parser.add_flag_argument('d', "do not add reference as additional path.");
//	argument_parser.add_optional_argument('a', "0", "sample subsets of paths of this size.");
	argument_parser.add_optional_argument('e', "3000000000", "size of hash used by jellyfish.");

	try {
		argument_parser.parse(argc, argv);
	} catch (const runtime_error& e) {
		argument_parser.usage();
		cerr << e.what() << endl;
		return 1;
	} catch (const exception& e) {
		return 0;
	}
	readfile = argument_parser.get_argument('i');
	reffile = argument_parser.get_argument('r');
	vcffile = argument_parser.get_argument('v');
	kmersize = stoi(argument_parser.get_argument('k'));
	outname = argument_parser.get_argument('o');
	sample_name = argument_parser.get_argument('s');
	nr_jellyfish_threads = stoi(argument_parser.get_argument('j'));
	nr_core_threads = stoi(argument_parser.get_argument('t'));

	bool genotyping_flag = argument_parser.get_flag('g');
	bool phasing_flag = argument_parser.get_flag('p');
	
	if (genotyping_flag && phasing_flag) {
		only_genotyping = false;
		only_genotyping = false;
	}
	if (!genotyping_flag && phasing_flag) {
		only_genotyping = false;
		only_phasing = true;
	}

//	effective_N = stold(argument_parser.get_argument('n'));
//	regularization = stold(argument_parser.get_argument('m'));
	count_only_graph = !argument_parser.get_flag('c');
	ignore_imputed = argument_parser.get_flag('u');
	add_reference = !argument_parser.get_flag('d');
//	sampling_size = stoi(argument_parser.get_argument('a'));
	istringstream iss(argument_parser.get_argument('e'));
	iss >> hash_size;

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// check if input files exist and are uncompressed
	check_input_file(reffile);
	check_input_file(vcffile);
	check_input_file(readfile);

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	VariantReader variant_reader (vcffile, reffile, kmersize, add_reference, sample_name);
	
	// TODO: only for analysis
	struct rusage r_usage00;
	getrusage(RUSAGE_SELF, &r_usage00);
	cerr << "#### Memory usage until now: " << (r_usage00.ru_maxrss / 1E6) << " GB ####" << endl;
	
	string segment_file = outname + "_path_segments.fasta";
	cerr << "Write path segments to file: " << segment_file << " ..." << endl;
	variant_reader.write_path_segments(segment_file);

	// determine chromosomes present in VCF
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);
	cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

	// TODO: only for analysis
	struct rusage r_usage0;
	getrusage(RUSAGE_SELF, &r_usage0);
	cerr << "#### Memory usage until now: " << (r_usage0.ru_maxrss / 1E6) << " GB ####" << endl;

	time_preprocessing = timer.get_interval_time();

	// UniqueKmers for each chromosome
	UniqueKmersMap unique_kmers_list;
        map<string,ProbabilityTable>  probabilities;

	{
		KmerCounter* read_kmer_counts = nullptr;
		// determine kmer copynumbers in reads
		if (readfile.substr(std::max(3, (int) readfile.size())-3) == std::string(".jf")) {
			cerr << "Read pre-computed read kmer counts ..." << endl;
			jellyfish::mer_dna::k(kmersize);
			read_kmer_counts = new JellyfishReader(readfile, kmersize);
		} else {
			cerr << "Count kmers in reads ..." << endl;
			if (count_only_graph) {
				read_kmer_counts = new JellyfishCounter(readfile, segment_file, kmersize, nr_jellyfish_threads, hash_size);
			} else {
				read_kmer_counts = new JellyfishCounter(readfile, kmersize, nr_jellyfish_threads, hash_size);
			}
		}

		size_t kmer_abundance_peak = read_kmer_counts->computeHistogram(10000, count_only_graph, outname + "_histogram.histo");
                cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

                map<string,size_t> kmer_abundance_peaks_per_sample;
                kmer_abundance_peaks_per_sample["sample"]=kmer_abundance_peak;
		// count kmers in allele + reference sequence
		cerr << "Count kmers in genome ..." << endl;
		JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);

		// TODO: only for analysis
		struct rusage r_usage1;
		getrusage(RUSAGE_SELF, &r_usage1);
		cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;


		time_kmer_counting = timer.get_interval_time();

		cerr << "Determine unique kmers ..." << endl;
		// determine number of cores to use
		size_t available_threads_uk = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size());
		size_t nr_cores_uk = min(nr_core_threads, available_threads_uk);
		if (nr_cores_uk < nr_core_threads) {
			cerr << "Warning: using " << nr_cores_uk << " for determining unique kmers." << endl;
		}

		// precompute probabilities
		probabilities["sample"] = ProbabilityTable(kmer_abundance_peak / 4, kmer_abundance_peak*4, 2*kmer_abundance_peak, regularization);

		{
			// create thread pool with at most nr_chromosomes threads
			pangenie::ThreadPool threadPool (nr_cores_uk);
			for (auto chromosome : chromosomes) {
				VariantReader* variants = &variant_reader;
				UniqueKmersMap* result = &unique_kmers_list;
				KmerCounter* genomic_counts = &genomic_kmer_counts;
                                map<string,ProbabilityTable>* probs = &probabilities;
				function<void()> f_unique_kmers = bind(prepare_unique_kmers, chromosome, genomic_counts, read_kmer_counts, variants, probs, result, kmer_abundance_peaks_per_sample);
				threadPool.submit(f_unique_kmers);
			}
		}

		// TODO: only for analysis
		struct rusage r_usage2;
		getrusage(RUSAGE_SELF, &r_usage2);
		cerr << "#### Memory usage until now: " << (r_usage2.ru_maxrss / 1E6) << " GB ####" << endl;

		delete read_kmer_counts;
		read_kmer_counts = nullptr;
		time_unique_kmers = timer.get_interval_time();
	}

	// TODO: only for analysis
	struct rusage r_usage3;
	getrusage(RUSAGE_SELF, &r_usage3);
	cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;

	// prepare subsets of paths to run on
	unsigned short nr_paths = variant_reader.nr_of_paths();
	// TODO: for too large panels, print waring
	if (nr_paths > 200) cerr << "Warning: panel is large and PanGenie might take a long time genotyping. Try reducing the panel size prior to genotyping." << endl;
	// handle case when sampling_size is not set
	if (sampling_size == 0) {
		if (nr_paths > 25) {
			sampling_size = 14;
		} else {
			sampling_size = nr_paths;		
		}
	}

	PathSampler path_sampler(nr_paths);
	vector<vector<unsigned short>> subsets;
	path_sampler.partition_samples(subsets, sampling_size);

	for (auto s : subsets) {
		for (auto b : s) {
			cout << b << endl;
		}
		cout << "-----" << endl;
	}

	if (!only_phasing) cerr << "Sampled " << subsets.size() << " subset(s) of paths each of size " << sampling_size << " for genotyping." << endl;

	// for now, run phasing only once on largest set of paths that can still be handled.
	// in order to use all paths, an iterative stradegie should be considered
	vector<unsigned short> phasing_paths;
	unsigned short nr_phasing_paths = min((unsigned short) nr_paths, (unsigned short) 30);
	path_sampler.select_single_subset(phasing_paths, nr_phasing_paths);
	if (!only_genotyping) cerr << "Sampled " << phasing_paths.size() << " paths to be used for phasing." << endl;
	time_path_sampling = timer.get_interval_time();
	
	// TODO: only for analysis
	struct rusage r_usage30;
	getrusage(RUSAGE_SELF, &r_usage30);
	cerr << "#### Memory usage until now: " << (r_usage30.ru_maxrss / 1E6) << " GB ####" << endl;

	cerr << "Construct HMM and run core algorithm ..." << endl;

	// determine max number of available threads for genotyping (at most one thread per chromosome and subsample possible)
	size_t available_threads = min(thread::hardware_concurrency(), (unsigned int) chromosomes.size() * (unsigned int) subsets.size());
	if (nr_core_threads > available_threads) {
		cerr << "Warning: using " << available_threads << " for genotyping." << endl;
		nr_core_threads = available_threads;
	}

	// run genotyping
	Results results;
	{
		// create thread pool
                pangenie::ThreadPool threadPool (nr_core_threads);
                for(auto uniqPair: unique_kmers_list.unique_kmers) {
                  string sampleName= uniqPair.first;
                  for (auto chromosome : chromosomes) {
                    vector<UniqueKmers *> *unique_kmers =
                        &unique_kmers_list.unique_kmers[sampleName][chromosome];
                    ProbabilityTable *probs = &(probabilities.find(sampleName)->second);
                    Results *r = &results;
                    // if requested, run phasing first
                    if (!only_genotyping) {
                      vector<unsigned short> *only_paths = &phasing_paths;
                      function<void()> f_genotyping =
                          bind(run_genotyping, chromosome, unique_kmers, probs,
                               false, true, effective_N, only_paths, r);
                      threadPool.submit(f_genotyping);
                    }

                    if (!only_phasing) {
                      // if requested, run genotying
                      for (size_t s = 0; s < subsets.size(); ++s) {
                        vector<unsigned short> *only_paths = &subsets[s];
                        function<void()> f_genotyping = bind(
                            run_genotyping, chromosome, unique_kmers, probs,
                            true, false, effective_N, only_paths, r);
                        threadPool.submit(f_genotyping);
                      }
                    }
                  }
                }
	}

	timer.get_interval_time();

	// output VCF
        for(auto result_pair:results.result) {

          string sampleName=result_pair.first;
          auto result= result_pair.second;

          // prepare output files
          if (! only_phasing)
            variant_reader.open_genotyping_outfile(outname+ "_" +sampleName + "_genotyping.vcf");
          if (! only_genotyping)
            variant_reader.open_phasing_outfile(outname+ "_" +sampleName  + "_phasing.vcf");


          cerr << "Write results to VCF ..."<<sampleName << endl;
          if (!(only_genotyping && only_phasing))
            assert(result.size() == chromosomes.size());
          // write VCF
          for (auto it = result.begin(); it != result.end();
               ++it) {
            if (!only_phasing) {
              // output genotyping results

              variant_reader.write_genotypes_of(
                  it->first, it->second,
                  &unique_kmers_list.unique_kmers[sampleName][it->first],
                  ignore_imputed);
            }
            if (!only_genotyping) {
              // output phasing results
              variant_reader.write_phasing_of(
                  it->first, it->second,
                  &unique_kmers_list.unique_kmers[sampleName][it->first],
                  ignore_imputed);
            }
          }

          if (!only_phasing)
            variant_reader.close_genotyping_outfile();
          if (!only_genotyping)
            variant_reader.close_phasing_outfile();
        }
	time_writing = timer.get_interval_time();
	time_total = timer.get_total_time();

	cerr << endl << "###### Summary ######" << endl;
	// output times
	cerr << "time spent reading input files:\t" << time_preprocessing << " sec" << endl;
	cerr << "time spent counting kmers: \t" << time_kmer_counting << " sec" << endl;
	cerr << "time spent selecting paths: \t" << time_path_sampling << " sec" << endl;
	cerr << "time spent determining unique kmers: \t" << time_unique_kmers << " sec" << endl; 
	// output per chromosome time
	double time_hmm = time_writing;
	for (auto chromosome : chromosomes) {
		double time_chrom = results.runtimes[chromosome] + unique_kmers_list.runtimes[chromosome];
		cerr << "time spent genotyping chromosome " << chromosome << ":\t" << time_chrom << endl;
		time_hmm += time_chrom;
	}
	cerr << "total running time:\t" << time_preprocessing + time_kmer_counting + time_path_sampling + time_unique_kmers +  time_hmm + time_writing << " sec"<< endl;
	cerr << "total wallclock time: " << time_total  << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	// destroy UniqueKmers
        for(auto uniq: unique_kmers_list.unique_kmers) {
          for (auto it = uniq.second.begin();
               it != uniq.second.end(); ++it) {
            for (size_t i = 0; i < it->second.size(); ++i) {
              delete it->second[i];
              it->second[i] = nullptr;
            }
          }
        }

	return 0;
}