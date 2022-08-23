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
#include "SamplesDatabase.h"
#include "omp.h"
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
    vector<
    map<string, vector<UniqueKmers*> >
    > unique_kmers;
    map<string, double> runtimes;
};

struct Results {
    mutex result_mutex;
    vector<
    map<string, vector<GenotypingResult>>
    > result;
    map<string, double> runtimes;
};



void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, SamplesDatabase* database, VariantReader* variant_reader,UniqueKmersMap* unique_kmers_map) {
    Timer timer;
    UniqueKmerComputer kmer_computer(genomic_kmer_counts, database,variant_reader, chromosome);
    std::vector< std::vector<UniqueKmers*> > unique_kmers;
    unsigned numSamples= database->getNumSamples();
    // this needs to be map or vector of vectors for each sample
    kmer_computer.compute_unique_kmers(&unique_kmers);
    // store the results
    {
        lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
        for(unsigned sampleID=0; sampleID<numSamples ;sampleID++) {
            auto tmp = pair<string, vector<UniqueKmers *>>(
                    chromosome, move(unique_kmers[sampleID]));
            unique_kmers_map->unique_kmers[sampleID].insert(tmp);
        }
        unique_kmers_map->runtimes[chromosome]+=timer.get_total_time();
    }
}

void run_genotyping(string chromosome,unsigned  sampleID, vector<UniqueKmers*>* unique_kmers, ProbabilityTable* probs, bool only_genotyping, bool only_phasing, long double effective_N, vector<unsigned short>* only_paths, Results* results) {
    Timer timer;
    // construct HMM and run genotyping/phasing
    HMM hmm(unique_kmers, probs, !only_phasing, !only_genotyping, 1.26, false, effective_N, only_paths, false);
    // store the results
    {
        lock_guard<mutex> lock_result (results->result_mutex);

        // combine the new results to the already existing ones (if present)
        if (results->result[sampleID].find(chromosome) == results->result[sampleID].end()) {
            results->result[sampleID].insert(pair<string, vector<GenotypingResult>> (chromosome, hmm.move_genotyping_result()));
        } else {
            // combine newly computed likelihoods with already exisiting ones
            size_t index = 0;
            vector<GenotypingResult> genotypes = hmm.move_genotyping_result();
            for (auto likelihoods : genotypes) {
                results->result[sampleID].at(chromosome).at(index).combine(likelihoods);
                index += 1;
            }
        }
        // normalize the likelihoods after they have been combined
        for (size_t i = 0; i < results->result[sampleID].at(chromosome).size(); ++i) {
            results->result[sampleID].at(chromosome).at(i).normalize();
        }

        if (results->runtimes.find(chromosome) == results->runtimes.end()) {
            results->runtimes.insert(pair<string,double>(chromosome, timer.get_total_time()));
        } else {
            results->runtimes[chromosome] += timer.get_total_time();
        }
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
    string graphFile = "";
    string annotFile = "";
    string descriptionFile= "";
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
    argument_parser.add_mandatory_argument('i', "Metagraph graph path.dbg");
    argument_parser.add_mandatory_argument('a', "Metagraph annotations containig kmer counts");
    argument_parser.add_mandatory_argument('f', "Description file .tsv");

    argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
    argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
    argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders.");
    //argument_parser.add_optional_argument('k', "31", "kmer size");
//	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
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
    graphFile = argument_parser.get_argument('i');
    annotFile = argument_parser.get_argument('a');
    descriptionFile= argument_parser.get_argument('f');
    reffile = argument_parser.get_argument('r');
    vcffile = argument_parser.get_argument('v');
    //kmersize = stoi(argument_parser.get_argument('k'));
    outname = argument_parser.get_argument('o');
    sample_name = argument_parser.get_argument('s');
    nr_jellyfish_threads = stoi(argument_parser.get_argument('j'));
    nr_core_threads = stoi(argument_parser.get_argument('t'));

    omp_set_num_threads(nr_core_threads);
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
    check_input_file(descriptionFile);
    check_input_file(graphFile);
    check_input_file(annotFile);

    cerr << "Load Database ..."<< endl;
    SamplesDatabase database(graphFile,annotFile,descriptionFile,regularization);
    kmersize=database.getKSize();
    unsigned numSamples= database.getNumSamples();
    vector<string> sampleNames=database.getSamplesName();
    struct rusage r_usageD;
    getrusage(RUSAGE_SELF, &r_usageD);
    cerr << "#### Memory usage until now: " << (r_usageD.ru_maxrss / 1E6) << " GB ####" << endl;

    // read allele sequences and unitigs inbetween, write them into file
    cerr << "Determine allele sequences ..." << endl;
    VariantReader variant_reader (vcffile, reffile, kmersize, add_reference, sampleNames);
    // TODO: only for analysis
    struct rusage r_usage00;
    getrusage(RUSAGE_SELF, &r_usage00);
    cerr << "#### Memory usage until now: " << (r_usage00.ru_maxrss / 1E6) << " GB ####" << endl;

    string segment_file = outname + "_path_segments.fasta";
    cerr << "Write path segments to file: " << segment_file << " ..." << endl;
    variant_reader.write_path_segments(segment_file);

    // determine chromosomes present in VCF
    vector<string> chromosomes= variant_reader.get_chromosomes_vcfSorted();
    cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;


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
    struct rusage r_usage0;
    getrusage(RUSAGE_SELF, &r_usage0);
    cerr << "#### Memory usage until now: " << (r_usage0.ru_maxrss / 1E6) << " GB ####" << endl;

    time_preprocessing = timer.get_interval_time();


    // count kmers in allele + reference sequence
    cerr << "Count kmers in genome ..." << endl;
    JellyfishCounter genomic_kmer_counts (segment_file, kmersize, nr_jellyfish_threads, hash_size);

    // TODO: only for analysis
    struct rusage r_usage1;
    getrusage(RUSAGE_SELF, &r_usage1);
    cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;


    time_kmer_counting = timer.get_interval_time();




    // UniqueKmers for each chromosome
    UniqueKmersMap unique_kmers_list;
    unique_kmers_list.unique_kmers.resize(numSamples);
    Results results;
    results.result.resize(numSamples);
    vector<VariantReader> outputVCFs(numSamples);

    for(unsigned sampleID=0; sampleID<numSamples ;sampleID++) {
        outputVCFs[sampleID] = variant_reader;
        string sampleName = database.getSampleName(sampleID);
        outputVCFs[sampleID].setSampleName(sampleName);
        // prepare output files
        if (!only_phasing)
            outputVCFs[sampleID].open_genotyping_outfile(outname + "_" + sampleName +
                                                         "_genotyping.vcf");
        if (!only_genotyping)
            outputVCFs[sampleID].open_phasing_outfile(outname + "_" + sampleName +
                                                      "_phasing.vcf");
    }


    for(auto chrom : chromosomes)
    {
        cerr << "Determine unique kmers for chromosome: "<< chrom << endl;
        VariantReader* variants = &variant_reader;
        UniqueKmersMap* result = &unique_kmers_list;
        KmerCounter* genomic_counts = &genomic_kmer_counts;
        timer.get_interval_time();
        prepare_unique_kmers( chrom, genomic_counts, &database, variants, result);
        time_unique_kmers += timer.get_interval_time();
        cerr<< "Finished Determining unique kmers for chromosome: "<< chrom << endl;
        struct rusage r_usage3;
        getrusage(RUSAGE_SELF, &r_usage3);
        cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;
        size_t available_threads = min(thread::hardware_concurrency(), numSamples);
        if (nr_core_threads > available_threads) {
            cerr << "Warning: using " << available_threads << " for genotyping." << endl;
            nr_core_threads = available_threads;
        }

        cerr << "Construct HMM and run core algorithm for chromosome: "<< chrom << endl;
        {
            // create thread pool
            pangenie::ThreadPool threadPool (nr_core_threads);
            for(unsigned sampleID=0; sampleID<numSamples ;sampleID++)  {
                vector<UniqueKmers *> *unique_kmers =
                        &unique_kmers_list.unique_kmers[sampleID][chrom];
                ProbabilityTable *probs = database.getSampleProbability(sampleID);
                Results *r = &results;
                // if requested, run phasing first
                if (!only_genotyping) {
                    vector<unsigned short> *only_paths = &phasing_paths;
                    function<void()> f_genotyping =
                            bind(run_genotyping, chrom, sampleID,unique_kmers, probs,
                                 false, true, effective_N, only_paths, r);
                    threadPool.submit(f_genotyping);
                }

                if (!only_phasing) {
                    // if requested, run genotying
                    for (size_t s = 0; s < subsets.size(); ++s) {
                        vector<unsigned short> *only_paths = &subsets[s];
                        function<void()> f_genotyping = bind(
                                run_genotyping, chrom, sampleID,unique_kmers, probs,
                                true, false, effective_N, only_paths, r);
                        threadPool.submit(f_genotyping);
                    }
                }
            }
        }

        cerr<< "Finished genotyping for chromosome: "<< chrom << endl;

        cerr<< "writing results for chromosome: "<< chrom << endl;
        timer.get_interval_time();
        for (unsigned sampleID = 0; sampleID < numSamples; sampleID++) {
            // write VCF

            if (!only_phasing) {
                // output genotyping results
                outputVCFs[sampleID].write_genotypes_of(
                        chrom, results.result[sampleID][chrom],
                        &unique_kmers_list.unique_kmers[sampleID][chrom],
                        ignore_imputed);
            }
            if (!only_genotyping) {
                // output phasing results
                outputVCFs[sampleID].write_phasing_of(
                        chrom, results.result[sampleID][chrom],
                        &unique_kmers_list.unique_kmers[sampleID][chrom],
                        ignore_imputed);
            }
        }
        for(auto uniq: unique_kmers_list.unique_kmers) {
            for (size_t i = 0; i < uniq[chrom].size(); ++i) {
                delete uniq[chrom][i];
                uniq[chrom][i] = nullptr;
            }
            uniq[chrom].clear();
        }
        results.result.clear();

        time_writing += timer.get_interval_time();
        cerr<< "Finished chromosome: "<< chrom << endl;


        getrusage(RUSAGE_SELF, &r_usage3);
        cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;


    }

    for(unsigned sampleID=0; sampleID<numSamples ;sampleID++) {
        if (!only_phasing)
            outputVCFs[sampleID].close_genotyping_outfile();
        if (!only_genotyping)
            outputVCFs[sampleID].close_phasing_outfile();
    }


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

    return 0;
}
