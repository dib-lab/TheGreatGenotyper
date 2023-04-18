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
#include "transitionprobabilitycomputer.hpp"
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
        ss << "File " << filename << " seems to be gzip-compressed. TheGreatGenotyper requires an uncompressed file." << endl;
        throw runtime_error(ss.str());
    }
}

struct UniqueKmersMap {
    mutex kmers_mutex;
    map<string, UniqueKmerComputer* > unique_kmers;
    map<string, double> runtimes;
};

struct Results {
    mutex result_mutex;
    map<string,
    map<string, vector<GenotypingResult>>
    > result;
    map<string, double> runtimes;
};



void prepare_unique_kmers(string chromosome, KmerCounter* genomic_kmer_counts, VariantReader* variant_reader,UniqueKmersMap* unique_kmers_map) {
    Timer timer;
    (*unique_kmers_map).unique_kmers[chromosome]=new UniqueKmerComputer(genomic_kmer_counts,variant_reader, chromosome);
    // this needs to be map or vector of vectors for each sample


    // store the results
    {
        lock_guard<mutex> lock_kmers (unique_kmers_map->kmers_mutex);
//        auto tmp = pair<string, vector<UniqueKmers *>>(
//        chromosome, move(unique_kmers));
        //unique_kmers_map->unique_kmers.insert(tmp);
        unique_kmers_map->runtimes[chromosome]+=timer.get_total_time();
    }
}

void run_genotyping(string chromosome,unsigned  sampleID,string sampleName, vector<UniqueKmers*>* unique_kmers,TransitionProbability* transitions, EmissionProbabilities* emissions, bool only_genotyping, bool only_phasing, vector<unsigned short>* only_paths, Results* results) {
    Timer timer;
    // construct HMM and run genotyping/phasing

    HMM hmm(unique_kmers,transitions,emissions, sampleID, !only_phasing, !only_genotyping, only_paths, false);

    // store the results
    {
        lock_guard<mutex> lock_result (results->result_mutex);
        if(results->result.find(chromosome) == results->result.end())
        {
            results->result[chromosome]=map<string, vector<GenotypingResult>>();
        }
        // combine the new results to the already existing ones (if present)
        if (results->result[chromosome].find(sampleName) == results->result[chromosome].end()) {
            results->result[chromosome].insert(pair<string, vector<GenotypingResult>> (sampleName, hmm.move_genotyping_result()));
        } else {
            // combine newly computed likelihoods with already exisiting ones
            size_t index = 0;
            vector<GenotypingResult> genotypes = hmm.move_genotyping_result();
            for (auto likelihoods : genotypes) {
                results->result[chromosome].at(sampleName).at(index).combine(likelihoods);
                index += 1;
            }
        }
        // normalize the likelihoods after they have been combined
        for (size_t i = 0; i < results->result[chromosome].at(sampleName).size(); ++i) {
            results->result[chromosome].at(sampleName).at(i).normalize();
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
    cerr << "program: The Great Genotyper for population genotyping by Moustafa Shokrof." << endl;
    cerr << "it is built on top of Pangenie tool by Jana Ebler" << endl << endl;
    string graphFolders = "";
    string graphFile = "";
    string annotFile = "";
    string descriptionFile= "";
    string reffile = "";
    string vcffile = "";
    string transitionsLoadFilePrefix = "";
    string transitionsSaveFilePrefix = "";
    string emissionsLoadFilePrefix = "";
    string emissionsSaveFilePrefix = "";
    string emissionsPrefix = "";
    size_t kmersize = 31;
    bool emissionOnly=false;
    string outname = "result";
    string sample_name = "sample";
    size_t nr_jellyfish_threads = 1;
    size_t nr_core_threads = 1;
    bool population_transitions = false;
    bool only_genotyping = true;
    bool only_phasing = false;
    long double effective_N = 0.00001L;
    long double regularization = 0.001L;
    bool count_only_graph = true;
    bool ignore_imputed = false;
    bool add_reference = true;
    size_t sampling_size = 0;
    uint64_t hash_size = 3000000000;
    bool log_scale=false;

    // parse the command line arguments
    CommandLineParser argument_parser;
    argument_parser.add_command("TheGreatGenotyper [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
    argument_parser.add_mandatory_argument('i', "Metagraph graph database folders");
   // argument_parser.add_mandatory_argument('a', "Metagraph annotations containig kmer counts");
   // argument_parser.add_mandatory_argument('f', "Description file .tsv");
    argument_parser.add_flag_argument('l', "the counts in the index are log scale.");

    argument_parser.add_mandatory_argument('r', "reference genome in FASTA format. NOTE: INPUT FASTA FILE MUST NOT BE COMPRESSED.");
    argument_parser.add_mandatory_argument('v', "variants in VCF format. NOTE: INPUT VCF FILE MUST NOT BE COMPRESSED.");
    argument_parser.add_optional_argument('o', "result", "prefix of the output files. NOTE: the given path must not include non-existent folders.");
    //argument_parser.add_optional_argument('k', "31", "kmer size");
//	argument_parser.add_optional_argument('s', "sample", "name of the sample (will be used in the output VCFs)");
    argument_parser.add_optional_argument('j', "1", "number of threads to use for kmer-counting");
    argument_parser.add_optional_argument('t', "1", "number of threads to use for core algorithm. Largest number of threads possible is the number of chromosomes given in the VCF");
    argument_parser.add_flag_argument('q', "use population transitions instead of LiStephens");
    argument_parser.add_flag_argument('a', "Genotype using kmers only");
    argument_parser.add_optional_argument('m', "", "Load the tranistions from this file");
    argument_parser.add_optional_argument('n', "", "compute the tranistions and save them to this file");

    argument_parser.add_optional_argument('x', "", "Load the emissions from this file");
    argument_parser.add_optional_argument('y', "", "compute the emissions and save them to this file");

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
    graphFolders = argument_parser.get_argument('i');
    //graphFile = argument_parser.get_argument('i');
    //annotFile = argument_parser.get_argument('a');
    //descriptionFile= argument_parser.get_argument('f');
    log_scale=argument_parser.get_flag('l');
    reffile = argument_parser.get_argument('r');
    vcffile = argument_parser.get_argument('v');
    //kmersize = stoi(argument_parser.get_argument('k'));
    outname = argument_parser.get_argument('o');
    sample_name = argument_parser.get_argument('s');
    nr_jellyfish_threads = stoi(argument_parser.get_argument('j'));
    nr_core_threads = stoi(argument_parser.get_argument('t'));
    emissionOnly = argument_parser.get_flag('a');
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

    population_transitions = argument_parser.get_flag('q');
    transitionsLoadFilePrefix= argument_parser.get_argument('m');
    transitionsSaveFilePrefix= argument_parser.get_argument('n');

    emissionsLoadFilePrefix= argument_parser.get_argument('x');
    emissionsSaveFilePrefix= argument_parser.get_argument('y');

    // print info
    cerr << "Files and parameters used:" << endl;
    argument_parser.info();

    // check if input files exist and are uncompressed
    check_input_file(reffile);
    check_input_file(vcffile);
    check_input_file(graphFolders);
//    check_input_file(descriptionFile);
//    check_input_file(graphFile);
//    check_input_file(annotFile);
    
    

    cerr << "Load Database ..."<< endl;
    vector<SamplesDatabase*> databases;
    ifstream databasesFolder(graphFolders);
    string prefix;
    unsigned numSamples= 0;
    vector<string> sampleNames;
    while(databasesFolder>>prefix)
    {
        graphFile = prefix + "graph.dbg";
        descriptionFile = prefix + "graph.desc.tsv";
        annotFile=prefix+"annotation.relaxed.row_diff_int_brwt.annodbg";
        databases.push_back(new SamplesDatabase(graphFile,annotFile,descriptionFile,regularization,log_scale));
        numSamples+= databases.back()->getNumSamples();
        vector<string> tmp=databases.back()->getSamplesName();
        for(auto s :tmp)
            sampleNames.push_back(s);
    }

    databases[0]->load_graph();
    kmersize=databases[0]->getKSize();


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
    if (nr_paths > 200) cerr << "Warning: panel is large and TheGreatGenotyper might take a long time genotyping. Try reducing the panel size prior to genotyping." << endl;
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
            cerr << b << endl;
        }
        cerr << "-----" << endl;
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
    JellyfishCounter* genomic_kmer_counts= new  JellyfishCounter(segment_file, kmersize, nr_jellyfish_threads, hash_size);

    // TODO: only for analysis
    struct rusage r_usage1;
    getrusage(RUSAGE_SELF, &r_usage1);
    cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;


    time_kmer_counting = timer.get_interval_time();




    // UniqueKmers for each chromosome
    UniqueKmersMap unique_kmers_list;

    for(auto chrom : chromosomes) {
        cerr << "Determine unique kmers for chromosome: " << chrom << endl;
        prepare_unique_kmers(chrom, genomic_kmer_counts, &variant_reader,&unique_kmers_list);
    }


    getrusage(RUSAGE_SELF, &r_usage1);
    cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

    //unique_kmers_list.unique_kmers.resize(numSamples);
    Results results;
    //results.result.resize(numSamples);
    variant_reader.open_genotyping_outfile(outname);
    map<string,vector<EmissionProbabilities*> > allEmissions;
    for(auto chrom : chromosomes)
    {
        allEmissions[chrom]=vector<EmissionProbabilities*>();
    }
    struct rusage r_usage3;
    cerr << "Calculate emissions  " << endl;
    if(emissionsSaveFilePrefix != "" && emissionsLoadFilePrefix== "")
    {
        emissionsPrefix=emissionsSaveFilePrefix;
    }
    else if(emissionsSaveFilePrefix == "" && emissionsLoadFilePrefix != "")
    {
        emissionsPrefix=emissionsLoadFilePrefix;
    }
    else{
        cerr<<"Exactly one of emissionsLoadFilePrefix and emissionsSaveFilePrefix has to be defined"<<endl;
        cerr<<"emissionsLoadFilePrefix(-x): "<<emissionsLoadFilePrefix<<endl;
        cerr<<"emissionsSaveFilePrefix(-y): "<<emissionsSaveFilePrefix<<endl;
        return -1;
    }


    for(unsigned i=0; i< databases.size(); i++)
    {
        cerr<< "Loading Database "<<i<<endl;
        if(i!=0)
            databases[i]->load_graph();
        for(auto chrom : chromosomes)
        {
            cerr << "Calculating emissions for chromosome: "<< chrom << endl;
            EmissionProbabilities* emissions=new EmissionProbabilities(databases[i],variant_reader.size_of(chrom));
            timer.get_interval_time();
            if(emissionsSaveFilePrefix != "") {
                unique_kmers_list.unique_kmers[chrom]->compute_emissions(databases[i], emissions);
                string filename=emissionsSaveFilePrefix+"."+chrom+ "."+to_string(i);
                emissions->save(filename);
                emissions->destroy();
            }
            allEmissions[chrom].push_back(emissions);
            time_unique_kmers += timer.get_interval_time();
            cerr<< "Finished Calculating emissions for chromosome: "<< chrom << endl;

            getrusage(RUSAGE_SELF, &r_usage3);
            cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;
        }
        databases[i]->delete_graph();
    }
    cerr <<"Finished computing emissions"<<endl;

    if(emissionOnly) {
        for (auto chrom: chromosomes) {
            map<string, vector<GenotypingResult>> res;
            for (unsigned i = 0; i < databases.size(); i++) {
                string filename=emissionsSaveFilePrefix+"."+chrom+ "."+to_string(i);
                allEmissions[chrom][i]->load(filename);
                allEmissions[chrom][i]->compute_most_likely_genotypes(&unique_kmers_list.unique_kmers[chrom]->uniqKmers);
                for (unsigned sampleID = 0; sampleID < databases[i]->getNumSamples(); sampleID++) {
                    string sampleName = databases[i]->getSampleName(sampleID);
                    res[sampleName]=allEmissions[chrom][i]->result[sampleID];
                }
                delete allEmissions[chrom][i];
            }
            variant_reader.write_genotypes_of(
                    chrom, res,
                    ignore_imputed);
        }
        time_total = timer.get_total_time();
        cerr<<"Finished GT using emissions only"<<endl;
        cerr<<"Total Time = "<< time_total /60.0 << " Minutes" << endl;
        return 0;
    }

    delete genomic_kmer_counts;

    for(auto chrom : chromosomes)
    {
        for(unsigned i=0; i< databases.size(); i++)
        {
            string filename=emissionsSaveFilePrefix+"."+chrom+ "."+to_string(i);
            allEmissions[chrom][i]->load(filename);
        }
        TransitionProbability* transitions;
        if(transitionsLoadFilePrefix != "")
        {
            cerr<<"Loading tranitions from "<<transitionsLoadFilePrefix+"."+chrom<<endl;
            transitions = new LiStephens(&variant_reader,chrom,1.26,effective_N);
//            vector<EmissionProbabilities*> tmp_emissions(4);
//            for(unsigned i=0;i <4 ;i++)
//            {
//                tmp_emissions[i]=new EmissionProbabilities();
//                string path=transitionsLoadFilePrefix +"."+chrom+"."+to_string(i);
//                tmp_emissions[i]->load(path);
//                tmp_emissions[i]->compute_most_likely_genotypes(&unique_kmers_list.unique_kmers[chrom]->uniqKmers);
//            }
//            transitions= new populationJointProbability(&variant_reader,chrom,tmp_emissions,&unique_kmers_list.unique_kmers[chrom]->uniqKmers);

//            for(unsigned i=0;i <4 ;i++)
//            {
//                delete tmp_emissions[i];
//            }
        }
        else {
            cerr<< "Calculating Transition probabilities for "<<chrom<<endl;
           // transitions= new populationJointProbability(&variant_reader,chrom,allEmissions[chrom],&(unique_kmers_list.unique_kmers[chrom]->uniqKmers));
            transitions= new LiStephens(&variant_reader,chrom,1.26,effective_N);
        }
        if(transitionsSaveFilePrefix != "")
        {
            cerr<< "Saving Transition probabilities to "<<transitionsSaveFilePrefix+"."+chrom<<endl;
            transitions->save(transitionsSaveFilePrefix+"."+chrom);
        }
        getrusage(RUSAGE_SELF, &r_usage3);
        cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;

        size_t available_threads = min(thread::hardware_concurrency(), numSamples);
        if (nr_core_threads > available_threads) {
            cerr << "Warning: using " << available_threads << " for genotyping." << endl;
            nr_core_threads = available_threads;
        }

        for(unsigned i=0; i< databases.size(); i++) {
            cerr<<"Processing database "<<i<<endl;
            cerr << "Construct HMM and run core algorithm for chromosome: " << chrom << endl;
            {
                // create thread pool
                pangenie::ThreadPool threadPool(nr_core_threads);
                for (unsigned sampleID = 0; sampleID < databases[i]->getNumSamples(); sampleID++) {
                    vector<UniqueKmers *> *unique_kmers =
                            &unique_kmers_list.unique_kmers[chrom]->uniqKmers;
                    Results *r = &results;
                    string sampleName= databases[i]->getSampleName(sampleID);
                    for (size_t s = 0; s < subsets.size(); ++s) {
                        vector<unsigned short> *only_paths = &subsets[s];
                        function<void()> f_genotyping = bind(
                                run_genotyping, chrom, sampleID,sampleName, unique_kmers, transitions, allEmissions[chrom][i],
                                true, false, only_paths, r);
                        threadPool.submit(f_genotyping);
                    }


                }
            }

            cerr << "Finished genotyping for chromosome: " << chrom << endl;

            cerr << "writing results for chromosome: " << chrom << endl;
            timer.get_interval_time();
            for (unsigned sampleID = 0; sampleID < databases[i]->getNumSamples(); sampleID++) {
                // write VCF
                // output phasing results
                string sampleName= databases[i]->getSampleName(sampleID);
                variant_reader.write_genotypes_of(
                        chrom, results.result[chrom],
                        ignore_imputed);

                results.result[sampleName][chrom].clear();
            }


        }
        for(auto uniq: unique_kmers_list.unique_kmers[chrom]->uniqKmers) {
            delete uniq;
        }
        delete transitions;
        for(unsigned i=0; i< databases.size(); i++)
        {
            delete allEmissions[chrom][i];
        }




        time_writing += timer.get_interval_time();
        cerr<< "Finished chromosome: "<< chrom << endl;


        getrusage(RUSAGE_SELF, &r_usage3);
        cerr << "#### Memory usage until now: " << (r_usage3.ru_maxrss / 1E6) << " GB ####" << endl;


    }

    variant_reader.close_genotyping_outfile();

    for(auto d: databases)
        delete d;

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
