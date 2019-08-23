#include <iostream>
#include <sstream>
#include <sys/resource.h>
#include "kmercounter.hpp"
#include "variantreader.hpp"
#include "kmercountcorrector.hpp"
#include "commandlineparser.hpp"

using namespace std;

int main (int argc, char* argv[])
{
	clock_t clock_start = clock();
	cerr << endl;
	cerr << "program: CountCorrector - compute corrected kmer counts using linear regression." << endl;
	cerr << "author: Jana Ebler" << endl << endl;
	string readfile = "";
	string reffile = "";
	string vcffile = "";
	size_t kmersize = 31;
	size_t small_kmersize = 5;
	string outname = "result";

	// parse the command line arguments
	CommandLineParser argument_parser;
	argument_parser.add_command("CountCorrector [options] -i <reads.fa/fq> -r <reference.fa> -v <variants.vcf>");
	argument_parser.add_mandatory_argument('i', "sequencing reads in FASTA/FASTQ format");
	argument_parser.add_mandatory_argument('r', "reference genome in FASTA format");
	argument_parser.add_mandatory_argument('v', "variants in VCF format");
	argument_parser.add_optional_argument('o', "result", "prefix of the output files");
	argument_parser.add_optional_argument('k', "31", "kmer size");
	argument_parser.add_optional_argument('s', "5", "small kmer size");
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
	small_kmersize = stoi(argument_parser.get_argument('s'));
	outname = argument_parser.get_argument('o');

	// print info
	cerr << "Files and parameters used:" << endl;
	argument_parser.info();

	// read allele sequences and unitigs inbetween, write them into file
	cerr << "Determine allele sequences ..." << endl;
	FastaReader reffile_reader(reffile);
	VariantReader variant_reader (vcffile, &reffile_reader, kmersize, "sample");
	string segment_file = outname + "_path_segments.fasta";
	cerr << "Write path segments to file: " << segment_file << " ..." << endl;
	variant_reader.write_path_segments(segment_file, true);

	// determine chromosomes present in VCF
	vector<string> chromosomes;
	variant_reader.get_chromosomes(&chromosomes);
	cerr << "Found " << chromosomes.size() << " chromosome(s) in the VCF." << endl;

	// TODO: only for analysis
	struct rusage r_usage0;
	getrusage(RUSAGE_SELF, &r_usage0);
	cerr << "#### Memory usage until now: " << (r_usage0.ru_maxrss / 1E6) << " GB ####" << endl;

	// determine kmer copynumbers in reads
	cerr << "Count kmers in reads ..." << endl;
	KmerCounter read_kmer_counts (readfile, kmersize);
	size_t kmer_abundance_peak = read_kmer_counts.computeHistogram(10000, outname + "_histogram.histo");
	cerr << "Computed kmer abundance peak: " << kmer_abundance_peak << endl;

	// count kmers in allele + reference sequence
	cerr << "Count kmers in genome ..." << endl;
	KmerCounter genomic_kmer_counts (segment_file, kmersize);

	cerr << "Compute corrected read kmer counts ..." << endl;
	// count correction
	string training_file = segment_file + ".train";
	read_kmer_counts.correct_read_counts(&genomic_kmer_counts, &reffile_reader, training_file, small_kmersize, 1/1000.0);

	size_t corrected_kmer_abundance_peak = read_kmer_counts.computeHistogram(10000, outname + "_corrected-histogram.histo");
	cerr << "Computed corrected kmer abundance peak: " << corrected_kmer_abundance_peak << endl;

	// TODO: only for analysis
	struct rusage r_usage1;
	getrusage(RUSAGE_SELF, &r_usage1);
	cerr << "#### Memory usage until now: " << (r_usage1.ru_maxrss / 1E6) << " GB ####" << endl;

	cerr << endl << "###### Summary ######" << endl;
	// total time
	double cpu_time = (double)(clock() - clock_start) / CLOCKS_PER_SEC;
	cerr << "Total CPU time: " << cpu_time << " sec" << endl;

	// memory usage
	struct rusage r_usage;
	getrusage(RUSAGE_SELF, &r_usage);
	cerr << "Total maximum memory usage: " << (r_usage.ru_maxrss / 1E6) << " GB" << endl;

	return 0;
}
