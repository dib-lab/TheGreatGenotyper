#include "kmercountcorrector.hpp"
#include <fstream>

using namespace std;

unsigned char encode (char base) {
	switch(base) {
		case 'A': return 0;
		case 'a': return 0;
		case 'C': return 1;
		case 'c': return 1;
		case 'G': return 2;
		case 'g': return 2;
		case 'T': return 3;
		case 't': return 3;
		default: return 4; 
	}
}

KmerCountCorrector::KmerCountCorrector(KmerCounter* genomic_kmers, FastaReader* fasta_reader, string& training_sequences, KmerCounter* read_kmers) {
	// read coordinates of training sequences from file
	ifstream infile(training_sequences);
	size_t start, end;
	string chromosome;
	while (infile >> chromosome >> start >> end) {
		// get corresponding sequence
		string result;
		fasta_reader.get_subsequence(chromosome, start, end, result);
		// enumerate all kmers
		// TODO get kmer_size
/**		size_t mask = (1 << (2*kmer_size)) - 1;
		size_t consecutive = 0;
		size_t current_kmer = 0;
		for (size_t i = 0; i < result.size(); ++i) {
			if (consecutive == kmer_size) {
				consecutive = 0;
				// TODO
				// check if kmer is unique
				// if so look up read count and store for regression
			} else {
				unsigned char code = encode(result.at(i));
				if (code > 3) {
					current_kmer = 0;
					consecutive = 0;
					continue;
				}
				consecutive += 1;
				current_kmer = ((current_kmer << 2) | code) & mask;
			}
			
		}
	}
**/	
}

size_t KmerCountCorrector::compute_corrected_count(string& kmer) {

}
