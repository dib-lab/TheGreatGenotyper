#include "sequenceutils.hpp"
#include <stdexcept>

using namespace std;

unsigned char encode (char base) {
	switch (base) {
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

unsigned char complement (unsigned char base) {
	switch(base) {
		case 0: return 3;
		case 1: return 2;
		case 2: return 1;
		case 3: return 0;
		case 4: return 4;
		default: throw runtime_error("complement: invalid base.");
	}
}

char decode (unsigned char number) {
	switch (number) {
		case 0: return 'A';
		case 1: return 'C';
		case 2: return 'G';
		case 3: return 'T';
		default: return 'N';
	}
}

void split_in_kmers (jellyfish::mer_dna& kmer, size_t kmer_size, size_t small_kmer_size, map<unsigned int,int>& kmer_to_count) {
	jellyfish::mer_dna::k(kmer_size);
	size_t current_kmer = 0;
	size_t shifts = small_kmer_size;
	size_t mask = (1 << 2*small_kmer_size) - 1;
	for (size_t i = 0; i != kmer_size; ++i) {
		if (shifts == 0) {
			kmer_to_count[current_kmer] += 1;
		}
		char base = kmer.shift_left('A');
		current_kmer = ((current_kmer << 2) | encode(base)) & mask;
		if (shifts > 0) shifts -= 1;
	}
	kmer_to_count[current_kmer] += 1;
}
