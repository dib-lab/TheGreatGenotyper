#include "kmercountcorrector.hpp"
#include <fstream>
#include <math.h>

using namespace std;

unsigned char code (char base) {
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


void split_in_kmers (jellyfish::mer_dna& kmer, size_t kmer_size, size_t small_kmer_size, map<unsigned int, unsigned int>& kmer_to_count) {
	unsigned int current_kmer = 0;
	unsigned int shifts = small_kmer_size;
	unsigned int mask = (1 << (2*small_kmer_size)) - 1;
	for (size_t i = 0; i != kmer_size; ++i) {
		if (shifts == 0) {
			kmer_to_count[current_kmer] += 1;
		}
		current_kmer = ((current_kmer << 2) | code(kmer.shift_left('N'))) & mask;
		if (shifts > 0) shifts -= 1;
	}
	kmer_to_count[current_kmer] += 1;
}

KmerCountCorrector::KmerCountCorrector(KmerCounter* genomic_kmers, FastaReader* fasta_reader, string& training_sequences, KmerCounter* read_kmers, size_t kmer_size, size_t small_kmer_size)
	: kmer_size(kmer_size),
	  small_kmer_size(small_kmer_size)
{
	// read coordinates of training sequences from file
	ifstream infile(training_sequences);
	size_t start, end;
	string chromosome;
	vector < pair<jellyfish::mer_dna, size_t> > unique_kmers;
	while (infile >> chromosome >> start >> end) {
		// get corresponding sequence
		DnaSequence result;
		fasta_reader->get_subsequence(chromosome, start, end, result);
		// enumerate all kmers and find unique ones
		size_t extra_shifts = kmer_size;
		jellyfish::mer_dna::k(kmer_size);
		jellyfish::mer_dna current_kmer("");
		for (size_t i = 0; i < result.size(); ++i) {
			char current_base = result[i];
			if (extra_shifts == 0) {
				if (genomic_kmers->getKmerAbundance(current_kmer) == 1) {
					unique_kmers.push_back(make_pair(current_kmer, read_kmers->getKmerAbundance(current_kmer)));
				}
			}
			if (( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
				extra_shifts = kmer_size + 1;
			}
			current_kmer.shift_left(current_base);
			if (extra_shifts > 0) extra_shifts -= 1;
		}

		// add last one
		if (genomic_kmers->getKmerAbundance(current_kmer) == 1) {
			unique_kmers.push_back(make_pair(current_kmer, read_kmers->getKmerAbundance(current_kmer)));
		}

	}
	// construct input for regression
	size_t nr_small_kmers = pow(4, small_kmer_size);
	// initialize matrix and response
	gsl_matrix* X = gsl_matrix_calloc(unique_kmers.size(), nr_small_kmers);
	gsl_vector* y = gsl_vector_alloc(unique_kmers.size());
	gsl_vector* c = gsl_vector_alloc(nr_small_kmers);
	gsl_matrix* cov = gsl_matrix_alloc(nr_small_kmers, nr_small_kmers);
	// break each kmer into small kmers and fill matrix with counts
	size_t row_number = 0;
	for (auto it = unique_kmers.begin(); it != unique_kmers.end(); ++it) {
		map<unsigned int, unsigned int> kmer_to_count;
		split_in_kmers (it->first, kmer_size, small_kmer_size, kmer_to_count);
		for (auto count = kmer_to_count.begin(); count != kmer_to_count.end(); ++count) {
			gsl_matrix_set (X, row_number, count->first, count->second);
		}
		gsl_vector_set (y, row_number, it->second);
		row_number += 1;
	}
	// linear regression
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(unique_kmers.size(), nr_small_kmers);
	double chisq;
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);

	// store the coefficients
	for (size_t i = 0; i < nr_small_kmers; ++i) {
		this->coefficients.push_back( gsl_vector_get(c, (i)));
	}

	// clean up
	gsl_matrix_free (X);
	gsl_vector_free (y);
	gsl_vector_free (c);
	gsl_matrix_free (cov);
}

size_t KmerCountCorrector::compute_corrected_count(jellyfish::mer_dna& kmer) {
	size_t corrected_count = 0;
	// split given kmer into small kmers
	map<unsigned int, unsigned int> kmer_to_count;
	split_in_kmers (kmer, this->kmer_size, this->small_kmer_size, kmer_to_count);

	// compute corrected count
	for (auto it = kmer_to_count.begin(); it != kmer_to_count.end(); ++it) {
		corrected_count += it->second * this->coefficients.at(it->first);
	}
	return corrected_count;
}
