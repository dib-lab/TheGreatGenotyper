#include "kmercounter.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <fstream>
#include "histogram.hpp"
#include <gsl/gsl_multifit.h>
#include "sequenceutils.hpp"

using namespace std;

KmerCounter::KmerCounter(size_t kmersize) 
	:kmer_size(kmersize),
	corrected(false)
{}

double KmerCounter::compute_scaling_factor(jellyfish::mer_dna kmer) const {
	double scaling_factor = 0.0;
	// split given kmer into small kmers
	map<unsigned int, int> kmer_to_count;
	split_in_kmers (kmer, this->kmer_size, this->small_kmer_size, kmer_to_count);

	// compute scaling factor
	for (auto it = kmer_to_count.begin(); it != kmer_to_count.end(); ++it) {
		scaling_factor += it->second * this->coefficients.at(it->first);
	}
	return scaling_factor;
}

size_t KmerCounter::compute_corrected_count(jellyfish::mer_dna kmer, size_t raw_count) const {
	double scaling_factor = compute_scaling_factor(kmer);
	if (scaling_factor > 0) {
		return round(raw_count / (double) scaling_factor);
	} else {
		return raw_count;
	}
}


void KmerCounter::correct_read_counts (KmerCounter* genomic_kmers, FastaReader* fasta_reader, string& training_sequences, size_t small_kmer_size, double train_frac, double coverage) {
	this->small_kmer_size = small_kmer_size;
	if (this->corrected) return;
	// read coordinates of training sequences from file
	ifstream infile(training_sequences);
	size_t start, end;
	string chromosome;
	vector < pair<jellyfish::mer_dna, size_t> > unique_kmers;
	while (infile >> chromosome >> start >> end) {
		// decide whether to use interval
		double r_number = (double)rand() / (double)RAND_MAX;
		if (r_number >= train_frac) continue;
		// get corresponding sequence
		DnaSequence result;
		if (end < start) {
			cerr << "Skipping region: " << chromosome << " " << start << " " << end << endl;
			continue;
		}
		fasta_reader->get_subsequence(chromosome, start, end, result);
		// enumerate all kmers and find unique ones
		size_t extra_shifts = kmer_size;
		jellyfish::mer_dna::k(kmer_size);
		jellyfish::mer_dna current_kmer("");
		size_t kmers_added = 0;
		for (size_t i = 0; i < result.size(); ++i) {
			if (kmers_added > 100) break;
			char current_base = result[i];
			if (extra_shifts == 0) {
				if (genomic_kmers->getKmerAbundance(current_kmer) == 1) {
					unique_kmers.push_back(make_pair(current_kmer, this->getKmerAbundance(current_kmer)));
					// do not consider overlapping kmers for training
					extra_shifts = kmer_size;
					kmers_added += 1;
				}
			}
			if (( current_base != 'A') && (current_base != 'C') && (current_base != 'G') && (current_base != 'T') ) {
				extra_shifts = kmer_size + 1;
			}
			current_kmer.shift_left(current_base);
			if (extra_shifts > 0) extra_shifts -= 1;
		}

		// add last one
		if ((genomic_kmers->getKmerAbundance(current_kmer) == 1) && (extra_shifts == 0)) {
			unique_kmers.push_back(make_pair(current_kmer, this->getKmerAbundance(current_kmer)));
		}
		// TODO: size of training set?
//		if (unique_kmers.size() > 100*pow(4,small_kmer_size)) break;
		if (unique_kmers.size() > 200000) break;
	}
	infile.close();
	cerr << "Identified " << unique_kmers.size() << " unique kmers for training." << endl;
	// total number of possible kmers of size small_kmer_size
	size_t nr_small_kmers = pow(4, small_kmer_size);
	// initialize matrices and vectors needed for linear regression
	gsl_matrix* X = gsl_matrix_calloc(unique_kmers.size(), nr_small_kmers);
	gsl_vector* y = gsl_vector_alloc(unique_kmers.size());
	gsl_vector* c = gsl_vector_alloc(nr_small_kmers);
	gsl_matrix* cov = gsl_matrix_alloc(nr_small_kmers, nr_small_kmers);
	// break each kmer into small kmers and fill matrix with counts
	size_t row_number = 0;
	for (auto it = unique_kmers.begin(); it != unique_kmers.end(); ++it) {
		map<unsigned int,int> kmer_to_count;
		// TODO: remove
//		cout << it->first << endl;
		split_in_kmers (it->first, kmer_size, small_kmer_size, kmer_to_count);
		for (auto count = kmer_to_count.begin(); count != kmer_to_count.end(); ++count) {
			gsl_matrix_set (X, row_number, count->first, count->second);
		}
		gsl_vector_set (y, row_number, it->second / coverage);
		row_number += 1;
	}

	// TODO: remove
//	for (size_t i = 0; i < unique_kmers.size(); ++i) {
//		cout << unique_kmers[i].first << " ";
//		for (size_t j = 0; j < nr_small_kmers; ++j) {
//			cout << gsl_matrix_get(X, i, j) << " ";
//		}
//		cout << gsl_vector_get(y, i) << endl;
//	}

	// linear regression
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(unique_kmers.size(), nr_small_kmers);
	double chisq;
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);

	// TODO remove print
	cerr << "regression coefficients:" << endl;
	// store the coefficients
	for (size_t i = 0; i < nr_small_kmers; ++i) {
		cout << gsl_vector_get(c, (i)) << endl;
		this->coefficients.push_back( gsl_vector_get(c, (i)));
	}
	// clean up
	gsl_matrix_free (X);
	gsl_vector_free (y);
	gsl_vector_free (c);
	gsl_matrix_free (cov);

	// TODO: remove
//	const auto jf_ary = this->jellyfish_hash->ary();
//	const auto endle = jf_ary->end();
//	for (auto it = jf_ary->begin(); it != endle; ++it) {
//		auto& key_val = *it;
//		cout << compute_corrected_count(key_val.first) << " original: " << key_val.second << endl;
//	}
	this->corrected = true;
}
