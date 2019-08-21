#include "kmercounter.hpp"
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <math.h>
#include <fstream>
#include "histogram.hpp"
#include <gsl/gsl_multifit.h>

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
		current_kmer = ((current_kmer << 2) | code(base)) & mask;
		if (shifts > 0) shifts -= 1;
	}
	kmer_to_count[current_kmer] += 1;
}

size_t KmerCounter::compute_corrected_count(jellyfish::mer_dna& kmer) const {
	double corrected_count = 0.0;
	// split given kmer into small kmers
	map<unsigned int, int> kmer_to_count;
	split_in_kmers (kmer, this->kmer_size, this->small_kmer_size, kmer_to_count);

	// compute corrected count
	for (auto it = kmer_to_count.begin(); it != kmer_to_count.end(); ++it) {
		corrected_count += it->second * this->coefficients.at(it->first);
	}
	if (corrected_count > 0) {
		return round(corrected_count);
	} else {
		return 0;
	}
}

KmerCounter::KmerCounter (string readfile, size_t kmer_size)
	:kmer_size(kmer_size),
	 corrected(false)
{
	jellyfish::mer_dna::k(kmer_size); // Set length of mers
	const uint64_t hash_size    = 10000000; // Initial size of hash.
	const uint32_t num_reprobes = 126;
	const uint32_t num_threads  = 1; // TODO Number of concurrent threads
	const uint32_t counter_len  = 7;  // Minimum length of counting field
	const bool canonical = true; // Use canonical representation

	// create the hash
	this->jellyfish_hash = new mer_hash_type(hash_size, jellyfish::mer_dna::k()*2, counter_len, num_threads, num_reprobes);

	// convert the string to char**
	vector<char*> args;
	istringstream iss (readfile);
	string token;
	while(iss >> token) {
		char* arg = new char [token.size() + 1];
		copy(token.begin(), token.end(), arg);
		arg[token.size()] = '\0';
		args.push_back(arg);
	}
	args.push_back(0);

	// count kmers
	mer_counter jellyfish_counter(num_threads, (*jellyfish_hash), &args[0], (&args[0])+1, canonical);
	jellyfish_counter.exec_join(num_threads);

	// delete the char**
	for(size_t i = 0; i < args.size(); i++)
		delete[] args[i];

}

size_t KmerCounter::getKmerAbundance(string kmer){

	jellyfish::mer_dna jelly_kmer(kmer);
	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	if (this->corrected) {
		val = compute_corrected_count(jelly_kmer);
	} else {
		jf_ary->get_val_for_key(jelly_kmer, &val);
	}

	const auto end = jf_ary->end();
	return val;
}

size_t KmerCounter::getKmerAbundance(jellyfish::mer_dna jelly_kmer){

	jelly_kmer.canonicalize();
	uint64_t val = 0;
	const auto jf_ary = this->jellyfish_hash->ary();
	if (this->corrected) {
		val = compute_corrected_count(jelly_kmer);
	} else {
		jf_ary->get_val_for_key(jelly_kmer, &val);
	}

	const auto end = jf_ary->end();
	return val;
}

size_t KmerCounter::computeKmerCoverage(size_t genome_kmers) const {
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto end = jf_ary->end();
	long double result = 0.0L;
	for (auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		long double count;
		if (this->corrected) {
			count = compute_corrected_count(key_val.first);
		} else {
			count = 1.0L * key_val.second;
		}
		long double genome = 1.0L * genome_kmers;
		result += (count/genome);
	}
	return (size_t) ceil(result);
}

size_t KmerCounter::computeHistogram(size_t max_count, string filename) const {
	Histogram histogram(max_count);
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto end = jf_ary->end();
	for (auto it = jf_ary->begin(); it != end; ++it) {
		auto& key_val = *it;
		size_t count;
		if (this->corrected) {
			count = compute_corrected_count(key_val.first);
		} else {
			count = key_val.second;
		}
		histogram.add_value(count);
	}
	// write histogram values to file
	if (filename != "") {
		histogram.write_to_file(filename);
	}
	// smooth the histogram
	histogram.smooth_histogram();
	// find peaks
	vector<size_t> peak_ids;
	vector<size_t> peak_values;
	histogram.find_peaks(peak_ids, peak_values);
	// identify the two largest peaks and return rightmost one
	if (peak_ids.size() < 2) {
		throw runtime_error("KmerCounter: less than 2 peaks found.");
	}
	size_t largest, second, largest_id, second_id;
	if (peak_values[0] < peak_values[1]){
		largest = peak_values[1];
		largest_id = peak_ids[1];
		second = peak_values[0];
		second_id = peak_ids[0];
	} else {
		largest = peak_values[0];
		largest_id = peak_ids[0];
		second = peak_values[1];
		second_id = peak_ids[1];
	}
	for (size_t i = 0; i < peak_values.size(); ++i) {
		if (peak_values[i] > largest) {
			second = largest;
			second_id = largest_id;
			largest = peak_values[i];
		} else if ((peak_values[i] > second) && (peak_values[i] != largest)) {
			second = peak_values[i];
			second_id = peak_ids[i];
		}
	}
	cerr << "Histogram peaks: " << largest_id << " (" << largest << "), " << second_id << " (" << second << ")" << endl;
	// add expected abundance counts to end of hist file
	if (filename != "") {
		ofstream histofile;
		histofile.open(filename, ios::app);
		histofile << "parameters\t" << 0.1 << '\t' << second_id/2.0 << '\t' << second_id << endl;
		histofile.close();
	}
	return second_id;
}

KmerCounter::~KmerCounter() {
	delete this->jellyfish_hash;
	this->jellyfish_hash = nullptr;
}

void KmerCounter::correct_read_counts (KmerCounter* genomic_kmers, FastaReader* fasta_reader, string& training_sequences, size_t small_kmer_size) {
	this->small_kmer_size = small_kmer_size;
	if (this->corrected) return;
	cout << "correct_read_counts: read input .." << endl;
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
					unique_kmers.push_back(make_pair(current_kmer, this->getKmerAbundance(current_kmer)));
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
			unique_kmers.push_back(make_pair(current_kmer, this->getKmerAbundance(current_kmer)));
		}
		// TODO remove, only for debugging
		if (unique_kmers.size() > 70) break;
	}
	infile.close();
	cout << "found " << unique_kmers.size() << " unique kmers for training." << endl;
	cout << "correct_read_counts: prepare regression" << endl;
	// construct input for regression
	size_t nr_small_kmers = pow(4, small_kmer_size);
	// initialize matrix and response
	gsl_matrix* X = gsl_matrix_calloc(unique_kmers.size(), nr_small_kmers);
	gsl_vector* y = gsl_vector_alloc(unique_kmers.size());
	gsl_vector* c = gsl_vector_alloc(nr_small_kmers);
	gsl_matrix* cov = gsl_matrix_alloc(nr_small_kmers, nr_small_kmers);
	// break each kmer into small kmers and fill matrix with counts
	size_t row_number = 0;
	cout << "before breaking" << endl;
	for (auto it = unique_kmers.begin(); it != unique_kmers.end(); ++it) {
		map<unsigned int,int> kmer_to_count;
		cout << it->first << endl;
		split_in_kmers (it->first, kmer_size, small_kmer_size, kmer_to_count);
		for (auto count = kmer_to_count.begin(); count != kmer_to_count.end(); ++count) {
			gsl_matrix_set (X, row_number, count->first, count->second);
		}
		gsl_vector_set (y, row_number, it->second);
		row_number += 1;
	}

	// TODO: remove
	for (size_t i = 0; i < unique_kmers.size(); ++i) {
		cout << unique_kmers[i].first << " ";
		for (size_t j = 0; j < nr_small_kmers; ++j) {
			cout << gsl_matrix_get(X, i, j) << " ";
		}
		cout << gsl_vector_get(y, i) << endl;
	}

	// linear regression
	gsl_multifit_linear_workspace * work = gsl_multifit_linear_alloc(unique_kmers.size(), nr_small_kmers);
	double chisq;
	gsl_multifit_linear (X, y, c, cov, &chisq, work);
	gsl_multifit_linear_free (work);
	cout << "store coefficients ... " << endl;

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
	const auto jf_ary = this->jellyfish_hash->ary();
	const auto endle = jf_ary->end();
	for (auto it = jf_ary->begin(); it != endle; ++it) {
		auto& key_val = *it;
		cout << compute_corrected_count(key_val.first) << " original: " << key_val.second << endl;
	}

	this->corrected = true;
}
