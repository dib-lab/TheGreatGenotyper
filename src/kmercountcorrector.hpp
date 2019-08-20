#ifndef KMERCOUNTCORRECTOR_HPP
#define KMERCOUNTCORRECTOR_HPP

#include <gsl/gsl_multifit.h>
#include <jellyfish/mer_dna.hpp>
#include <vector>
#include <string>
#include "kmercounter.hpp"
#include "fastareader.hpp"

class KmerCountCorrector {
public:
	KmerCountCorrector(KmerCounter* genomic_kmers, FastaReader* fasta_reader, std::string& training_sequences, KmerCounter* read_kmers, size_t kmer_size, size_t small_kmer_size);
	size_t compute_corrected_count(jellyfish::mer_dna& kmer);
private:
	size_t kmer_size;
	size_t small_kmer_size;
	std::vector<double> coefficients;
};

#endif // KMERCOUNTCORRECTOR_HPP
