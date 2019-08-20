#ifndef KMERCOUNTCORRECTOR_HPP
#define KMERCOUNTCORRECTOR_HPP

#include<gsl/gsl_multifit.h>
#include "kmercounter.hpp"
#include "fastareader.hpp"

class KmerCountCorrector {
	KmerCountCorrector(KmerCounter* genomic_kmers, FastaReader* fasta_reader, std::string& training_sequences, KmerCounter* read_kmers);
	size_t compute_corrected_count(std::string& kmer);
};

#endif // KMERCOUNTCORRECTOR_HPP
