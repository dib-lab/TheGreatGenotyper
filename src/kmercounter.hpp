#ifndef KMERCOUNTER_HPP
#define KMERCOUNTER_HPP

#include <map>
#include <vector>
#include <string>
#include <jellyfish/mer_dna.hpp>
#include "fastareader.hpp"

class KmerCounter {
public:
	KmerCounter(size_t kmersize);

	/** get the abundance of given kmer (string) **/
	virtual size_t getKmerAbundance(std::string kmer) = 0;

	/** get the abundance of given kmer (jellyfish kmer) **/
	virtual size_t getKmerAbundance(jellyfish::mer_dna jelly_kmer) = 0;

	/** compute the kmer coverage relative to the number of kmers in the genome **/
	virtual size_t computeKmerCoverage(size_t genome_kmers) = 0;

	/** computes kmer abundance histogram and returns the three highest peaks **/
	virtual size_t computeHistogram(size_t max_count, std::string filename = "") = 0;

	/** compute corrected read kmer counts **/
	void correct_read_counts (KmerCounter* genomic_kmers, FastaReader* fasta_reader, std::string& training_sequences, size_t small_kmer_size, double train_frac);

	virtual ~KmerCounter () {};

protected:
	std::vector<double> coefficients;
	size_t kmer_size;
	size_t small_kmer_size;
	bool corrected;
	size_t compute_corrected_count(jellyfish::mer_dna& kmer) const;
};
#endif // KMERCOUNTER_HPP
