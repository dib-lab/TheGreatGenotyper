#ifndef VARIANT_READER_HPP
#define VARIANT_READER_HPP

#include <string>
#include <map>
#include <vector>
#include <fstream>
#include <numeric>
#include <algorithm>
#include <cassert>
#include "fastareader.hpp"
#include "variant.hpp"
#include "genotypingresult.hpp"
#include "uniquekmers.hpp"

//std::vector<unsigned char> construct_index(std::vector<DnaSequence>& alleles, bool reference_added);
//std::vector<unsigned char> construct_index(std::vector<std::string>& alleles, bool reference_added);

template<class T>
std::vector<unsigned char> construct_index(std::vector<T>& alleles, bool reference_added) {
	size_t length = alleles.size();
	unsigned char offset = 0;
	if (reference_added) {
		assert(length > 0);
		length -= 1;
		offset += 1;
	}
	std::vector<unsigned char> index(length);
	std::iota(index.begin(), index.end(), 0);
	std::sort(index.begin(), index.end(), [&](unsigned char a, unsigned char b) { return alleles[a+offset] < alleles[b+offset]; });
	return index;
}

class VariantReader {
public:
	VariantReader (std::string filename, std::string reference_filename, size_t kmer_size, bool add_reference, std::vector<std::string> sample );
	/**  writes all path segments (allele sequences + reference sequences in between)
	*    to the given file.
	**/
    VariantReader (const VariantReader &old_obj);
    VariantReader (){}
    VariantReader& operator=(const VariantReader& other);
	size_t get_kmer_size() const;
	void write_path_segments(std::string filename) const;
	void get_chromosomes(std::vector<std::string>* result) const;
    std::vector<std::string> get_chromosomes_vcfSorted() const;
	size_t size_of(std::string chromosome) const;
	const Variant& get_variant(std::string chromosome, size_t index) const;
	const std::vector<Variant>& get_variants_on_chromosome(std::string chromosome) const;
	void open_genotyping_outfile(std::string outfile_name);
	void open_phasing_outfile(std::string outfile_name);
	void write_genotypes_of(std::string chromosome,std::string sample_name, const std::vector<GenotypingResult>& genotyping_result, bool ignore_imputed = false);
	void write_phasing_of(std::string chromosome, const std::vector<GenotypingResult>& genotyping_result, std::vector<UniqueKmers*>* unique_kmers, bool ignore_imputed = false);
	void close_genotyping_outfile();
	void close_phasing_outfile();
	size_t nr_of_genomic_kmers() const;
	size_t nr_of_paths() const;
	void get_left_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
	void get_right_overhang(std::string chromosome, size_t index, size_t length, DnaSequence& result) const;
    void setSampleName(std::string sampleName);
    void addVariantStat(unsigned int variantID, std::string sampleName,std::string chromosome, std::vector<VariantStats> & stat);
private:
	FastaReader fasta_reader;
	size_t kmer_size;
	size_t nr_paths;
	size_t nr_variants;
    std::vector<std::string> chromsomes;
	bool add_reference;
	std::vector<std::string> samples;
	std::vector<std::ofstream> genotyping_outfile;
	std::ofstream phasing_outfile;
	bool genotyping_outfile_open;
	bool phasing_outfile_open;
	static std::map< std::string, std::vector<Variant> > variants_per_chromosome;
	static std::map< std::string, std::vector<std::vector<std::string>>> variant_ids;
    std::map<std::string, std::map<std::string, std::vector<std::vector<VariantStats> >  > > variantsStatsPerSample;
	void add_variant_cluster(std::string& chromosome, std::vector<Variant>* cluster);
	void insert_ids(std::string& chromosome, std::vector<DnaSequence>& alleles, std::vector<std::string>& variant_ids, bool reference_added);
	std::string get_ids(std::string chromosome, std::vector<std::string>& alleles, size_t variant_index, bool reference_added);
};

#endif // VARIANT_READER_HPP
