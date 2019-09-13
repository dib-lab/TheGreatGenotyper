#include "catch.hpp"
#include "utils.hpp"
#include "../src/jellyfishcounter.hpp"
#include "../src/sequenceutils.hpp"
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <jellyfish/mer_dna.hpp>

using namespace std;

TEST_CASE("Corrector", "[Corrector]") {
	string reads = "../tests/data/corrector.reads";
	string reference = "../tests/data/corrector.fasta";
	string train = "../tests/data/corrector.train";

	// read kmers
	JellyfishCounter read_kmers (reads, 7);
	// genomic kmers
	JellyfishCounter genomic_kmers (reference, 7);
	// fasta reader
	FastaReader fasta_reader (reference);
	// correct counts
//	read_kmers.correct_read_counts (&genomic_kmers, &fasta_reader, train, 2);
}

TEST_CASE("Corrector2", "[Corrector2]") {
	string reads = "../tests/data/corrector2.reads";
	string reference = "../tests/data/corrector2.fasta";
	string train = "../tests/data/corrector2.train";

	// read kmers
	JellyfishCounter read_kmers (reads, 16);
	// check uncorrected counts
	jellyfish::mer_dna::k(16);
	jellyfish::mer_dna current_kmer("AAAAAAAAAAAAAAAA");
	for (size_t i = 0; i < 17; ++i) {
		REQUIRE(read_kmers.getKmerAbundance(current_kmer) == 4);
		current_kmer.shift_left('C');
	}
	// genomic kmers
	JellyfishCounter genomic_kmers (reference, 16);
	// fasta reader
	FastaReader fasta_reader (reference);
	// correct counts
	read_kmers.correct_read_counts (&genomic_kmers, &fasta_reader, train, 2, 1.0);
	// check corrected counts
	current_kmer = jellyfish::mer_dna("AAAAAAAAAAAAAAAA");
	for (size_t i = 0; i < 17; ++i) {
		REQUIRE(read_kmers.getKmerAbundance(current_kmer) == 4);
		current_kmer.shift_left('C');
	}
}

TEST_CASE("split_in_kmers 1", "[split_in_kmers 1]") {
	string test = "ATGCGTGATCGATT";
	size_t small_kmer_size = 5;
	jellyfish::mer_dna::k(test.size());
	jellyfish::mer_dna kmer(test);
	map<unsigned int, int> kmer_to_count;
	split_in_kmers (kmer, test.size(), small_kmer_size, kmer_to_count);

	// kmers: "ATCGA", "ATGCG", "CGATT", "CGTGA", "GATCG", "GCGTG", "GTGAT", "TCGAT", "TGATC", "TGCGT"
	vector<unsigned int> expected_kmers = {216, 230, 399, 440, 566, 622, 739, 867, 909, 923};
	size_t index = 0;
	for (auto it = kmer_to_count.begin(); it != kmer_to_count.end(); ++it) {
		REQUIRE (it->first == expected_kmers[index]);
		REQUIRE (it->second == 1);
		index += 1;
	}
}

TEST_CASE("split_in_kmers 2", "[split_in_kmers 2]") {
	string test = "AAAAAAAA";
	size_t small_kmer_size = 5;
	jellyfish::mer_dna::k(8);
	jellyfish::mer_dna kmer (test);
	map<unsigned int, int> kmer_to_count;
	split_in_kmers (kmer, test.size(), small_kmer_size, kmer_to_count);

	// kmers: AAAAA
	unsigned int expected_kmer = 0;
	REQUIRE(kmer_to_count.size() == 1);
	REQUIRE(kmer_to_count[0] == 4);
}
