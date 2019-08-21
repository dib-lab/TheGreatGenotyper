#include "catch.hpp"
#include "utils.hpp"
#include "../src/kmercounter.hpp"
#include <vector>
#include <string>
#include <iostream>

using namespace std;

TEST_CASE("Corrector", "[Corrector]") {
	string reads = "../tests/data/corrector.reads";
	string reference = "../tests/data/corrector.fasta";
	string train = "../tests/data/corrector.train";

	// read kmers
	KmerCounter read_kmers (reads, 7);
	// genomic kmers
	KmerCounter genomic_kmers (reference, 7);
	// fasta reader
	FastaReader fasta_reader (reference);
	// correct counts
	read_kmers.correct_read_counts (&genomic_kmers, &fasta_reader, train, 2);
}
