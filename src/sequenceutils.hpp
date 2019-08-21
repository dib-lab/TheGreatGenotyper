#ifndef SEQUENCE_UTILS_HPP
#define SEQUENCE_UTILS_HPP

#include <jellyfish/mer_dna.hpp>
#include <map>

unsigned char encode (char base);

unsigned char complement (unsigned char base);

char decode (unsigned char number);

void split_in_kmers (jellyfish::mer_dna&, size_t kmer_size, size_t small_kmer_size, std::map<unsigned int,int>& kmer_to_count);

#endif // SEQUENCE_UTILS_HPP
