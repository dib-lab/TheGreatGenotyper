#include <iostream>
#include <fstream>
#include "fastareader.hpp"

using namespace std;

FastaReader::FastaReader(string filename, size_t kmer_size)
	:kmer_size(kmer_size)
{
	parse_file(filename);
}


// std::map<std::string, DnaSequence*> name_to_sequence;
void FastaReader::parse_file(string filename) {
	ifstream file(filename);
	string line;
	DnaSequence* dna_seq = nullptr;

	// non-N characters read since last N
	size_t current_len = 0;
	// was last character N?
	bool last_N = false;
	string name = "";
	while (getline(file, line)) {
		if (line.size() == 0) continue;
		size_t start = line.find_first_not_of(" \t\r\n");
		size_t end = line.find_last_not_of(" \t\r\n");
		line = line.substr(start, end-start+1);
		if (line.size() == 0) continue;
		// sequence name ?
		if (line[0] == '>') {
			start = line.find_first_not_of(" \t", 1);
			end = line.find_first_of(" \t", start);
			if (end == string::npos) end = line.size();

			// update kmer number of previous sequence
			if( (name != "") && (!last_N) && (current_len >= this->kmer_size) ) {
				this->name_to_kmers[name] += (current_len - this->kmer_size + 1);
			}
			current_len = 0;
			last_N = false;

			// get new sequence name
			name = line.substr(start,end-start);
			if (this->name_to_sequence.find(name) != this->name_to_sequence.end()) {
				// sequence with same name already seen, replace it
				delete this->name_to_sequence.at(name);
			}
			dna_seq = new DnaSequence;
			this->name_to_sequence[name] = dna_seq;
			this->name_to_kmers[name] = 0;

		} else {
			if ((dna_seq == nullptr) || (name == "")) {
				throw runtime_error("FastaReader::parse_file: file is malformatted.");
			} else {
				for (auto b : line) {
					if ( (b == 'N') or (b == 'n') ) {
						if ((!last_N) && (current_len >= this->kmer_size)) {
							this->name_to_kmers[name] += (current_len - this->kmer_size + 1);
						}
						last_N = true;
						current_len = 0;
					} else {
						current_len += 1;
						last_N = false;
					}
					string s(1,b);
					dna_seq->append(s);	
				}
			}
		}
	}

	if( (name != "") && (!last_N) && (current_len >= this->kmer_size) ) {
		this->name_to_kmers[name] += (current_len - this->kmer_size + 1);
	}
}

FastaReader::~FastaReader() {
	for (auto it = this->name_to_sequence.begin(); it != this->name_to_sequence.end(); ++it) {
		delete it->second;
	}
	this->name_to_sequence.clear();
}

bool FastaReader::contains_name(string name) const {
	auto it = this->name_to_sequence.find(name);
	return (it != this->name_to_sequence.end());
}

size_t FastaReader::get_size_of(string name) const {
	if (this->contains_name(name)) {
		return this->name_to_sequence.at(name)->length();
	} else {
		throw runtime_error("FastaReader::get_size_of: chromosome " + name + " is not present in FASTA-file.");
	}
}

size_t FastaReader::get_total_kmers() const {
	size_t total_kmers = 0;
	for (auto it = this->name_to_kmers.begin(); it != this->name_to_kmers.end(); ++it) {
		total_kmers += (it->second);
	}
	return total_kmers;
}

void FastaReader::get_subsequence(string name, size_t start, size_t end, string& result) const {
	if (this->contains_name(name)) {
		this->name_to_sequence.at(name)->substr(start, end, result);
	} else {
		throw runtime_error("FastaReader::get_subsequence: chromosome " + name + " is not present in FASTA-file.");
	}
}
