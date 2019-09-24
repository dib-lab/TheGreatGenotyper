#include <iostream>
#include <cmath>
#include <sstream>
#include "columnindexer.hpp"

using namespace std;

ColumnIndexer::ColumnIndexer(const UniqueKmers* unique_kmers)
	:nr_paths(unique_kmers->get_nr_paths()),
	 index_to_allele(this->nr_paths),
	 index_to_path(this->nr_paths)
{
	vector<size_t> path_ids;
	vector<unsigned char> alleles;
	unique_kmers->get_path_ids(path_ids, alleles);

	if (path_ids.size() == 0) {
		ostringstream oss;
		oss << "ColumnIndexer: column for variant " << unique_kmers->get_variant_position() << " is not covered by any paths." << endl;
		throw runtime_error(oss.str()); 
	}	
	for (size_t i = 0; i < this->nr_paths; ++i) {
		// fill path map
		this->path_to_index[path_ids[i]] = i;
		this->index_to_path[i] = path_ids[i];
		this->index_to_allele[i] = alleles[i];
	}
}

size_t ColumnIndexer::get_index_of(size_t path_id1, size_t path_id2) const {
	size_t index1 = this->path_to_index.at(path_id1);
	size_t index2 = this->path_to_index.at(path_id2);
	// compute index in column
	return index1*this->nr_paths + index2;
}

size_t ColumnIndexer::get_size() const {
	return this->nr_paths*this->nr_paths;
}

size_t ColumnIndexer::get_nr_paths() const {
	return this->nr_paths;
}

void ColumnIndexer::get_alleles(size_t index, unsigned char& allele1, unsigned char& allele2) const {
/**	if (index >= this->get_size()) {
		throw runtime_error("ColumnIndexer::get_alleles: index out of range.");
	}
**/
	// lookup alleles
	size_t id0 = index % this->nr_paths;
	size_t id1 = (index / this->nr_paths) % this->nr_paths;
	allele1 = this->index_to_allele[id1];
	allele2 = this->index_to_allele[id0];
}

void ColumnIndexer::get_paths(size_t index, size_t& path_id1, size_t& path_id2) const {
/**	if (index >= this->get_size()) {
		throw runtime_error("ColumnIndexer::get_paths: index out of range.");
	}
**/
	// lookup corresponding paths
	size_t id0 = index % this->nr_paths;
	size_t id1 = (index / this->nr_paths) % this->nr_paths;
	path_id1 = this->index_to_path[id1];
	path_id2 = this->index_to_path[id0];
}

map<size_t,size_t>::const_iterator ColumnIndexer::paths_begin() const {
	return this->path_to_index.begin();
}

map<size_t,size_t>::const_iterator ColumnIndexer::paths_end() const {
	return this->path_to_index.end();
}

bool ColumnIndexer::exists_path(size_t path_id) const {
	if (this->path_to_index.find(path_id) == this->path_to_index.end()) {
		return false;
	} else {
		return true;
	}
}
