#ifndef COLUMNINDEXER_HPP
#define COLUMNINDEXER_HPP

#include <utility>
#include <vector>
#include <map>
#include "uniquekmers.hpp"

/** 
* Assign column index to pair of paths it represents.
**/

typedef std::pair<unsigned char, unsigned char> State;

class ColumnIndexer {
public:
	ColumnIndexer(const UniqueKmers* unique_kmers);
	/** get column index of a pair of paths **/
	size_t get_index_of(size_t path_id1, size_t path_id2) const;
	/** get size of the column **/
	size_t get_size() const;
	/** get number of paths **/
	size_t get_nr_paths() const;
	/** get alleles at index **/
	void get_alleles(size_t index, unsigned char& allele1, unsigned char& allele2) const;
	/** get paths at index **/
	void get_paths(size_t index, size_t& path_id1, size_t& path_id2) const;
	/** return iterator for paths,index pairs **/
	std::map<size_t,size_t>::const_iterator paths_begin() const;
	std::map<size_t,size_t>::const_iterator paths_end() const;
	/** check if column is covered by a path **/
	bool exists_path(size_t path_id) const;

private:
	size_t variant_id;
	size_t nr_paths;
	/** maps path to its index **/
	std::map<size_t,size_t> path_to_index;
	/** maps index to its path **/
	std::vector<size_t> index_to_path;
	/** maps index to allele **/
	std::vector<unsigned char> index_to_allele;
};

#endif // COLUMNINDEXER_HPP
