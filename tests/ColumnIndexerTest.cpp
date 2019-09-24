#include "catch.hpp"
#include "../src/columnindexer.hpp"
#include <vector>
#include <string>

using namespace std;

TEST_CASE("ColumnIndexer", "[ColumnIndexer]") {
	UniqueKmers u (0,100);
	u.insert_path(0,0);
	u.insert_path(1,2);
	ColumnIndexer c(&u);
	REQUIRE(c.get_size() == 4);

	size_t path1, path2;
	c.get_paths(0, path1, path2);
	REQUIRE( ((path1 == 0) && (path2 == 0)) );

	c.get_paths(1, path1, path2);
	REQUIRE(((path1 == 0) && (path2 == 1)));

	c.get_paths(2, path1, path2);
	REQUIRE(((path1 == 1) && (path2 == 0)));
	
	c.get_paths(3, path1, path2);
	REQUIRE(((path1 == 1) && (path2 == 1)));

	unsigned char allele1, allele2;
	c.get_alleles(0, allele1, allele2);
	REQUIRE(((allele1 == 0) && (allele2 == 0)));

	c.get_alleles(1, allele1, allele2);
	REQUIRE(((allele1 == 0) && (allele2 == 2)));

	c.get_alleles(2, allele1, allele2);
	REQUIRE(((allele1 == 2) && (allele2 == 0)));

	c.get_alleles(3, allele1, allele2);
	REQUIRE(((allele1 == 2) && (allele2 == 2)));
}

TEST_CASE("ColumnIndexer test_invalid", "[ColumnIndexer test_invalid]") {
	UniqueKmers u (0,100);
	REQUIRE_THROWS(ColumnIndexer(&u));

	u.insert_path(3,0);
	ColumnIndexer c(&u);
	REQUIRE(c.get_size() == 1);

	size_t path1, path2;
	c.get_paths(0, path1, path2);
	REQUIRE(((path1 == 3) && (path2 == 3)));

	unsigned char allele1, allele2;
	c.get_alleles(0, allele1, allele2);
	REQUIRE(((allele1 == 0) && (allele2 == 0)));

	REQUIRE_THROWS(c.get_paths(1, path1, path2));
	REQUIRE_THROWS(c.get_alleles(1, allele1, allele2));
}
