#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include "BIOMReader.hpp"

TEST_CASE("BIOM empty table", "[BIOMReader]") {
	auto path = "data/biom/empty.biom";
	miint::BIOMReader reader(path);
	auto table = reader.read();

	REQUIRE((table.COOFeatureIndices().size() == 0));
	REQUIRE((table.COOSampleIndices().size() == 0));
	REQUIRE((table.COOValues().size() == 0));
	REQUIRE((table.nnz() == 0));
}

TEST_CASE("BIOM basic table", "[BIOMReader]") {
	auto path = "data/biom/test.biom";
	miint::BIOMReader reader(path);
	auto table = reader.read();

	REQUIRE((table.COOFeatureIndices().size() == 15));
	REQUIRE((table.COOSampleIndices().size() == 15));
	REQUIRE((table.COOValues().size() == 15));

	REQUIRE((table.nnz() == 15));
}

