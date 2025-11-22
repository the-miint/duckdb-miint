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

	std::vector<std::string> exp_sample_ids = {"Sample1", "Sample1", "Sample2", "Sample2", "Sample2",
	                                           "Sample3", "Sample3", "Sample3", "Sample3", "Sample4",
	                                           "Sample4", "Sample5", "Sample6", "Sample6", "Sample6"};

	std::vector<std::string> exp_feature_ids = {"GG_OTU_2", "GG_OTU_4", "GG_OTU_2", "GG_OTU_4", "GG_OTU_5",
	                                            "GG_OTU_1", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5", "GG_OTU_2",
	                                            "GG_OTU_3", "GG_OTU_2", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4"};

	std::vector<double> exp_values = {5.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 4.0, 3.0, 1.0, 2.0, 1.0};

	// Convert indices to IDs for comparison
	const auto &sample_indices = table.COOSampleIndices();
	const auto &feature_indices = table.COOFeatureIndices();
	const auto &sample_ids_ordered = table.SampleIDs();
	const auto &feature_ids_ordered = table.FeatureIDs();

	std::vector<std::string> actual_sample_ids;
	std::vector<std::string> actual_feature_ids;
	for (size_t i = 0; i < sample_indices.size(); i++) {
		actual_sample_ids.push_back(sample_ids_ordered[sample_indices[i]]);
		actual_feature_ids.push_back(feature_ids_ordered[feature_indices[i]]);
	}

	REQUIRE((actual_sample_ids == exp_sample_ids));
	REQUIRE((actual_feature_ids == exp_feature_ids));
	REQUIRE((table.COOValues() == exp_values));
}

TEST_CASE("BIOM IsBIOM", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error printing
	auto path1 = "data/biom/test.biom";
	auto path2 = "data/biom/empty.biom";
	auto path3 = "data/sam/foo_no_header.sam";
	auto path4 = "data/biom/doesntexist";
	auto path5 = "data/biom/notbiom.h5";
	REQUIRE(miint::BIOMReader::IsBIOM(path1));
	REQUIRE(miint::BIOMReader::IsBIOM(path2));
	REQUIRE(!miint::BIOMReader::IsBIOM(path3));
	REQUIRE(!miint::BIOMReader::IsBIOM(path4));
}

TEST_CASE("BIOM file error recovery", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error printing
	// Try to open invalid file (should throw)
	REQUIRE_THROWS_AS(miint::BIOMReader("data/biom/notbiom.h5"), std::runtime_error);

	// Verify we can still open a valid file after the error
	REQUIRE_NOTHROW(miint::BIOMReader("data/biom/test.biom"));

	// And use it
	miint::BIOMReader reader("data/biom/test.biom");
	auto table = reader.read();
	REQUIRE((table.nnz() == 15));
}
