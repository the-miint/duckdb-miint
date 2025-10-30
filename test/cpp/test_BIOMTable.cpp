#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include "BIOMTable.hpp"

TEST_CASE("BIOM from COO", "[BIOMTable]") {
	SECTION("Simple") {
		std::vector<std::string> features = {"F1", "F2", "F3"};
		std::vector<std::string> samples = {"S1", "S2", "S3"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		std::vector<size_t> exp_feature_indices = {0, 1, 2};
		std::vector<size_t> exp_sample_indices = {0, 1, 2};
		std::vector<double> exp_values = {1.0, 2.0, 3.0};
		REQUIRE((table.COOFeatureIndices() == exp_feature_indices));
		REQUIRE((table.COOSampleIndices() == exp_sample_indices));
		REQUIRE((table.COOValues() == exp_values));
		REQUIRE((table.nnz() == 3));
	}

	SECTION("Duplicates") {
		std::vector<std::string> features = {"F1", "F2", "F1"};
		std::vector<std::string> samples = {"S1", "S2", "S1"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		std::vector<size_t> exp_feature_indices = {0, 1};
		std::vector<size_t> exp_sample_indices = {0, 1};
		std::vector<double> exp_values = {4.0, 2.0};
		REQUIRE((table.COOFeatureIndices() == exp_feature_indices));
		REQUIRE((table.COOSampleIndices() == exp_sample_indices));
		REQUIRE((table.COOValues() == exp_values));
		REQUIRE((table.nnz() == 2));
	}

	SECTION("Complex with zeros") {
		std::vector<std::string> features = {"F1", "F2", "F1", "F3", "F1", "F2", "F3"};
		std::vector<std::string> samples = {"S1", "S2", "S1", "S1", "S2", "S2", "S1"};
		std::vector<double> values = {1.0, 2.0, 3.0, 4.0, 5.0, 0.0, 6.0};
		miint::BIOMTable table(features, samples, values);

		// note: compression sorts by sample then feature, so our resulting order
		// may differ from input, while preserving the index <-> name association
		std::vector<size_t> exp_feature_indices = {0, 2, 0, 1};
		std::vector<size_t> exp_sample_indices = {0, 0, 1, 1};
		std::vector<double> exp_values = {1.0 + 3.0, 4.0 + 6.0, 5.0, 2.0 + 0.0};

		REQUIRE((table.COOFeatureIndices() == exp_feature_indices));
		REQUIRE((table.COOSampleIndices() == exp_sample_indices));
		REQUIRE((table.COOValues() == exp_values));
		REQUIRE((table.nnz() == 4));
	}
}

TEST_CASE("Unique ID ordering", "[BIOMTable]") {
	SECTION("Empty") {
		std::vector<std::string> names;
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		REQUIRE((result.size() == 0));
	}

	SECTION("Single element") {
		std::vector<std::string> names = {"foo"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		REQUIRE((result.size() == 1));
		REQUIRE((result[0] == "foo"));
	}

	SECTION("All same element") {
		std::vector<std::string> names = {"foo", "foo", "foo", "foo"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		REQUIRE((result.size() == 1));
		REQUIRE((result[0] == "foo"));
	}

	SECTION("Already sorted") {
		std::vector<std::string> names = {"aaa", "bbb", "ccc", "ddd"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		REQUIRE((result == names));
	}

	SECTION("Reverse sorted") {
		std::vector<std::string> names = {"zzz", "yyy", "xxx"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		std::vector<std::string> exp = {"zzz", "yyy", "xxx"};
		REQUIRE((result == exp));
	}

	SECTION("Duplicates at boundaries") {
		std::vector<std::string> names = {"aaa", "aaa", "bbb", "ccc", "ccc"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		std::vector<std::string> exp = {"aaa", "bbb", "ccc"};
		REQUIRE((result == exp));
	}

	SECTION("No duplicates") {
		std::vector<std::string> names = {"foo", "bar", "baz"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		REQUIRE((result.size() == 3));
		REQUIRE((result == names));
	}

	SECTION("Complex") {
		std::vector<std::string> names = {"foo", "bar", "bar", "foo", "bar", "foo", "baz", "foo", "bar", "foo"};
		std::vector<std::string> result = miint::unique_ids_in_order(names);
		std::vector<std::string> exp = {"foo", "bar", "baz"};
		REQUIRE((result.size() == 3));
		REQUIRE((result == exp));
	}
}
TEST_CASE("Permutation with cycle following", "[BIOMTable]") {

	SECTION("Identity permutation") {
		std::vector<size_t> rows = {0, 1, 2, 3};
		std::vector<size_t> cols = {0, 1, 2, 3};
		std::vector<double> values = {1.0, 2.0, 3.0, 4.0};
		std::vector<size_t> indices = {0, 1, 2, 3};

		miint::apply_permutation(rows, cols, values, indices);

		REQUIRE((rows == std::vector<size_t> {0, 1, 2, 3}));
		REQUIRE((cols == std::vector<size_t> {0, 1, 2, 3}));
		REQUIRE((values == std::vector<double> {1.0, 2.0, 3.0, 4.0}));
	}

	SECTION("Simple swap") {
		std::vector<size_t> rows = {0, 1};
		std::vector<size_t> cols = {0, 1};
		std::vector<double> values = {1.0, 2.0};
		std::vector<size_t> indices = {1, 0};

		miint::apply_permutation(rows, cols, values, indices);

		REQUIRE((rows == std::vector<size_t> {1, 0}));
		REQUIRE((cols == std::vector<size_t> {1, 0}));
		REQUIRE((values == std::vector<double> {2.0, 1.0}));
	}

	SECTION("Cycle of three") {
		std::vector<size_t> rows = {0, 1, 2};
		std::vector<size_t> cols = {10, 11, 12};
		std::vector<double> values = {100.0, 200.0, 300.0};
		std::vector<size_t> indices = {1, 2, 0}; // 0->1, 1->2, 2->0

		miint::apply_permutation(rows, cols, values, indices);

		REQUIRE((rows == std::vector<size_t> {1, 2, 0}));
		REQUIRE((cols == std::vector<size_t> {11, 12, 10}));
		REQUIRE((values == std::vector<double> {200.0, 300.0, 100.0}));
	}

	SECTION("Multiple cycles") {
		std::vector<size_t> rows = {0, 1, 2, 3};
		std::vector<size_t> cols = {10, 11, 12, 13};
		std::vector<double> values = {100.0, 200.0, 300.0, 400.0};
		std::vector<size_t> indices = {1, 0, 3, 2}; // Two swaps: (0,1) and (2,3)

		miint::apply_permutation(rows, cols, values, indices);

		REQUIRE((rows == std::vector<size_t> {1, 0, 3, 2}));
		REQUIRE((cols == std::vector<size_t> {11, 10, 13, 12}));
		REQUIRE((values == std::vector<double> {200.0, 100.0, 400.0, 300.0}));
	}

	SECTION("Reverse order") {
		std::vector<size_t> rows = {0, 1, 2, 3, 4};
		std::vector<size_t> cols = {10, 11, 12, 13, 14};
		std::vector<double> values = {100.0, 200.0, 300.0, 400.0, 500.0};
		std::vector<size_t> indices = {4, 3, 2, 1, 0};

		miint::apply_permutation(rows, cols, values, indices);

		REQUIRE((rows == std::vector<size_t> {4, 3, 2, 1, 0}));
		REQUIRE((cols == std::vector<size_t> {14, 13, 12, 11, 10}));
		REQUIRE((values == std::vector<double> {500.0, 400.0, 300.0, 200.0, 100.0}));
	}
}

TEST_CASE("BIOMTable ToCSR conversion", "[BIOMTable]") {
	SECTION("Simple 2x3 matrix") {
		// Matrix (2 features × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  2.0  4.0
		//   F1    3.0  0.0  0.0
		std::vector<std::string> features = {"F0", "F0", "F0", "F1"};
		std::vector<std::string> samples = {"S0", "S1", "S2", "S0"};
		std::vector<double> values = {1.0, 2.0, 4.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		// Expected CSR (scipy verified):
		std::vector<double> exp_data = {1.0, 2.0, 4.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 1, 2, 0};
		std::vector<int32_t> exp_indptr = {0, 3, 4};

		REQUIRE((csr.data == exp_data));
		REQUIRE((csr.indices == exp_indices));
		REQUIRE((csr.indptr == exp_indptr));
	}

	SECTION("Empty matrix") {
		std::vector<std::string> features;
		std::vector<std::string> samples;
		std::vector<double> values;
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		REQUIRE(csr.data.empty());
		REQUIRE(csr.indices.empty());
		REQUIRE((csr.indptr == std::vector<int32_t> {0}));
	}

	SECTION("Single row") {
		// Matrix (1 feature × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  2.0  3.0
		std::vector<std::string> features = {"F0", "F0", "F0"};
		std::vector<std::string> samples = {"S0", "S1", "S2"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 1, 2};
		std::vector<int32_t> exp_indptr = {0, 3};

		REQUIRE((csr.data == exp_data));
		REQUIRE((csr.indices == exp_indices));
		REQUIRE((csr.indptr == exp_indptr));
	}

	SECTION("Single column") {
		// Matrix (3 features × 1 sample):
		//         S0
		//   F0    1.0
		//   F1    2.0
		//   F2    3.0
		std::vector<std::string> features = {"F0", "F1", "F2"};
		std::vector<std::string> samples = {"S0", "S0", "S0"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 0, 0};
		std::vector<int32_t> exp_indptr = {0, 1, 2, 3};

		REQUIRE((csr.data == exp_data));
		REQUIRE((csr.indices == exp_indices));
		REQUIRE((csr.indptr == exp_indptr));
	}

	SECTION("3x3 diagonal") {
		// Matrix (3 features × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  0.0  0.0
		//   F1    0.0  2.0  0.0
		//   F2    0.0  0.0  3.0
		std::vector<std::string> features = {"F0", "F1", "F2"};
		std::vector<std::string> samples = {"S0", "S1", "S2"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 1, 2};
		std::vector<int32_t> exp_indptr = {0, 1, 2, 3};

		REQUIRE((csr.data == exp_data));
		REQUIRE((csr.indices == exp_indices));
		REQUIRE((csr.indptr == exp_indptr));
	}

	SECTION("Single value") {
		std::vector<std::string> features = {"F0"};
		std::vector<std::string> samples = {"S0"};
		std::vector<double> values = {5.0};
		miint::BIOMTable table(features, samples, values);

		auto csr = table.ToCSR();

		std::vector<double> exp_data = {5.0};
		std::vector<int32_t> exp_indices = {0};
		std::vector<int32_t> exp_indptr = {0, 1};

		REQUIRE((csr.data == exp_data));
		REQUIRE((csr.indices == exp_indices));
		REQUIRE((csr.indptr == exp_indptr));
	}
}

TEST_CASE("BIOMTable ToCSC conversion", "[BIOMTable]") {
	SECTION("Simple 2x3 matrix") {
		// Matrix (2 features × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  2.0  4.0
		//   F1    3.0  0.0  0.0
		std::vector<std::string> features = {"F0", "F0", "F0", "F1"};
		std::vector<std::string> samples = {"S0", "S1", "S2", "S0"};
		std::vector<double> values = {1.0, 2.0, 4.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		// Expected CSC (scipy verified):
		std::vector<double> exp_data = {1.0, 3.0, 2.0, 4.0};
		std::vector<int32_t> exp_indices = {0, 1, 0, 0};
		std::vector<int32_t> exp_indptr = {0, 2, 3, 4};

		REQUIRE((csc.data == exp_data));
		REQUIRE((csc.indices == exp_indices));
		REQUIRE((csc.indptr == exp_indptr));
	}

	SECTION("Empty matrix") {
		std::vector<std::string> features;
		std::vector<std::string> samples;
		std::vector<double> values;
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		REQUIRE(csc.data.empty());
		REQUIRE(csc.indices.empty());
		REQUIRE((csc.indptr == std::vector<int32_t> {0}));
	}

	SECTION("Single row") {
		// Matrix (1 feature × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  2.0  3.0
		std::vector<std::string> features = {"F0", "F0", "F0"};
		std::vector<std::string> samples = {"S0", "S1", "S2"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 0, 0};
		std::vector<int32_t> exp_indptr = {0, 1, 2, 3};

		REQUIRE((csc.data == exp_data));
		REQUIRE((csc.indices == exp_indices));
		REQUIRE((csc.indptr == exp_indptr));
	}

	SECTION("Single column") {
		// Matrix (3 features × 1 sample):
		//         S0
		//   F0    1.0
		//   F1    2.0
		//   F2    3.0
		std::vector<std::string> features = {"F0", "F1", "F2"};
		std::vector<std::string> samples = {"S0", "S0", "S0"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 1, 2};
		std::vector<int32_t> exp_indptr = {0, 3};

		REQUIRE((csc.data == exp_data));
		REQUIRE((csc.indices == exp_indices));
		REQUIRE((csc.indptr == exp_indptr));
	}

	SECTION("3x3 diagonal") {
		// Matrix (3 features × 3 samples):
		//         S0   S1   S2
		//   F0    1.0  0.0  0.0
		//   F1    0.0  2.0  0.0
		//   F2    0.0  0.0  3.0
		std::vector<std::string> features = {"F0", "F1", "F2"};
		std::vector<std::string> samples = {"S0", "S1", "S2"};
		std::vector<double> values = {1.0, 2.0, 3.0};
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		std::vector<double> exp_data = {1.0, 2.0, 3.0};
		std::vector<int32_t> exp_indices = {0, 1, 2};
		std::vector<int32_t> exp_indptr = {0, 1, 2, 3};

		REQUIRE((csc.data == exp_data));
		REQUIRE((csc.indices == exp_indices));
		REQUIRE((csc.indptr == exp_indptr));
	}

	SECTION("Single value") {
		std::vector<std::string> features = {"F0"};
		std::vector<std::string> samples = {"S0"};
		std::vector<double> values = {5.0};
		miint::BIOMTable table(features, samples, values);

		auto csc = table.ToCSC();

		std::vector<double> exp_data = {5.0};
		std::vector<int32_t> exp_indices = {0};
		std::vector<int32_t> exp_indptr = {0, 1};

		REQUIRE((csc.data == exp_data));
		REQUIRE((csc.indices == exp_indices));
		REQUIRE((csc.indptr == exp_indptr));
	}
}
