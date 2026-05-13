#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "unifrac_metadata.hpp"

#include <stdexcept>
#include <string>
#include <vector>

using miint::unifrac::BuildGroupings;
using miint::unifrac::MetadataRow;
using miint::unifrac::NamedGrouping;

// PERMANOVA requires each sample in the distance matrix to have a categorical
// value for the variable being tested. Tests below pin the contract:
//   1) value-driven labeling so PERMANOVA equality across identical groupings
//      holds (Rule-7 invariant in the plan)
//   2) order-stable factorization (first-appearance in canonical sample order)
//      so seeded reproducibility survives metadata shuffling
//   3) loud failure on missing rows (silent NaN would corrupt the F-stat)
//   4) silent skip of extra rows (metadata may oversample relative to the
//      distance-matrix sample set)

TEST_CASE("BuildGroupings throws naming sample and variable on missing row", "[unifrac][metadata]") {
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"}, {"S2", "treatment", "trt"},
	    // S3 has no value for "treatment"
	};
	std::vector<std::string> samples = {"S1", "S2", "S3"};

	REQUIRE_THROWS_WITH(BuildGroupings(rows, samples, {"treatment"}),
	                    Catch::Matchers::ContainsSubstring("S3") && Catch::Matchers::ContainsSubstring("treatment"));
}

TEST_CASE("BuildGroupings single-group variable yields all-zero labels", "[unifrac][metadata]") {
	// PERMANOVA on a single-group variable produces degenerate stats but the
	// call must still succeed. n_groups=1 documents that to downstream.
	std::vector<MetadataRow> rows = {
	    {"S1", "site", "A"},
	    {"S2", "site", "A"},
	    {"S3", "site", "A"},
	};
	std::vector<std::string> samples = {"S1", "S2", "S3"};

	auto groupings = BuildGroupings(rows, samples, {"site"});
	REQUIRE(groupings.size() == 1);
	REQUIRE(groupings[0].variable == "site");
	REQUIRE(groupings[0].n_groups == 1);
	REQUIRE(groupings[0].labels == std::vector<uint32_t> {0, 0, 0});
}

TEST_CASE("BuildGroupings two variables with identical partition yield identical labels", "[unifrac][metadata]") {
	// WHY (Rule-7): factorization is *value-driven, not order-driven*. If two
	// variables partition the samples identically (just with different value
	// names), their grouping arrays must be byte-equal — PERMANOVA's F-stat
	// is a function of the partition, not the label strings.
	std::vector<MetadataRow> rows = {
	    {"S1", "v1", "alpha"}, {"S2", "v1", "beta"},  {"S3", "v1", "alpha"},
	    {"S1", "v2", "north"}, {"S2", "v2", "south"}, {"S3", "v2", "north"},
	};
	std::vector<std::string> samples = {"S1", "S2", "S3"};

	auto groupings = BuildGroupings(rows, samples, {"v1", "v2"});
	REQUIRE(groupings.size() == 2);
	REQUIRE(groupings[0].labels == groupings[1].labels);
	REQUIRE(groupings[0].n_groups == groupings[1].n_groups);
}

TEST_CASE("BuildGroupings is stable under input row shuffling", "[unifrac][metadata]") {
	// Labels follow first-appearance in canonical sample order (the
	// distance-matrix order), NOT input row order. Otherwise the seeded
	// PERMANOVA F/p reproducibility test breaks across re-runs of the same
	// metadata table in different orders.
	std::vector<MetadataRow> base = {
	    {"S1", "treatment", "ctrl"},
	    {"S2", "treatment", "trt"},
	    {"S3", "treatment", "ctrl"},
	    {"S4", "treatment", "trt"},
	};
	std::vector<MetadataRow> shuffled = {
	    {"S4", "treatment", "trt"},
	    {"S2", "treatment", "trt"},
	    {"S3", "treatment", "ctrl"},
	    {"S1", "treatment", "ctrl"},
	};
	std::vector<std::string> samples = {"S1", "S2", "S3", "S4"};

	auto a = BuildGroupings(base, samples, {"treatment"});
	auto b = BuildGroupings(shuffled, samples, {"treatment"});
	REQUIRE(a[0].labels == b[0].labels);
	// First-appearance in sample order: S1=ctrl -> 0, S2=trt -> 1.
	REQUIRE(a[0].labels == std::vector<uint32_t> {0, 1, 0, 1});
	REQUIRE(a[0].n_groups == 2);
}

TEST_CASE("BuildGroupings throws naming the unknown variable", "[unifrac][metadata]") {
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"},
	    {"S2", "treatment", "trt"},
	};
	std::vector<std::string> samples = {"S1", "S2"};

	REQUIRE_THROWS_WITH(BuildGroupings(rows, samples, {"site"}), Catch::Matchers::ContainsSubstring("site"));
}

TEST_CASE("BuildGroupings silently ignores metadata samples not in the distance matrix", "[unifrac][metadata]") {
	// Upstream convention: metadata tables often cover more samples than the
	// current analysis. Trying to error on "extra" metadata rows would force
	// users to pre-filter; libssu/skbb don't, so neither do we.
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"}, {"S2", "treatment", "trt"}, {"S_EXTRA", "treatment", "trt"}, // not in sample set
	};
	std::vector<std::string> samples = {"S1", "S2"};

	auto groupings = BuildGroupings(rows, samples, {"treatment"});
	REQUIRE(groupings.size() == 1);
	REQUIRE(groupings[0].labels.size() == 2);
	REQUIRE(groupings[0].n_groups == 2);
}

TEST_CASE("BuildGroupings with empty variables list tests all distinct variables in metadata", "[unifrac][metadata]") {
	// User-facing default: don't make the caller enumerate variables when they
	// want them all.
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"},
	    {"S2", "treatment", "trt"},
	    {"S1", "site", "A"},
	    {"S2", "site", "B"},
	};
	std::vector<std::string> samples = {"S1", "S2"};

	auto groupings = BuildGroupings(rows, samples, {});
	REQUIRE(groupings.size() == 2);

	// Order of returned variables: sorted lexicographically for stability.
	REQUIRE(groupings[0].variable == "site");
	REQUIRE(groupings[1].variable == "treatment");
}

TEST_CASE("BuildGroupings throws on empty sample list", "[unifrac][metadata]") {
	// Caller bug guard: no samples means no distance matrix, no PERMANOVA.
	std::vector<MetadataRow> rows = {{"S1", "v", "a"}};
	REQUIRE_THROWS_AS(BuildGroupings(rows, {}, {"v"}), std::invalid_argument);
}

TEST_CASE("BuildGroupings throws on duplicate (sample, variable) with conflicting values", "[unifrac][metadata]") {
	// WHY: silent last-write-wins on a botched join would change PERMANOVA's
	// F-stat without any user-visible signal. The error names both keys so
	// the user can fix upstream.
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"},
	    {"S2", "treatment", "trt"},
	    {"S1", "treatment", "trt"}, // conflicting with the earlier row
	};
	std::vector<std::string> samples = {"S1", "S2"};

	REQUIRE_THROWS_WITH(BuildGroupings(rows, samples, {"treatment"}),
	                    Catch::Matchers::ContainsSubstring("S1") && Catch::Matchers::ContainsSubstring("treatment") &&
	                        Catch::Matchers::ContainsSubstring("ctrl") && Catch::Matchers::ContainsSubstring("trt"));
}

TEST_CASE("BuildGroupings tolerates idempotent duplicate metadata rows", "[unifrac][metadata]") {
	// Duplicates that agree are a no-op (common when concatenating overlapping
	// metadata exports); only conflicts are errors.
	std::vector<MetadataRow> rows = {
	    {"S1", "treatment", "ctrl"}, {"S2", "treatment", "trt"}, {"S1", "treatment", "ctrl"}, // identical to the first
	};
	std::vector<std::string> samples = {"S1", "S2"};

	auto groupings = BuildGroupings(rows, samples, {"treatment"});
	REQUIRE(groupings.size() == 1);
	REQUIRE(groupings[0].n_groups == 2);
	REQUIRE(groupings[0].labels == std::vector<uint32_t> {0, 1});
}
