#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <vector>

#include "simulate_resemblance.hpp"

using Catch::Matchers::ContainsSubstring;

// Pins the DETERMINISTIC / structural guarantees of the cluster model. Exact
// perturbed values depend on the (non-portable) std RNG stream, so those are
// covered statistically by the SQL parity test; here we pin what must hold on
// every platform: determinism, depth, cluster labeling, and the multiplicative
// zero-preservation that distinguishes '*sample'/'clip' from additive noise.

TEST_CASE("SimulateCluster is deterministic for a fixed seed", "[simulate]") {
	const std::vector<double> abund = {0.4, 0.25, 0.2, 0.1, 0.05};
	const std::vector<int32_t> sizes = {2, 3};
	auto a = miint::SimulateCluster(abund, sizes, 500, 1.0, 0.5, "*sample", "clip", 123);
	auto b = miint::SimulateCluster(abund, sizes, 500, 1.0, 0.5, "*sample", "clip", 123);
	REQUIRE(a.size() == b.size());
	CHECK(a.sample_id == b.sample_id);
	CHECK(a.otu_id == b.otu_id);
	CHECK(a.count == b.count);
	CHECK(a.ground_truth == b.ground_truth);
}

TEST_CASE("SimulateCluster labels samples by cluster and enforces depth", "[simulate]") {
	const std::vector<double> abund = {0.4, 0.25, 0.2, 0.1, 0.05};
	const std::vector<int32_t> sizes = {2, 3}; // 5 samples: clusters 0,0,1,1,1
	auto coo = miint::SimulateCluster(abund, sizes, 500, 1.0, 0.5, "*sample", "clip", 7);

	std::vector<int64_t> depth(5, 0);
	std::vector<int> cluster_of(5, -1);
	for (size_t i = 0; i < coo.size(); i++) {
		CHECK(coo.count[i] >= 1); // nonzero-only COO
		auto sid = coo.sample_id[i];
		REQUIRE(sid >= 0);
		REQUIRE(sid < 5);
		depth[static_cast<size_t>(sid)] += coo.count[i];
		cluster_of[static_cast<size_t>(sid)] = static_cast<int>(coo.ground_truth[i]);
	}
	for (int s = 0; s < 5; s++) {
		CHECK(depth[static_cast<size_t>(s)] == 500);
	}
	// cluster_sizes [2,3] => sample->cluster mapping 0,0,1,1,1.
	CHECK(cluster_of == std::vector<int> {0, 0, 1, 1, 1});
}

TEST_CASE("SimulateCluster multiplicative/clip preserves structural zeros", "[simulate]") {
	// Middle species has zero base abundance: 0 * noise -> 0, clip keeps 0, renorm
	// keeps 0, never sampled. So otu_id 1 must be entirely absent from the output.
	const std::vector<double> abund = {0.5, 0.0, 0.5};
	const std::vector<int32_t> sizes = {3, 3};
	auto coo = miint::SimulateCluster(abund, sizes, 500, 1.0, 0.5, "*sample", "clip", 7);
	for (auto o : coo.otu_id) {
		CHECK(o != 1);
	}
}

TEST_CASE("SimulateCluster fails loud on all-zero abundances (NaN-safe)", "[simulate]") {
	// All-zero base would make root = 0/0 = NaN; the guard must reject it cleanly
	// (the wrapper blocks this too, but the pure core must be self-safe).
	const std::vector<double> zero = {0.0, 0.0, 0.0};
	const std::vector<int32_t> sizes = {2, 2};
	CHECK_THROWS_WITH(miint::SimulateCluster(zero, sizes, 100, 1.0, 0.5, "*sample", "clip", 7),
	                  ContainsSubstring("positive"));
}

TEST_CASE("SimulateCluster rescale collapses a single-species node -> throws", "[simulate]") {
	// 'rescale' subtracts the minimum; with one species that element IS the min,
	// so it becomes exactly 0 -> sum 0 -> degenerate. Deterministic regardless of seed.
	const std::vector<double> one_sp = {1.0};
	const std::vector<int32_t> sizes = {2, 2};
	CHECK_THROWS_WITH(miint::SimulateCluster(one_sp, sizes, 100, 1.0, 0.5, "*sample", "rescale", 7),
	                  ContainsSubstring("degenerate"));
}
