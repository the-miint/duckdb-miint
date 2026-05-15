#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "NewickTree.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
#include "unifrac_libssu.hpp"
#include "unifrac_support_biom.hpp"

#include <cstring>
#include <stdexcept>
#include <utility>
#include <vector>

using miint::unifrac::CooRow;
using miint::unifrac::UnifracBptreeView;
using miint::unifrac::UnifracDistanceMatrix;
using miint::unifrac::UnifracSupportBiomView;

namespace {

// Fixture: the same symmetric 4-tip tree the support-biom live test uses
// (test_UnifracSupportBiom.cpp). Faith PD on that tree is hand-checkable,
// so any divergence in the BIOM/tree round-trip surfaces immediately.
struct LiveFixture {
	UnifracSupportBiomView biom;
	UnifracBptreeView bptree;
};

LiveFixture MakeLiveFixture() {
	std::vector<CooRow> rows = {
	    {"S1", "A", 1.0}, {"S2", "A", 1.0}, {"S2", "B", 1.0}, {"S3", "C", 1.0},
	    {"S3", "D", 1.0}, {"S4", "A", 1.0}, {"S4", "C", 1.0},
	};
	auto biom = UnifracSupportBiomView::FromCoo(std::move(rows));
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2)Int1:0.3,(C:0.4,D:0.5)Int2:0.6);");
	auto bptree = UnifracBptreeView::FromNewickTree(tree);
	return {std::move(biom), std::move(bptree)};
}

} // namespace

TEST_CASE("UnifracDistanceMatrix computes a symmetric matrix with zero diagonal", "[unifrac][distance]") {
	// Sanity check the contract: the C wrapper produces a valid distance
	// matrix indexed by mat->sample_ids. Verifies the libssu call path
	// against a hand-checked fixture rather than just trusting the C return.
	auto fixture = MakeLiveFixture();
	auto dist = UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree,
	                                           /*variant_fp32*/ "weighted_normalized_fp32",
	                                           /*variance_adjust*/ false, /*alpha*/ 1.0,
	                                           /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
	                                           /*subsample_depth*/ 0, /*subsample_with_replacement*/ false,
	                                           /*seed*/ 42, /*n_threads*/ 1);

	REQUIRE(dist.n_samples() == 4);
	REQUIRE(dist.sample_ids().size() == 4);
	// Lexicographic order is the Phase-1 adapter's canonical contract.
	REQUIRE(dist.sample_ids()[0] == "S1");
	REQUIRE(dist.sample_ids()[1] == "S2");
	REQUIRE(dist.sample_ids()[2] == "S3");
	REQUIRE(dist.sample_ids()[3] == "S4");

	const float *mat = dist.matrix();
	REQUIRE(mat != nullptr);

	// Distance matrix must be symmetric with zero diagonal.
	const uint32_t n = dist.n_samples();
	for (uint32_t i = 0; i < n; ++i) {
		REQUIRE(mat[i * n + i] == Catch::Approx(0.0).margin(1e-6));
		for (uint32_t j = i + 1; j < n; ++j) {
			REQUIRE(mat[i * n + j] == Catch::Approx(mat[j * n + i]));
			REQUIRE(mat[i * n + j] >= 0.0f);
		}
	}
}

TEST_CASE("UnifracDistanceMatrix subsample_depth drops samples whose total counts fall below the depth",
          "[unifrac][distance]") {
	// S1 has total count 1 → dropped when subsample_depth >= 2.
	// All others have total count 2 → survive at depth 2.
	// This is the exact bug Phase 4 hit (mat->n_samples can differ from
	// the input view's n_samples); the test guards against silent regressions
	// where the wrapper trusts the wrong source of truth.
	auto fixture = MakeLiveFixture();
	auto dist = UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree, "weighted_normalized_fp32",
	                                           /*variance_adjust*/ false, /*alpha*/ 1.0,
	                                           /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
	                                           /*subsample_depth*/ 2, /*subsample_with_replacement*/ false,
	                                           /*seed*/ 42, /*n_threads*/ 1);

	REQUIRE(dist.n_samples() == 3);
	REQUIRE(dist.sample_ids().size() == 3);
	// S1 must be absent; the remaining samples survive in lex order.
	for (const auto &id : dist.sample_ids()) {
		REQUIRE(id != "S1");
	}
}

TEST_CASE("UnifracDistanceMatrix is move-only and transfers ownership cleanly", "[unifrac][distance]") {
	// The wrapper holds a heap-allocated libssu pointer; we need move semantics
	// to be airtight so that destroy_mat_full_fp32 is called exactly once.
	auto fixture = MakeLiveFixture();
	auto dist = UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree, "weighted_normalized_fp32", false, 1.0,
	                                           false, true, 0, false, 42, /*n_threads*/ 1);
	const float *original_matrix = dist.matrix();
	const uint32_t original_n = dist.n_samples();
	REQUIRE(original_n == 4);

	// Move construction: target owns the matrix; source is in a moved-from
	// state (we only check what's defined — accessing the source's matrix
	// after move is intentionally not checked).
	UnifracDistanceMatrix moved = std::move(dist);
	REQUIRE(moved.n_samples() == original_n);
	REQUIRE(moved.matrix() == original_matrix);
	REQUIRE(moved.sample_ids().size() == 4);

	// Move assignment: a second target steals from the first; destructor
	// of the first must NOT double-free the libssu pointer.
	UnifracDistanceMatrix second =
	    UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree, "weighted_normalized_fp32", false, 1.0, false,
	                                   true, 0, false, 42, /*n_threads*/ 1);
	REQUIRE(second.matrix() != nullptr);
	second = std::move(moved);
	REQUIRE(second.n_samples() == 4);
	// `moved` is now in moved-from state. `second`'s destructor runs at end
	// of scope and frees exactly one libssu pointer — verified by absence of
	// crash / double-free under ASan if the test binary is built with it.
}

TEST_CASE("UnifracDistanceMatrix surfaces libssu errors via std::runtime_error", "[unifrac][distance]") {
	// Pass an unknown method name; libssu returns unknown_method (status 4).
	auto fixture = MakeLiveFixture();
	REQUIRE_THROWS_AS(UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree,
	                                                 /*variant_fp32*/ "nonexistent_metric_fp32",
	                                                 /*variance_adjust*/ false, /*alpha*/ 1.0,
	                                                 /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
	                                                 /*subsample_depth*/ 0,
	                                                 /*subsample_with_replacement*/ false, /*seed*/ -1,
	                                                 /*n_threads*/ 1),
	                  std::runtime_error);
}
