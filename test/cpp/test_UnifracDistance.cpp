#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "NewickTree.hpp"
#include "concurrent_test_util.hpp"
#include "unifrac_bptree.hpp"
#include "unifrac_distance.hpp"
#include "api.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_subsample_bridge.hpp"
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

// A fixture whose subsampling is genuinely random: every sample totals 40 counts
// spread over 4 features, so a draw to depth 8 has many possible outcomes.
// MakeLiveFixture cannot serve here — its samples total exactly 2, so a depth-2
// draw is forced and would hide any RNG difference rather than expose it.
LiveFixture MakeRandomizableFixture() {
	std::vector<CooRow> rows;
	for (const char *s : {"S1", "S2", "S3", "S4"}) {
		for (const char *f : {"A", "B", "C", "D"}) {
			rows.push_back({s, f, 10.0});
		}
	}
	auto biom = UnifracSupportBiomView::FromCoo(std::move(rows));
	auto tree = miint::NewickTree::parse("((A:0.1,B:0.2)Int1:0.3,(C:0.4,D:0.5)Int2:0.6);");
	auto bptree = UnifracBptreeView::FromNewickTree(tree);
	return {std::move(biom), std::move(bptree)};
}

// A view's CSR payload, flattened so two subsample draws can be compared cell
// for cell.
std::vector<double> BiomCells(const UnifracSupportBiomView &v) {
	const support_biom_t *b = v.support_biom();
	std::vector<double> out;
	out.push_back(static_cast<double>(b->n_obs));
	out.push_back(static_cast<double>(b->n_samples));
	const uint32_t nnz = b->indptr[b->n_obs];
	for (uint32_t i = 0; i < nnz; ++i) {
		out.push_back(static_cast<double>(b->indices[i]));
		out.push_back(b->data[i]);
	}
	return out;
}

} // namespace

TEST_CASE("UnifracDistanceMatrix::Compute leaves libssu's process-global RNG untouched", "[unifrac][distance]") {
	// WHY THIS MATTERS: Compute used to resolve its seed by calling
	// ssu_set_random_seed, which mutates state shared by the entire process — and
	// per upstream's README ("Two traps worth knowing") that one call reseeds
	// scikit-bio-binaries' global generator as well. Two consequences followed:
	// an unrelated later unseeded draw silently changed, and two concurrent
	// Computes clobbered each other's seeding. The latter is precisely why the
	// call needed a process-wide mutex, so this property is the PRECONDITION for
	// dropping that lock — not a nicety. Passing the seed per call to
	// one_off_matrix_inmem_fp32_v4 is what buys it.
	//
	// The probe is an UNSEEDED BridgeSubsample, which by contract draws from that
	// global generator. Seed the global, draw, then require an interposed seeded
	// Compute not to change what the same draw produces.
	auto fixture = MakeRandomizableFixture();

	// Hold the OpenMP width fixed across both draws. A seeded libssu draw is
	// distributed across the OpenMP team, so team size selects the result; pinning
	// it here isolates the RNG effect this test is actually about.
	miint::unifrac::OmpThreadPin pin(1);

	const auto unseeded_draw = [&]() {
		return BiomCells(miint::unifrac::BridgeSubsample(fixture.biom, /*subsample_depth*/ 8,
		                                                 /*with_replacement*/ false, /*seed*/ -1));
	};

	ssu_set_random_seed(1234);
	const std::vector<double> reference = unseeded_draw();

	// Control: the probe must be sound before it can indict anything. Re-seeding
	// the global and drawing again with nothing in between has to reproduce, or a
	// failure below would just mean "the draw is noisy", not "Compute reseeded it".
	ssu_set_random_seed(1234);
	REQUIRE(unseeded_draw() == reference);

	ssu_set_random_seed(1234);
	auto dist = UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree, "weighted_normalized_fp32",
	                                           /*variance_adjust*/ false, /*alpha*/ 1.0,
	                                           /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
	                                           /*subsample_depth*/ 0, /*subsample_with_replacement*/ false,
	                                           /*seed*/ 42, /*n_threads*/ 1);
	REQUIRE(dist.n_samples() == 4); // the compute really ran

	REQUIRE(unseeded_draw() == reference);
}

TEST_CASE("concurrent seeded libssu UniFrac reproduces the serial result exactly", "[unifrac][distance][omp]") {
	// THE CAPABILITY CLAIM, and the reason UnifracDistanceMatrix::Compute no longer
	// takes a process-wide lock: several UniFrac computes may run at once in one
	// process. That lock was not paranoia — before upstream's a9cea63 every compute
	// bracketed itself with register_report_status()/remove_report_status(), which
	// calloc'd then free+NULL'd a process-global pointer the inner loop dereferenced
	// unchecked, so two concurrent computes were a use-after-free (ASan reproduced
	// SEGV on 0x0). The flags are a static atomic array now.
	//
	// Bit-identical, not merely close: anything less would mean concurrent calls
	// perturbed each other's numerics, which is exactly what a shared RNG or a
	// fought-over thread count would cause. Subsampling is ON so the seeded draw is
	// exercised too — that is the part that would silently diverge if the seed were
	// ever resolved globally again.
	auto fixture = MakeRandomizableFixture();
	const int seed = 7;

	const auto compute_once = [&]() {
		auto dist = UnifracDistanceMatrix::Compute(fixture.biom, fixture.bptree, "weighted_normalized_fp32",
		                                           /*variance_adjust*/ false, /*alpha*/ 1.0,
		                                           /*bypass_tips*/ false, /*normalize_sample_counts*/ true,
		                                           /*subsample_depth*/ 8, /*subsample_with_replacement*/ false, seed,
		                                           /*n_threads*/ 1);
		const uint32_t n = dist.n_samples();
		return std::vector<float>(dist.matrix(), dist.matrix() + static_cast<size_t>(n) * n);
	};

	const std::vector<float> serial = compute_once();
	REQUIRE(serial.size() == 16); // 4 samples survive depth 8, so the draw really ran

	const size_t workers = 8;
	std::vector<std::vector<float>> got(workers);
	miint::test::RunConcurrently(workers, [&](size_t t) { got[t] = compute_once(); });

	for (size_t t = 0; t < workers; ++t) {
		REQUIRE(got[t] == serial);
	}
}

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
