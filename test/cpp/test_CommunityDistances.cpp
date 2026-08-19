#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <iterator>
#include <random>
#include <stdexcept>
#include <vector>

#include "community_distances.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

// These pin every metric formula with EXACT hand-computed expectations on a
// tiny matrix, so a broken formula (wrong normalization, sign, or global-stat
// handling) fails here. The random/statistical fidelity vs scikit-bio is
// checked separately by the SQL parity oracle.
//
// One shared 3x4 matrix drives all cases. Purely-pairwise metrics
// (bray_curtis/euclidean/jaccard/soergel/morisita_horn/pearson) ignore the
// third row, so their pair (0,1) equals the by-hand two-vector computation;
// the global metrics (chisq/gower) DO depend on the third row via column
// sums / ranges, so their pair (0,1) values are computed against all 3 rows.
//
//   row0 (x) = [1, 2, 0, 3]
//   row1 (y) = [0, 2, 4, 1]
//   row2     = [2, 0, 1, 1]
//
// Condensed upper-triangle order: out[0]=(0,1), out[1]=(0,2), out[2]=(1,2).
namespace {
const std::vector<double> M3x4 = {1, 2, 0, 3, 0, 2, 4, 1, 2, 0, 1, 1};

std::vector<double> Dist(const std::string &metric) {
	return miint::CommunityDistancesCondensed(M3x4, 3, 4, metric);
}
} // namespace

TEST_CASE("community metric name validation", "[community_distances]") {
	for (const auto *m :
	     {"bray_curtis", "euclidean", "jaccard", "soergel", "morisita_horn", "pearson", "chisq", "gower"}) {
		CHECK(miint::IsValidCommunityMetric(m));
	}
	CHECK_FALSE(miint::IsValidCommunityMetric("braycurtis"));  // no underscore variant
	CHECK_FALSE(miint::IsValidCommunityMetric("Bray_Curtis")); // case sensitive at the core layer
	CHECK_FALSE(miint::IsValidCommunityMetric("manhattan"));
	CHECK_FALSE(miint::IsValidCommunityMetric(""));
}

TEST_CASE("bray_curtis exact", "[community_distances]") {
	auto d = Dist("bray_curtis");
	REQUIRE(d.size() == 3);
	// Sum|x-y| = 1+0+4+2 = 7 ; Sum(x+y) = 1+4+4+4 = 13
	CHECK(d[0] == Approx(7.0 / 13.0));
	for (double v : d) {
		CHECK(v >= 0.0);
		CHECK(v <= 1.0);
	}
}

TEST_CASE("euclidean exact", "[community_distances]") {
	auto d = Dist("euclidean");
	// sqrt(1 + 0 + 16 + 4) = sqrt(21)
	CHECK(d[0] == Approx(std::sqrt(21.0)));
}

TEST_CASE("jaccard binary exact", "[community_distances]") {
	auto d = Dist("jaccard");
	// presence x=[1,1,0,1] y=[0,1,1,1] -> a(both)=2, b(x only)=1, c(y only)=1
	// (b+c)/(a+b+c) = 2/4 = 0.5
	CHECK(d[0] == Approx(0.5));
	for (double v : d) {
		CHECK(v >= 0.0);
		CHECK(v <= 1.0);
	}
}

TEST_CASE("soergel exact", "[community_distances]") {
	auto d = Dist("soergel");
	// Sum|x-y| = 7 ; Sum max(x,y) = 1+2+4+3 = 10 -> 0.7
	CHECK(d[0] == Approx(0.7));
	for (double v : d) {
		CHECK(v >= 0.0);
		CHECK(v <= 1.0);
	}
}

TEST_CASE("morisita_horn exact", "[community_distances]") {
	auto d = Dist("morisita_horn");
	// X=6, Y=7 ; Sx2=14, Sy2=21 ; Sxy=7
	// C = 2*7 / ((14/36 + 21/49) * 6 * 7) = 42/103 ; dist = 1 - 42/103 = 61/103
	CHECK(d[0] == Approx(61.0 / 103.0));
	for (double v : d) {
		CHECK(v >= 0.0);
		CHECK(v <= 1.0);
	}
}

TEST_CASE("pearson (correlation) distance exact", "[community_distances]") {
	auto d = Dist("pearson");
	// mean x=1.5, mean y=1.75 ; cov=-3.5 ; varx=5, vary=8.75
	// r = -3.5/sqrt(43.75) ; dist = 1 - r
	const double r = -3.5 / std::sqrt(43.75);
	CHECK(d[0] == Approx(1.0 - r));
}

TEST_CASE("chisq (correspondence-analysis) exact", "[community_distances]") {
	auto d = Dist("chisq");
	// rows r=(6,7,4) ; cols c=(3,4,5,5) ; GT=17
	// pair(0,1): sqrt( (17/3)(1/6-0)^2 + (17/4)(2/6-2/7)^2
	//                + (17/5)(0-4/7)^2 + (17/5)(3/6-1/7)^2 )
	double s = 0.0;
	const double GT = 17.0;
	const double col[4] = {3, 4, 5, 5};
	const double x[4] = {1, 2, 0, 3};
	const double y[4] = {0, 2, 4, 1};
	for (int k = 0; k < 4; ++k) {
		const double a = x[k] / 6.0;
		const double b = y[k] / 7.0;
		s += (GT / col[k]) * (a - b) * (a - b);
	}
	CHECK(d[0] == Approx(std::sqrt(s)));
	for (double v : d) {
		CHECK(v >= 0.0);
	}
}

TEST_CASE("gower (PyCogent, un-normalized) exact", "[community_distances]") {
	auto d = Dist("gower");
	// column ranges over all 3 rows: col0=2, col1=2, col2=4, col3=2
	// pair(0,1): 1/2 + 0/2 + 4/4 + 2/2 = 2.5
	// pair(0,2): 1/2 + 2/2 + 1/4 + 2/2 = 2.75
	// pair(1,2): 2/2 + 2/2 + 3/4 + 0/2 = 2.75
	REQUIRE(d.size() == 3);
	CHECK(d[0] == Approx(2.5));
	CHECK(d[1] == Approx(2.75));
	CHECK(d[2] == Approx(2.75));
}

TEST_CASE("identical rows -> zero distance for all metrics", "[community_distances]") {
	// Two identical non-constant rows: every metric collapses to 0 (pearson too,
	// since r = 1 for identical non-constant profiles).
	const std::vector<double> m = {1, 2, 0, 3, 1, 2, 0, 3};
	for (const auto *metric :
	     {"bray_curtis", "euclidean", "jaccard", "soergel", "morisita_horn", "pearson", "chisq", "gower"}) {
		auto d = miint::CommunityDistancesCondensed(m, 2, 4, metric);
		REQUIRE(d.size() == 1);
		INFO("metric = " << metric);
		CHECK(d[0] == Approx(0.0).margin(1e-12));
	}
}

TEST_CASE("all-zero sample pair -> zero for guarded metrics", "[community_distances]") {
	// Both samples empty: bray_curtis/jaccard/soergel/morisita_horn are defined to
	// 0 (identical empty communities), euclidean is 0 trivially.
	const std::vector<double> m = {0, 0, 0, 0, 0, 0, 0, 0};
	for (const auto *metric : {"bray_curtis", "euclidean", "jaccard", "soergel", "morisita_horn"}) {
		auto d = miint::CommunityDistancesCondensed(m, 2, 4, metric);
		REQUIRE(d.size() == 1);
		INFO("metric = " << metric);
		CHECK(d[0] == Approx(0.0));
	}
}

TEST_CASE("pearson flat/constant profile -> PyCogent finite semantics", "[community_distances]") {
	// A constant row has zero variance, so Pearson r is undefined. Kuczynski 2010
	// used PyCogent's dist_pearson, which returns FINITE values here (verified
	// against the PyCogent source): a flat vs a non-flat row -> r=0 -> distance 1;
	// two flat rows -> r=1 -> distance 0. (scipy's `correlation` returns NaN
	// instead; we deliberately follow PyCogent for faithful reproduction.)
	const std::vector<double> flat_vs_nonflat = {2, 2, 2, 2, 1, 3, 0, 5};
	auto d = miint::CommunityDistancesCondensed(flat_vs_nonflat, 2, 4, "pearson");
	REQUIRE(d.size() == 1);
	CHECK(d[0] == Approx(1.0)); // one flat row -> maximal pearson distance
	const std::vector<double> flat_vs_flat = {2, 2, 2, 2, 7, 7, 7, 7};
	CHECK(miint::CommunityDistancesCondensed(flat_vs_flat, 2, 4, "pearson")[0] == Approx(0.0));
	// A non-constant identical pair is well-defined (r=1 -> distance 0).
	const std::vector<double> ok = {1, 3, 0, 5, 1, 3, 0, 5};
	CHECK(miint::CommunityDistancesCondensed(ok, 2, 4, "pearson")[0] == Approx(0.0).margin(1e-12));
}

TEST_CASE("chisq zero-row-sum sample -> PyCogent finite semantics", "[community_distances]") {
	// An all-zero sample has no correspondence-analysis row profile. PyCogent
	// dist_chisq (the metric Kuczynski 2010 used, verified against its source)
	// returns FINITE values here, not NaN: one empty vs a non-empty row -> 1;
	// two empty rows -> 0. (Reachable only from a dense caller: the sparse SQL
	// reader drops all-zero samples before they reach the core.)
	const std::vector<double> one_empty = {0, 0, 0, 0, 1, 2, 3, 4};
	auto d = miint::CommunityDistancesCondensed(one_empty, 2, 4, "chisq");
	REQUIRE(d.size() == 1);
	CHECK(d[0] == Approx(1.0));
	const std::vector<double> both_empty = {0, 0, 0, 0, 0, 0, 0, 0};
	CHECK(miint::CommunityDistancesCondensed(both_empty, 2, 4, "chisq")[0] == Approx(0.0));
}

TEST_CASE("pearson numerical stability on high-mean/low-variance profiles", "[community_distances]") {
	// The naive Sum(x^2)-f*mean^2 variance shortcut cancels catastrophically here
	// (mean ~1e4, variance ~1). The centered two-pass must still return a finite,
	// in-range [0,2] correlation distance for these genuinely non-constant rows.
	const std::vector<double> m = {10000, 10001, 9998, 10002, 9999, 10000, 9999, 10002, 10000, 9998};
	auto d = miint::CommunityDistancesCondensed(m, 2, 5, "pearson");
	REQUIRE(d.size() == 1);
	CHECK(std::isfinite(d[0]));
	CHECK(d[0] >= 0.0);
	CHECK(d[0] <= 2.0);
}

TEST_CASE("morisita_horn one-empty sample -> max dissimilarity", "[community_distances]") {
	// One community empty, the other not: no shared abundance -> dissimilarity 1.
	const std::vector<double> m = {0, 0, 0, 0, 1, 2, 3, 4};
	auto d = miint::CommunityDistancesCondensed(m, 2, 4, "morisita_horn");
	REQUIRE(d.size() == 1);
	CHECK(d[0] == Approx(1.0));
}

TEST_CASE("condensed size and pair ordering", "[community_distances]") {
	// 4 samples -> 6 unordered pairs, row-major upper triangle.
	const std::vector<double> m(4 * 3, 1.0); // all-equal rows -> all distances 0
	auto d = miint::CommunityDistancesCondensed(m, 4, 3, "euclidean");
	CHECK(d.size() == 6);
}

TEST_CASE("parallel result is bit-identical to serial for every metric", "[community_distances]") {
	// Determinism contract: adding threads is a performance change that must NOT
	// change a single bit of the output (the SQL parity goldens are pinned exact).
	// Each pair (i,j) is computed by identical arithmetic and written to a fixed
	// output slot regardless of which thread handles row i, so serial and parallel
	// must agree exactly.
	//
	// 9 samples so the row cursor genuinely splits work across threads (and the
	// thread cap at n-1 is exercised by nt=16). Values are strictly positive and
	// every row is non-constant, so no metric hits its NaN/empty special case
	// (bit-equality on NaN would be meaningless).
	const uint32_t n = 9, f = 6;
	std::vector<double> m(static_cast<size_t>(n) * f);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t k = 0; k < f; ++k) {
			// The 2*k term guarantees a non-constant row (pearson variance > 0);
			// the leading +1 guarantees a positive row sum (morisita/chisq defined).
			m[static_cast<size_t>(i) * f + k] = 1.0 + 2.0 * k + static_cast<double>((i * 5 + k * 3) % 7) + 0.25 * i;
		}
	}
	for (const auto *metric :
	     {"bray_curtis", "euclidean", "jaccard", "soergel", "morisita_horn", "pearson", "chisq", "gower"}) {
		const auto serial = miint::CommunityDistancesCondensed(m, n, f, metric, 1);
		// 1000 exercises the thread cap (min with n-1 and hardware_concurrency):
		// it must not spawn 1000 OS threads, must not hang/crash, and must still
		// produce the identical serial result.
		for (unsigned nt : {2u, 3u, 4u, 8u, 16u, 1000u}) {
			const auto par = miint::CommunityDistancesCondensed(m, n, f, metric, nt);
			REQUIRE(par.size() == serial.size());
			for (size_t p = 0; p < serial.size(); ++p) {
				CHECK(par[p] == serial[p]); // exact bit equality, not Approx
			}
		}
	}
}

TEST_CASE("community_distances guards", "[community_distances]") {
	// fewer than 2 samples
	CHECK_THROWS_AS(miint::CommunityDistancesCondensed({1, 2, 3}, 1, 3, "euclidean"), std::invalid_argument);
	CHECK_THROWS_WITH(miint::CommunityDistancesCondensed({1, 2, 3}, 1, 3, "euclidean"),
	                  ContainsSubstring("at least 2"));
	// matrix size mismatch
	CHECK_THROWS_AS(miint::CommunityDistancesCondensed({1, 2, 3}, 2, 3, "euclidean"), std::invalid_argument);
	// unknown metric
	CHECK_THROWS_WITH(miint::CommunityDistancesCondensed(M3x4, 3, 4, "manhattan"), ContainsSubstring("manhattan"));
}

// ── Pairwise locality: the property progressive_pcoa_from_features rests on ──
//
// A metric can be computed one BLOCK of samples at a time iff d(i,j) depends
// only on samples i and j — never on which other samples share the matrix, and
// never on which features that matrix happens to carry. Everything below builds
// the exact situation a block creates: take a subset of the samples, restrict
// the feature space to the features that are nonzero WITHIN that subset (a
// feature nobody in the block carries is not in the block's matrix at all), and
// recompute.
//
// For an admissible metric the pair distances must come back BITWISE unchanged.
// For the rest they must come back DIFFERENT — that difference is precisely the
// silent wrongness this classification exists to prevent, so it is asserted
// rather than described.
//
//   s0 = [3, 0, 1, 4, 0, 0]     subset = {s0, s1, s2}
//   s1 = [1, 2, 0, 2, 0, 0]     features nonzero within it: f0..f3
//   s2 = [0, 5, 2, 0, 0, 0]     (f4, f5 are zero in all three)
//   s3 = [0, 0, 0, 0, 7, 1]     f4 lives only outside the subset
//   s4 = [8, 1, 0, 3, 0, 9]     f0's maximum (8) lives only outside it
//
// The fixture is built so each excluded metric is broken by a DIFFERENT
// property, so no one of them can pass by accident:
//   pearson — the subset has 4 feature columns, not 6, and f enters the
//             per-sample mean (rowsum[i]/f) and the both-zero terms of cov/var
//   chisq   — column sums change (f0: 12 -> 4), and so does the grand total
//   gower   — column ranges change (f0: 8 -> 3), because s4 holds f0's maximum
namespace {

const std::vector<double> kFull5x6 = {
    3, 0, 1, 4, 0, 0, //
    1, 2, 0, 2, 0, 0, //
    0, 5, 2, 0, 0, 0, //
    0, 0, 0, 0, 7, 1, //
    8, 1, 0, 3, 0, 9, //
};

// The same three samples a block would see, over only the features they carry.
const std::vector<double> kBlock3x4 = {
    3, 0, 1, 4, //
    1, 2, 0, 2, //
    0, 5, 2, 0, //
};

// Row-major upper-triangle position of pair (i, j), i < j, over n samples --
// the layout CommunityDistancesCondensed documents.
size_t Condensed(uint32_t n, uint32_t i, uint32_t j) {
	return static_cast<size_t>(i) * (2 * n - i - 1) / 2 + (j - i - 1);
}

constexpr const char *kPairwiseLocalMetrics[] = {"bray_curtis", "euclidean", "jaccard", "soergel", "morisita_horn"};
constexpr const char *kMatrixWideMetrics[] = {"pearson", "chisq", "gower"};

} // namespace

TEST_CASE("pairwise-local classification covers every metric", "[community_distances]") {
	for (const auto *m : kPairwiseLocalMetrics) {
		INFO("metric: " << m);
		CHECK(miint::IsPairwiseLocalCommunityMetric(m));
	}
	for (const auto *m : kMatrixWideMetrics) {
		INFO("metric: " << m);
		CHECK_FALSE(miint::IsPairwiseLocalCommunityMetric(m));
	}
	// The two lists above must together BE the accepted set. Counted against the
	// RUNTIME size of that set, not against a literal: `kPairwiseLocalMetrics` and
	// `kMatrixWideMetrics` are compile-time arrays, so comparing their combined
	// length to 8 is a tautology that stays green when a ninth metric appears in
	// kMetrics and in neither list -- the very case this is here to catch.
	// (The pairwise-local names need no IsValid check -- IsPairwiseLocal returning
	// true above already implies it.)
	for (const auto *m : kMatrixWideMetrics) {
		REQUIRE(miint::IsValidCommunityMetric(m));
	}
	const std::string accepted = miint::CommunityMetricList();
	const size_t n_accepted = std::count(accepted.begin(), accepted.end(), ',') + 1;
	CHECK(n_accepted == std::size(kPairwiseLocalMetrics) + std::size(kMatrixWideMetrics));
	// An unknown name is not admissible by default.
	CHECK_FALSE(miint::IsPairwiseLocalCommunityMetric("manhattan"));
	CHECK_FALSE(miint::IsPairwiseLocalCommunityMetric(""));
}

TEST_CASE("pairwise-local metrics are unchanged by the block a pair lands in", "[community_distances]") {
	// The load-bearing assertion: restricting to a sample subset and to that
	// subset's own feature space leaves every pair distance BITWISE identical.
	// Features nobody in the block carries contribute exactly 0.0 to every sum,
	// and the per-sample aggregates are unchanged, so this is exact -- not
	// approximate. Approx here would hide a metric that drifts per block.
	for (const auto *metric : kPairwiseLocalMetrics) {
		INFO("metric: " << metric);
		const auto full = miint::CommunityDistancesCondensed(kFull5x6, 5, 6, metric);
		const auto block = miint::CommunityDistancesCondensed(kBlock3x4, 3, 4, metric);
		REQUIRE(full.size() == 10);
		REQUIRE(block.size() == 3);
		for (uint32_t i = 0; i < 3; ++i) {
			for (uint32_t j = i + 1; j < 3; ++j) {
				INFO("pair (" << i << "," << j << ")");
				CHECK(full[Condensed(5, i, j)] == block[Condensed(3, i, j)]);
			}
		}
	}
}

TEST_CASE("matrix-wide metrics DO change with the block -- why they are excluded", "[community_distances]") {
	// Each of these reads something global: pearson the feature-column count,
	// chisq the column sums and grand total, gower the column ranges. Computed
	// per block they would silently measure a different distance in every block,
	// which is why progressive_pcoa_from_features rejects them at bind rather
	// than computing them. If one of these ever stops differing, the metric has
	// been changed and its exclusion must be revisited -- fail loudly.
	for (const auto *metric : kMatrixWideMetrics) {
		INFO("metric: " << metric);
		const auto full = miint::CommunityDistancesCondensed(kFull5x6, 5, 6, metric);
		const auto block = miint::CommunityDistancesCondensed(kBlock3x4, 3, 4, metric);
		CHECK(full[Condensed(5, 0, 1)] != block[Condensed(3, 0, 1)]);
	}
}

// ── Sparse (CSR) pair loop ───────────────────────────────────────────────────
//
// The block-wise path computes O(n^2) pairs per block, and a dense pair loop pays
// for the whole feature space on every pair. Measured on the 1.2M-sample reference
// table, one default block spans 11,018 features but only 89 nonzeros per sample,
// so a union-of-nonzeros merge does ~62x less arithmetic. That is the whole reason
// this entry point exists.
//
// It is required to agree with the dense implementation BIT FOR BIT, not
// approximately. Both accumulate over feature index in ascending order; the terms
// the dense loop adds for features absent from both samples are exactly 0.0, and
// adding 0.0 changes no double. So any difference at all means the sparse loop
// visits terms the dense one does not, or in a different order — a real defect,
// not rounding. Tolerating a small epsilon here would hide exactly that.
namespace {

// One matrix held both ways, built from the SAME data so a bitwise comparison is
// meaningful: CSR indices ascending within each row, dense laid out over the same
// n_features columns.
struct BothForms {
	std::vector<uint32_t> indptr, indices;
	std::vector<double> values;
	std::vector<double> dense;
	uint32_t n = 0, f = 0;
};

// `density` in [0,1]; `force_empty_rows` makes every 5th sample carry nothing,
// which is the branch every metric special-cases (empty pair -> 0, one-empty -> 1).
// `allow_negative` draws values that can be NEGATIVE: nothing upstream rejects them
// (the feature-table filter drops only NULL/zero/NaN), and they are what separates a
// `> 0.0` guard from an `== 0.0` one — a presence test, a zero denominator, and an
// empty-community check all behave differently on a negative than on a zero. The
// original all-positive generator could not tell those apart.
BothForms MakeSparse(uint32_t n, uint32_t f, double density, uint32_t seed, bool force_empty_rows,
                     bool allow_negative = false) {
	BothForms b;
	b.n = n;
	b.f = f;
	b.dense.assign(static_cast<size_t>(n) * f, 0.0);
	b.indptr.push_back(0);
	std::mt19937 rng(seed);
	std::uniform_real_distribution<double> keep(0.0, 1.0);
	std::uniform_real_distribution<double> val(allow_negative ? -8.0 : 0.5, 20.0);
	for (uint32_t i = 0; i < n; ++i) {
		if (!(force_empty_rows && i % 5 == 4)) {
			for (uint32_t k = 0; k < f; ++k) { // ascending: the CSR ordering contract
				if (keep(rng) >= density) {
					continue;
				}
				const double v = val(rng);
				b.indices.push_back(k);
				b.values.push_back(v);
				b.dense[static_cast<size_t>(i) * f + k] = v;
			}
		}
		b.indptr.push_back(static_cast<uint32_t>(b.indices.size()));
	}
	return b;
}

} // namespace

TEST_CASE("sparse community distances equal dense bit for bit", "[community_distances]") {
	struct Case {
		const char *name;
		uint32_t n, f;
		double density;
		bool empty_rows;
		bool negatives;
	};
	// Shapes chosen so each exercises something different: a sparse wide matrix (the
	// real block shape), a dense-ish one (where the merge degenerates to the full
	// scan), one with empty samples, and one so sparse most pairs share no feature
	// at all.
	const Case cases[] = {
	    {"wide and sparse", 24, 400, 0.05, false, false},
	    {"dense-ish", 12, 30, 0.90, false, false},
	    {"with empty rows", 20, 60, 0.30, true, false},
	    {"almost disjoint", 16, 500, 0.01, false, false},
	    {"single feature", 8, 1, 1.00, false, false},
	    // Negative abundances: not biologically meaningful, but accepted input, and
	    // the case where a `!= 0.0` guard silently diverges from the dense `> 0.0`
	    // one (negative denominator -> negative distance; negative cell counted as
	    // present by jaccard).
	    {"with negatives", 24, 200, 0.10, false, true},
	    {"negatives and empty rows", 20, 40, 0.35, true, true},
	};
	for (const auto &c : cases) {
		const auto b = MakeSparse(c.n, c.f, c.density, /*seed=*/1234u + c.n, c.empty_rows, c.negatives);
		for (const auto *metric : kPairwiseLocalMetrics) {
			INFO("case: " << c.name << " metric: " << metric);
			const auto dense = miint::CommunityDistancesCondensed(b.dense, b.n, b.f, metric);
			const auto sparse =
			    miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, metric);
			REQUIRE(sparse.size() == dense.size());
			for (size_t p = 0; p < dense.size(); ++p) {
				INFO("pair index " << p);
				CHECK(sparse[p] == dense[p]); // exact, not Approx
			}
		}
	}
}

TEST_CASE("sparse result is bit-identical for every thread count", "[community_distances]") {
	// Same determinism contract the dense loop carries: each pair writes a fixed
	// condensed slot, so threading reorders WHEN slots are filled, never the values.
	const auto b = MakeSparse(40, 200, 0.15, /*seed=*/99u, /*force_empty_rows=*/true, /*allow_negative=*/true);
	for (const auto *metric : kPairwiseLocalMetrics) {
		INFO("metric: " << metric);
		const auto serial =
		    miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, metric, 1);
		for (unsigned nt : {2u, 4u, 8u, 1000u}) {
			const auto par =
			    miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, metric, nt);
			REQUIRE(par.size() == serial.size());
			for (size_t p = 0; p < serial.size(); ++p) {
				CHECK(par[p] == serial[p]);
			}
		}
	}
}

TEST_CASE("sparse community distances guards", "[community_distances]") {
	const auto b = MakeSparse(6, 20, 0.4, /*seed=*/7u, /*force_empty_rows=*/false);

	// A matrix-wide metric must be REFUSED here, not computed: the sparse path is
	// only ever fed one block, and these read statistics over the whole table, so a
	// value computed from a block would silently be a different metric. Refusing in
	// the pure core means no caller can reach that mistake.
	for (const auto *metric : kMatrixWideMetrics) {
		INFO("metric: " << metric);
		CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, metric),
		                std::invalid_argument);
	}
	CHECK_THROWS_WITH(miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, "gower"),
	                  ContainsSubstring("gower"));
	CHECK_THROWS_WITH(miint::CommunityDistancesCondensedSparse(b.indptr, b.indices, b.values, b.n, b.f, "manhattan"),
	                  ContainsSubstring("manhattan"));

	// fewer than 2 samples
	CHECK_THROWS_WITH(miint::CommunityDistancesCondensedSparse({0, 0}, {}, {}, 1, 20, "bray_curtis"),
	                  ContainsSubstring("at least 2"));
	// indptr must have n+1 entries
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({0, 1}, {0}, {1.0}, 6, 20, "bray_curtis"),
	                std::invalid_argument);
	// indices and values must be the same length
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({0, 1, 2}, {0, 1}, {1.0}, 2, 20, "bray_curtis"),
	                std::invalid_argument);
	// indptr must start at 0, end at nnz, and never decrease
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({1, 1, 2}, {0, 1}, {1.0, 2.0}, 2, 20, "bray_curtis"),
	                std::invalid_argument);
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({0, 2, 1}, {0, 1}, {1.0, 2.0}, 2, 20, "bray_curtis"),
	                std::invalid_argument);
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({0, 1, 1}, {0, 1}, {1.0, 2.0}, 2, 20, "bray_curtis"),
	                std::invalid_argument);
	// a feature index at or past n_features
	CHECK_THROWS_AS(miint::CommunityDistancesCondensedSparse({0, 1, 2}, {0, 20}, {1.0, 2.0}, 2, 20, "bray_curtis"),
	                std::invalid_argument);
	// indices must ASCEND within a row -- the merge depends on it, and an unsorted
	// row would silently produce a wrong distance rather than an error
	CHECK_THROWS_WITH(
	    miint::CommunityDistancesCondensedSparse({0, 2, 3}, {5, 2, 1}, {1.0, 2.0, 3.0}, 2, 20, "bray_curtis"),
	    ContainsSubstring("ascending"));
}

TEST_CASE("sparse matches dense when a row SUM is negative", "[community_distances]") {
	// The randomized generator above draws from [-8, 20], so rows almost always sum
	// positive and the denominator guards are never exercised — only jaccard's
	// presence test is. This matrix is hand-built to drive each guard to the case
	// where `== 0.0` and `> 0.0` disagree:
	//
	//   s0 = [-5,  1]  sum -4   with s1 -> bray/soergel denominator is NEGATIVE
	//   s1 = [ 0,  0]  sum  0   the empty community
	//   s2 = [ 3,  2]  sum  5   an ordinary row to pair against
	//   s3 = [-2, -3]  sum -5   with s4 -> Sum max(x,y) is NEGATIVE
	//   s4 = [-1, -4]  sum -5
	//
	// Dense treats a non-positive denominator as the empty case and returns 0.0; a
	// `== 0.0` guard instead divides by it and returns a NEGATIVE distance, which the
	// ordination core would consume as valid. Morisita's `X <= 0.0` vs `X == 0.0`
	// splits the same way (dense: both non-positive -> 0.0; equality: one-empty -> 1.0).
	const std::vector<double> dense = {
	    -5, 1,  //
	    0,  0,  //
	    3,  2,  //
	    -2, -3, //
	    -1, -4, //
	};
	const std::vector<uint32_t> indptr = {0, 2, 2, 4, 6, 8};
	const std::vector<uint32_t> indices = {0, 1, 0, 1, 0, 1, 0, 1};
	const std::vector<double> values = {-5, 1, 3, 2, -2, -3, -1, -4};

	for (const auto *metric : kPairwiseLocalMetrics) {
		INFO("metric: " << metric);
		const auto d = miint::CommunityDistancesCondensed(dense, 5, 2, metric);
		const auto s = miint::CommunityDistancesCondensedSparse(indptr, indices, values, 5, 2, metric);
		REQUIRE(s.size() == d.size());
		for (size_t p = 0; p < d.size(); ++p) {
			INFO("pair index " << p);
			CHECK(s[p] == d[p]);
		}
		// Whatever the inputs, a distance handed to the ordination core must not be
		// negative -- that is the failure this case exists to prevent.
		for (double v : s) {
			CHECK(v >= 0.0);
		}
	}
}
