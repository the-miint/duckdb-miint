#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <cstdint>
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

TEST_CASE("pearson constant profile -> NaN (undefined correlation)", "[community_distances]") {
	// A constant row has zero variance, so Pearson r is undefined (scipy returns
	// nan). We surface NaN so a downstream reader rejects it as "not provided".
	const std::vector<double> m = {2, 2, 2, 2, 1, 3, 0, 5};
	auto d = miint::CommunityDistancesCondensed(m, 2, 4, "pearson");
	REQUIRE(d.size() == 1);
	CHECK(std::isnan(d[0]));
	// A non-constant identical pair is well-defined (r=1 -> distance 0), i.e. NaN
	// is specific to the zero-variance case, not to pearson generally.
	const std::vector<double> ok = {1, 3, 0, 5, 1, 3, 0, 5};
	CHECK(miint::CommunityDistancesCondensed(ok, 2, 4, "pearson")[0] == Approx(0.0).margin(1e-12));
}

TEST_CASE("chisq zero-row-sum sample -> NaN (no row profile)", "[community_distances]") {
	// An all-zero sample has no correspondence-analysis row profile; the distance
	// is undefined and surfaced as NaN. (Reachable only from a dense caller: the
	// sparse SQL reader drops all-zero samples before they reach the core.)
	const std::vector<double> m = {0, 0, 0, 0, 1, 2, 3, 4};
	auto d = miint::CommunityDistancesCondensed(m, 2, 4, "chisq");
	REQUIRE(d.size() == 1);
	CHECK(std::isnan(d[0]));
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
