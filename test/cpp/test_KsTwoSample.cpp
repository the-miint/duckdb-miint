#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <KsTwoSample.hpp>
#include <alignment_functions_internal.hpp>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <numeric>
#include <random>
#include <utility>
#include <vector>

using namespace miint;

// Two-sided two-sample Kolmogorov-Smirnov test (issue #218). micov's last SciPy
// dependency: every published curve comparison in Weng, Guccione, McDonald et al.
// 2025 is backed by scipy.stats.ks_2samp, so the numbers must not move.
//
// The p-value is derived cleanroom from the lattice-path formulation
// (Hodges 1958, following Drion / Gnedenko-Korolyuk), NOT transcribed from SciPy.
// Under H0 all C(m+n,m) interleavings of the two samples are equally likely. Encode
// one as a monotone lattice path (0,0) -> (m,n): step +x when the next pooled order
// statistic comes from sample A, +y when from B. After i A's and j B's,
//
//     F_a - F_b = i/m - j/n = (i*ng - j*mg) / L
//
// with g = gcd(m,n), mg = m/g, ng = n/g, L = lcm(m,n). So D = max|i*ng - j*mg| / L
// over the path's vertices, an integer h over L, and
//
//     P(D >= h/L) = P(the path leaves the band |i*ng - j*mg| < h)
//
// which this file pins against three independent oracles.

namespace {

// ---------------------------------------------------------------------------
// Oracle A -- the statistic, by literal ECDF evaluation.
// Deliberately the dumbest possible implementation: for every pooled value, count
// how many observations in each sample are <= it. O(N^2) and shares no logic with
// the merge walk under test.
// ---------------------------------------------------------------------------
double BruteForceStatistic(const std::vector<double> &a, const std::vector<double> &b) {
	std::vector<double> pooled;
	pooled.insert(pooled.end(), a.begin(), a.end());
	pooled.insert(pooled.end(), b.begin(), b.end());
	double best = 0.0;
	for (const double v : pooled) {
		size_t ca = 0;
		size_t cb = 0;
		for (const double x : a) {
			if (x <= v) {
				ca++;
			}
		}
		for (const double x : b) {
			if (x <= v) {
				cb++;
			}
		}
		const double diff = std::fabs(static_cast<double>(ca) / static_cast<double>(a.size()) -
		                              static_cast<double>(cb) / static_cast<double>(b.size()));
		best = std::max(best, diff);
	}
	return best;
}

// ---------------------------------------------------------------------------
// Oracle A2 -- the same literal ECDF sweep, but in EXACT integer arithmetic.
// Maximise |ca*nb - cb*na| over the pooled values and divide once at the end. A
// single IEEE division of two exactly-representable integers is correctly rounded,
// so this is the true D to the bit -- which Oracle A is not: two divisions and a
// subtraction per pooled value leave its answer 1-3 ULP off (issue #256).
// ---------------------------------------------------------------------------
double ExactStatistic(const std::vector<double> &a, const std::vector<double> &b) {
	const int64_t na = static_cast<int64_t>(a.size());
	const int64_t nb = static_cast<int64_t>(b.size());
	std::vector<double> pooled;
	pooled.insert(pooled.end(), a.begin(), a.end());
	pooled.insert(pooled.end(), b.begin(), b.end());
	int64_t best = 0;
	for (const double v : pooled) {
		int64_t ca = 0;
		int64_t cb = 0;
		for (const double x : a) {
			if (x <= v) {
				ca++;
			}
		}
		for (const double x : b) {
			if (x <= v) {
				cb++;
			}
		}
		best = std::max<int64_t>(best, std::llabs(ca * nb - cb * na));
	}
	return static_cast<double>(best) / static_cast<double>(na * nb);
}

// ---------------------------------------------------------------------------
// Oracle B -- the exact p-value, by enumerating EVERY interleaving.
// This is the definition itself: walk all C(m+n,m) lattice paths, take each one's
// max|i*ng - j*mg|, and count the fraction reaching h. Ground truth, no recursion
// with the implementation. Feasible only for m+n <= ~20, which is the point.
// ---------------------------------------------------------------------------
double EnumerateExactP(int m, int n, int64_t h) {
	const int64_t g = std::gcd(static_cast<int64_t>(m), static_cast<int64_t>(n));
	const int64_t mg = m / g;
	const int64_t ng = n / g;
	const int total = m + n;
	uint64_t paths = 0;
	uint64_t left_band = 0;
	for (uint64_t mask = 0; mask < (1ull << total); mask++) {
		if (__builtin_popcountll(mask) != static_cast<unsigned>(m)) {
			continue;
		}
		paths++;
		int64_t i = 0;
		int64_t j = 0;
		int64_t best = 0;
		for (int s = 0; s < total; s++) {
			if ((mask >> s) & 1ull) {
				i++;
			} else {
				j++;
			}
			const int64_t offset = i * ng - j * mg;
			best = std::max(best, offset < 0 ? -offset : offset);
		}
		if (best >= h) {
			left_band++;
		}
	}
	return static_cast<double>(left_band) / static_cast<double>(paths);
}

// ---------------------------------------------------------------------------
// Oracle C -- the exact p-value over the FULL rectangle, no banding.
// Same absorbing-boundary probability recursion as the implementation but visiting
// every one of the (m+1)*(n+1) cells, so it cannot share a band-index bug. The
// implementation only visits a window per column; this is what proves the window is
// wide enough.
// ---------------------------------------------------------------------------
double FullRectangleExactP(int64_t m, int64_t n, int64_t h) {
	if (h <= 0) {
		return 1.0;
	}
	const int64_t g = std::gcd(m, n);
	const int64_t mg = m / g;
	const int64_t ng = n / g;
	auto inside = [&](int64_t x, int64_t y) {
		return std::llabs(x * ng - y * mg) < h;
	};

	std::vector<std::vector<double>> f(static_cast<size_t>(m + 1),
	                                   std::vector<double>(static_cast<size_t>(n + 1), 0.0));
	f[0][0] = 1.0; // h >= 1 so (0,0) is inside
	double escape = 0.0;
	for (int64_t x = 0; x <= m; x++) {
		for (int64_t y = 0; y <= n; y++) {
			if (x == 0 && y == 0) {
				continue;
			}
			const double denom = static_cast<double>(m + n - x - y + 1);
			double inflow = 0.0;
			if (x > 0) {
				inflow +=
				    f[static_cast<size_t>(x - 1)][static_cast<size_t>(y)] * static_cast<double>(m - x + 1) / denom;
			}
			if (y > 0) {
				inflow +=
				    f[static_cast<size_t>(x)][static_cast<size_t>(y - 1)] * static_cast<double>(n - y + 1) / denom;
			}
			if (inside(x, y)) {
				f[static_cast<size_t>(x)][static_cast<size_t>(y)] = inflow;
			} else {
				escape += inflow;
			}
		}
	}
	return std::min(std::max(escape, 0.0), 1.0);
}

// 2 / C(m+n, m), as a cancellation-free product: 1/C(m+n,m) = prod_{k=1..m} k/(n+k).
// At D == 1 the only interleavings reaching |i*ng - j*mg| == L are the two fully
// separated orders (all of A then all of B, and the reverse), so this is the exact
// p-value there -- an independent closed form, and it stays accurate where the
// p-value is tiny.
double TwoOverBinom(int64_t m, int64_t n) {
	double p = 2.0;
	for (int64_t k = 1; k <= m; k++) {
		p *= static_cast<double>(k) / static_cast<double>(n + k);
	}
	return p;
}

int64_t Lcm(int64_t a, int64_t b) {
	return a / std::gcd(a, b) * b;
}

} // namespace

// ===========================================================================
// The statistic
// ===========================================================================

TEST_CASE("KsStatistic is zero for a sample against itself", "[ks_2samp]") {
	// REPEATS ARE THE POINT. A walk that compares the ECDFs mid-tie -- after
	// consuming one of the two 1.0s rather than both -- reports 1/5 here instead of
	// 0. Distinct values could not distinguish the two implementations.
	std::vector<double> a {1.0, 1.0, 2.0, 2.0, 3.0};
	std::vector<double> b {1.0, 1.0, 2.0, 2.0, 3.0};
	REQUIRE(KsStatistic(a, b) == 0.0);
}

TEST_CASE("KsStatistic is 1 for disjoint samples", "[ks_2samp]") {
	std::vector<double> a {1.0, 2.0, 3.0};
	std::vector<double> b {10.0, 20.0};
	REQUIRE(KsStatistic(a, b) == 1.0);
	// Symmetric in its arguments: b entirely BELOW a exercises the other branch of
	// the min() in the merge walk.
	std::vector<double> c {10.0, 20.0};
	std::vector<double> d {1.0, 2.0, 3.0};
	REQUIRE(KsStatistic(c, d) == 1.0);
}

TEST_CASE("KsStatistic on a hand-computed case", "[ks_2samp]") {
	// a = {1,2,3,4}, b = {2.5}.
	//   v=1.0: F_a=1/4, F_b=0   -> 0.25
	//   v=2.0: F_a=2/4, F_b=0   -> 0.50
	//   v=2.5: F_a=2/4, F_b=1   -> 0.50
	//   v=3.0: F_a=3/4, F_b=1   -> 0.25
	//   v=4.0: F_a=1,   F_b=1   -> 0.00
	std::vector<double> a {4.0, 1.0, 3.0, 2.0}; // unsorted on purpose
	std::vector<double> b {2.5};
	REQUIRE(KsStatistic(a, b) == 0.5);
}

TEST_CASE("KsStatistic handles ties across the two samples", "[ks_2samp]") {
	// Every value is shared; only the multiplicities differ.
	//   v=1: F_a=2/3, F_b=1/3 -> 1/3
	//   v=2: F_a=1,   F_b=1   -> 0
	std::vector<double> a {1.0, 1.0, 2.0};
	std::vector<double> b {1.0, 2.0, 2.0};
	REQUIRE_THAT(KsStatistic(a, b), Catch::Matchers::WithinAbs(1.0 / 3.0, 1e-15));
}

TEST_CASE("KsStatistic with one observation per sample", "[ks_2samp]") {
	std::vector<double> a {1.0};
	std::vector<double> b {2.0};
	REQUIRE(KsStatistic(a, b) == 1.0);
	std::vector<double> c {5.0};
	std::vector<double> d {5.0};
	REQUIRE(KsStatistic(c, d) == 0.0);
}

TEST_CASE("KsStatistic tolerates infinities", "[ks_2samp]") {
	// +/-inf are legitimate order statistics: they sort and compare, so the ECDF is
	// still well defined. Only NaN is rejected (it would break the sort).
	const double inf = std::numeric_limits<double>::infinity();
	std::vector<double> a {1.0, inf};
	std::vector<double> b {2.0, 3.0};
	REQUIRE(KsStatistic(a, b) == 0.5);
	std::vector<double> c {-inf, 0.0};
	std::vector<double> d {-inf, 0.0};
	REQUIRE(KsStatistic(c, d) == 0.0);
}

TEST_CASE("KsStatistic matches the literal ECDF oracle on random inputs", "[ks_2samp]") {
	std::mt19937_64 rng(20260817);
	std::uniform_int_distribution<int> size_dist(1, 12);
	int trials = 0;
	for (int t = 0; t < 500; t++) {
		const int na = size_dist(rng);
		const int nb = size_dist(rng);
		// Draw from a SMALL integer alphabet so ties are frequent -- continuous draws
		// would never exercise the tie handling, which is the part that is easy to
		// get wrong.
		std::uniform_int_distribution<int> val_dist(0, 4);
		std::vector<double> a(static_cast<size_t>(na));
		std::vector<double> b(static_cast<size_t>(nb));
		for (auto &v : a) {
			v = static_cast<double>(val_dist(rng));
		}
		for (auto &v : b) {
			v = static_cast<double>(val_dist(rng));
		}
		const double expected = BruteForceStatistic(a, b);
		const double got = KsStatistic(a, b);
		REQUIRE_THAT(got, Catch::Matchers::WithinAbs(expected, 1e-15));
		trials++;
	}
	REQUIRE(trials == 500);
}

TEST_CASE("KsStatistic returns the exact D, not the raw sweep value", "[ks_2samp]") {
	// The merge walk accumulates i/na - j/nb in floating point -- two divisions and a
	// subtraction, three rounding steps -- so its running maximum is not in general the
	// correctly-rounded D. The smallest case where it is not: na=1, nb=3 with the pooled
	// order b < a < b < b. The gap peaks at (i,j) = (1,1), and D = 2/3 (issue #256).
	//
	// The first REQUIRE is load-bearing, not decoration: it is the pre-snap value, and
	// if it ever coincided with 2.0/3.0 the assertion below would pin nothing.
	REQUIRE(1.0 / 1.0 - 1.0 / 3.0 != 2.0 / 3.0);
	std::vector<double> a {2.0};
	std::vector<double> b {1.0, 3.0, 4.0};
	REQUIRE(KsStatistic(a, b) == 2.0 / 3.0);
}

TEST_CASE("KsStatistic is bit-exact against the integer-arithmetic oracle", "[ks_2samp]") {
	// Every achievable D is an exact multiple of 1/lcm(na,nb), so "close enough" is the
	// wrong bar here: there is one right double and the function should return it. Sizes
	// reach 30 because micov's n is the number of samples in a metadata group, 10-100.
	std::mt19937_64 rng(20260825);
	std::uniform_int_distribution<int> size_dist(1, 30);
	int inexact = 0;
	for (int t = 0; t < 500; t++) {
		const int na = size_dist(rng);
		const int nb = size_dist(rng);
		// Small integer alphabet, as above: cumulative coverage curves are nothing but
		// ties, and the maximum then lands on a lattice point reachable in several ways.
		std::uniform_int_distribution<int> val_dist(0, 8);
		std::vector<double> a(static_cast<size_t>(na));
		std::vector<double> b(static_cast<size_t>(nb));
		for (auto &v : a) {
			v = static_cast<double>(val_dist(rng));
		}
		for (auto &v : b) {
			v = static_cast<double>(val_dist(rng));
		}
		const double expected = ExactStatistic(a, b);
		// Oracle A evaluates the same expression the merge walk does, in the same order,
		// so it reproduces the UNSNAPPED answer. Counting where it disagrees proves this
		// grid actually reaches the cases that distinguish the two (issue #256) -- an
		// alphabet or a size range that never did would make the check below vacuous.
		if (BruteForceStatistic(a, b) != expected) {
			inexact++;
		}
		REQUIRE(KsStatistic(a, b) == expected);
	}
	REQUIRE(inexact > 50);
}

// ===========================================================================
// The exact p-value
// ===========================================================================

TEST_CASE("KsExactPValue is 1 when the statistic is 0", "[ks_2samp]") {
	REQUIRE(KsExactPValue(5, 5, 0.0) == 1.0);
	REQUIRE(KsExactPValue(7, 13, 0.0) == 1.0);
}

TEST_CASE("KsExactPValue matches enumeration of every interleaving", "[ks_2samp]") {
	// Ground truth straight from the definition. Covers equal sizes, coprime sizes
	// (g == 1, so the band has a non-integer slope), and a common factor (g > 1),
	// at every achievable h.
	int checked = 0;
	for (const auto &mn : std::vector<std::pair<int, int>> {{1, 1},
	                                                        {1, 4},
	                                                        {2, 2},
	                                                        {2, 3},
	                                                        {3, 3},
	                                                        {3, 4},
	                                                        {4, 6},
	                                                        {5, 5},
	                                                        {3, 7},
	                                                        {6, 9},
	                                                        {4, 4},
	                                                        {2, 8},
	                                                        {7, 7},
	                                                        {5, 9}}) {
		const int m = mn.first;
		const int n = mn.second;
		const int64_t lcm = Lcm(m, n);
		for (int64_t h = 1; h <= lcm; h++) {
			const double d = static_cast<double>(h) / static_cast<double>(lcm);
			const double expected = EnumerateExactP(m, n, h);
			const double got = KsExactPValue(m, n, d);
			REQUIRE_THAT(got,
			             Catch::Matchers::WithinRel(expected, 1e-12) || Catch::Matchers::WithinAbs(expected, 1e-300));
			checked++;
		}
	}
	// Exactly one check per achievable h, so the total is the sum of the lcms above:
	// 1+4+2+6+3+12+12+5+21+18+4+8+7+45. Pinned exactly rather than as a floor, so
	// that editing the grid has to be deliberate instead of quietly losing coverage.
	REQUIRE(checked == 148);
}

TEST_CASE("KsExactPValue banded window agrees with the full rectangle", "[ks_2samp]") {
	// The implementation visits only a window of rows per column. A window that is
	// one row too narrow loses escaping probability mass, and a p-value that is too
	// small is exactly the failure a user would never notice. The input class that
	// breaks a too-narrow window is a LOPSIDED pair (m << n), where the band spans
	// many rows per column, plus coprime sizes so the band edges land between rows.
	int checked = 0;
	for (const auto &mn : std::vector<std::pair<int64_t, int64_t>> {{2, 97},
	                                                                {3, 100},
	                                                                {97, 2},
	                                                                {100, 3},
	                                                                {11, 89},
	                                                                {89, 11},
	                                                                {50, 50},
	                                                                {49, 51},
	                                                                {64, 48},
	                                                                {37, 91},
	                                                                {1, 60}}) {
		const int64_t m = mn.first;
		const int64_t n = mn.second;
		const int64_t lcm = Lcm(m, n);
		for (const double frac : {0.01, 0.05, 0.13, 0.25, 0.5, 0.77, 1.0}) {
			const int64_t h = static_cast<int64_t>(std::llround(frac * static_cast<double>(lcm)));
			if (h <= 0) {
				continue;
			}
			const double d = static_cast<double>(h) / static_cast<double>(lcm);
			const double expected = FullRectangleExactP(m, n, h);
			const double got = KsExactPValue(m, n, d);
			REQUIRE_THAT(got,
			             Catch::Matchers::WithinRel(expected, 1e-12) || Catch::Matchers::WithinAbs(expected, 1e-300));
			checked++;
		}
	}
	// Exactly 11 (m,n) pairs x 7 fractions, with nothing skipped: the smallest lcm in the
	// list is 50, so even frac 0.01 rounds to h >= 1. Pinned exactly rather than as a
	// floor -- a floor of 60 would stay green after deleting all four lopsided pairs,
	// which the comment above calls the only input class that breaks a narrow window.
	REQUIRE(checked == 77);
}

TEST_CASE("KsExactPValue at D == 1 equals the closed form 2/C(m+n,m)", "[ks_2samp]") {
	// An independent closed form, and the one that pins PRECISION: at n=m=200 the
	// answer is 1.9e-119. Computing the p-value as 1 - P(stay inside) would cancel
	// away every significant digit here (measured: 4.2e-06 relative error at
	// n=m=20), which is why the implementation accumulates escaping mass instead.
	for (const auto &mn : std::vector<std::pair<int64_t, int64_t>> {
	         {1, 1}, {2, 3}, {5, 5}, {7, 13}, {20, 20}, {50, 50}, {64, 48}, {100, 100}, {200, 200}, {37, 91}}) {
		const int64_t m = mn.first;
		const int64_t n = mn.second;
		const double expected = TwoOverBinom(m, n);
		const double got = KsExactPValue(m, n, 1.0);
		REQUIRE_THAT(got, Catch::Matchers::WithinRel(expected, 1e-12));
	}
}

TEST_CASE("KsExactPValue clamps a rounding overshoot to 1", "[ks_2samp]") {
	// The escape mass is a sum of positive terms, so rounding can carry it a step
	// above 1 and the clamp is what turns it back into a probability.
	//
	// These two cases are where that actually happens: n1=7, n2=13 at h=9/91 and
	// h=10/91 accumulate to 1.0000000000000002. Found by replicating the recursion's
	// exact operation order in double precision and searching -- 684 of the
	// (n1,n2,h) combinations in 1..40 x 1..40 overshoot, so it is ordinary rather
	// than exotic. An earlier version of this test used the two cases below instead,
	// where the raw sum is already exactly 1.0, and mutation testing showed it could
	// not fail: removing the clamp left it green.
	REQUIRE(KsExactPValue(7, 13, 9.0 / 91.0) == 1.0);
	REQUIRE(KsExactPValue(7, 13, 10.0 / 91.0) == 1.0);

	// n1 = n2 = 5 with D = 0.2 does NOT exercise the clamp -- its raw sum is exactly
	// 1.0. It is kept for a different reason: this is the case where SciPy's OWN
	// exact routine returns 1.0000000000000002, declares itself unsuccessful, and
	// silently falls back to its asymptotic branch before clipping to 1.0. We land on
	// the same answer without needing that fallback.
	REQUIRE(KsExactPValue(5, 5, 0.2) == 1.0);
	// h == 1 admits only vertices with i*ng == j*mg exactly. For equal sizes that is
	// i == j, which no path satisfies after its first step, so every path escapes and
	// the p-value is exactly 1.
	REQUIRE(KsExactPValue(9, 9, 1.0 / 9.0) == 1.0);
}

TEST_CASE("KsExactPValue snaps an unachievable d to the nearest achievable one", "[ks_2samp]") {
	// D is always an exact multiple of 1/lcm(n1,n2), so KsStatistic can never produce
	// anything else -- but KsExactPValue is public and documents that it ROUNDS rather
	// than rejects. At n1 = n2 = 5, lcm is 5 and the achievable statistics are k/5, so
	// 0.3 is not one of them: llround(0.3 * 5) = 2, which is D = 0.4.
	//
	// Without this, replacing llround with a truncating cast would leave the whole suite
	// green -- 0.3 would silently become D = 0.2 (p = 1.0) instead of D = 0.4
	// (p = 0.873) -- because every other caller passes an already-achievable d.
	REQUIRE(KsExactPValue(5, 5, 0.3) == KsExactPValue(5, 5, 0.4));
	REQUIRE(KsExactPValue(5, 5, 0.3) != KsExactPValue(5, 5, 0.2));
	// Rounds down as well as up: 0.29 is nearer 0.2.
	REQUIRE(KsExactPValue(5, 5, 0.29) == KsExactPValue(5, 5, 0.2));
}

TEST_CASE("KsExactPValue is symmetric in its sample sizes, to the bit", "[ks_2samp]") {
	// D is symmetric, so the p-value must be too. Asserted as EXACT equality, which
	// holds only because the implementation canonicalises to n1 >= n2 before running
	// the recursion: without that, transposing reorders the arithmetic and the two
	// answers differ by 1 ULP (measured 0.2311145510835913 vs 0.23111455108359133 at
	// 7 vs 13). So this pins a deliberate design decision, not a tautology.
	for (const auto &mn : std::vector<std::pair<int64_t, int64_t>> {{3, 11}, {7, 13}, {12, 30}, {37, 91}, {2, 97}}) {
		const int64_t lcm = Lcm(mn.first, mn.second);
		for (const int64_t h : {static_cast<int64_t>(1), static_cast<int64_t>(3), lcm / 4, lcm / 2, lcm}) {
			if (h <= 0) {
				continue;
			}
			const double d = static_cast<double>(h) / static_cast<double>(lcm);
			REQUIRE(KsExactPValue(mn.first, mn.second, d) == KsExactPValue(mn.second, mn.first, d));
		}
	}
}

TEST_CASE("KsExactPValue is non-increasing in the statistic", "[ks_2samp]") {
	// A larger observed difference cannot be more probable under H0. This is NOT
	// structurally guaranteed by the recursion -- widening the band changes which
	// cells absorb, so a band-edge sign error breaks monotonicity while leaving the
	// endpoints (h=1 -> 1.0, h=lcm -> 2/C) correct.
	for (const auto &mn : std::vector<std::pair<int64_t, int64_t>> {{10, 10}, {7, 13}, {12, 30}, {5, 9}}) {
		const int64_t lcm = Lcm(mn.first, mn.second);
		double previous = 1.0;
		for (int64_t h = 1; h <= lcm; h++) {
			const double d = static_cast<double>(h) / static_cast<double>(lcm);
			const double p = KsExactPValue(mn.first, mn.second, d);
			REQUIRE(p <= previous + 1e-15);
			REQUIRE(p >= 0.0);
			REQUIRE(p <= 1.0);
			previous = p;
		}
	}
}

// ===========================================================================
// The whole test, and its input validation
// ===========================================================================

TEST_CASE("KsTwoSample combines the statistic and the exact p-value", "[ks_2samp]") {
	const auto r = KsTwoSample({1.0, 2.0, 3.0, 4.0, 5.0}, {6.0, 7.0, 8.0, 9.0, 10.0});
	REQUIRE(r.statistic == 1.0);
	// 2/C(10,5) = 2/252
	REQUIRE_THAT(r.pvalue, Catch::Matchers::WithinRel(2.0 / 252.0, 1e-13));

	const auto same = KsTwoSample({1.0, 2.0, 3.0}, {1.0, 2.0, 3.0});
	REQUIRE(same.statistic == 0.0);
	REQUIRE(same.pvalue == 1.0);
}

TEST_CASE("KsTwoSample rejects an empty sample", "[ks_2samp]") {
	REQUIRE_THROWS_AS(KsTwoSample({}, {1.0, 2.0}), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(KsTwoSample({}, {1.0, 2.0}), Catch::Matchers::ContainsSubstring("must not be empty"));
	REQUIRE_THROWS_WITH(KsTwoSample({1.0, 2.0}, {}), Catch::Matchers::ContainsSubstring("must not be empty"));
}

TEST_CASE("KsTwoSample rejects NaN", "[ks_2samp]") {
	const double nan = std::numeric_limits<double>::quiet_NaN();
	REQUIRE_THROWS_AS(KsTwoSample({1.0, nan}, {2.0, 3.0}), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(KsTwoSample({1.0, nan}, {2.0, 3.0}), Catch::Matchers::ContainsSubstring("NaN"));
	// Both argument positions are checked, not just the first.
	REQUIRE_THROWS_WITH(KsTwoSample({1.0, 2.0}, {nan}), Catch::Matchers::ContainsSubstring("NaN"));
}

TEST_CASE("KsTwoSample rejects a sample above the exact ceiling", "[ks_2samp]") {
	// The ceiling applies to max(n1, n2), so an oversized sample must be rejected
	// even when the OTHER sample is tiny -- that is the case a `n1 > cap && n2 > cap`
	// check would wrongly admit.
	std::vector<double> big(static_cast<size_t>(KS_MAX_EXACT_N) + 1, 1.0);
	REQUIRE_THROWS_AS(KsTwoSample(big, {1.0, 2.0}), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(KsTwoSample(big, {1.0, 2.0}), Catch::Matchers::ContainsSubstring("10000"));
	// The message has to say what to do about it, not just that it failed.
	REQUIRE_THROWS_WITH(KsTwoSample(big, {1.0, 2.0}), Catch::Matchers::ContainsSubstring("asymp"));
	REQUIRE_THROWS_WITH(KsTwoSample({1.0, 2.0}, big), Catch::Matchers::ContainsSubstring("10000"));

	// Exactly at the ceiling is accepted -- pins the boundary against an off-by-one.
	std::vector<double> at_cap(static_cast<size_t>(KS_MAX_EXACT_N), 1.0);
	const auto r = KsTwoSample(at_cap, at_cap);
	REQUIRE(r.statistic == 0.0);
	REQUIRE(r.pvalue == 1.0);
}

TEST_CASE("ValidateKsMethod accepts auto and exact, case-insensitively", "[ks_2samp]") {
	REQUIRE_NOTHROW(ValidateKsMethod("auto"));
	REQUIRE_NOTHROW(ValidateKsMethod("exact"));
	REQUIRE_NOTHROW(ValidateKsMethod("AUTO"));
	REQUIRE_NOTHROW(ValidateKsMethod("Exact"));
}

TEST_CASE("ValidateKsMethod rejects asymp with its own explanation", "[ks_2samp]") {
	// 'asymp' is a legitimate SciPy method that miint does not implement, so it must
	// not be reported as though the caller made a typo.
	REQUIRE_THROWS_AS(ValidateKsMethod("asymp"), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(ValidateKsMethod("asymp"), Catch::Matchers::ContainsSubstring("not implemented"));
	REQUIRE_THROWS_WITH(ValidateKsMethod("asymp"), Catch::Matchers::ContainsSubstring("218"));
}

TEST_CASE("ValidateKsMethod rejects an unknown method", "[ks_2samp]") {
	REQUIRE_THROWS_AS(ValidateKsMethod("bogus"), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(ValidateKsMethod("bogus"), Catch::Matchers::ContainsSubstring("unknown method"));
	// Names the accepted values so the caller can fix it without reading the docs.
	REQUIRE_THROWS_WITH(ValidateKsMethod("bogus"), Catch::Matchers::ContainsSubstring("'auto'"));
	REQUIRE_THROWS_WITH(ValidateKsMethod(""), Catch::Matchers::ContainsSubstring("unknown method"));
}

TEST_CASE("KsTwoSample end to end against enumeration on random samples", "[ks_2samp]") {
	// Ties everywhere (small integer alphabet), unequal sizes, and the p-value
	// checked against literal path enumeration -- so statistic and p-value are
	// verified together rather than each against its own convenient fixture.
	std::mt19937_64 rng(987654321);
	std::uniform_int_distribution<int> size_dist(1, 8);
	std::uniform_int_distribution<int> val_dist(0, 3);
	int trials = 0;
	for (int t = 0; t < 200; t++) {
		const int na = size_dist(rng);
		const int nb = size_dist(rng);
		std::vector<double> a(static_cast<size_t>(na));
		std::vector<double> b(static_cast<size_t>(nb));
		for (auto &v : a) {
			v = static_cast<double>(val_dist(rng));
		}
		for (auto &v : b) {
			v = static_cast<double>(val_dist(rng));
		}
		const auto r = KsTwoSample(a, b);
		const int64_t lcm = Lcm(na, nb);
		const int64_t h = static_cast<int64_t>(std::llround(r.statistic * static_cast<double>(lcm)));
		const double expected = h <= 0 ? 1.0 : EnumerateExactP(na, nb, h);
		REQUIRE_THAT(r.pvalue,
		             Catch::Matchers::WithinRel(expected, 1e-12) || Catch::Matchers::WithinAbs(expected, 1e-300));
		trials++;
	}
	REQUIRE(trials == 200);
}
