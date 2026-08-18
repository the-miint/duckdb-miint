#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <CumulativeCoverageAccumulator.hpp>
#include <alignment_functions_internal.hpp>

#include <algorithm>
#include <random>
#include <set>
#include <vector>

using namespace miint;

// cumulative_coverage answers: for each rank k, how many bases are covered by the
// UNION of every interval with rank <= k (issue #214). micov recompresses the whole
// accumulated set at every rank, which is O(n^2) in intervals and dominates a micov
// run. This computes all ranks in one O(n log n) sweep.
//
// The trick is the min-covering-rank identity: the union of ranks 0..k is exactly the
// set of bases whose LOWEST covering rank is <= k. So histogram bases by their minimum
// covering rank, and the curve is that histogram's prefix sum.
//
// Intervals are 1-based half-open, matching read_alignments and compress_intervals.

namespace {

// Literal set-union oracle. Deliberately the dumbest possible implementation --
// materializes every covered base into a std::set -- so it shares no logic with the
// sweep under test. Only usable on small coordinate ranges, which is the point.
std::vector<int64_t> BruteForceCurve(const std::vector<std::tuple<int32_t, int64_t, int64_t>> &obs, int32_t n_ranks) {
	std::vector<int64_t> curve;
	std::set<int64_t> accumulated;
	for (int32_t k = 0; k < n_ranks; k++) {
		for (const auto &o : obs) {
			if (std::get<0>(o) != k) {
				continue;
			}
			for (int64_t p = std::get<1>(o); p < std::get<2>(o); p++) {
				accumulated.insert(p);
			}
		}
		curve.push_back(static_cast<int64_t>(accumulated.size()));
	}
	return curve;
}

std::vector<int64_t> CoveredOf(const std::vector<CumulativePoint> &pts) {
	std::vector<int64_t> out;
	out.reserve(pts.size());
	for (const auto &p : pts) {
		out.push_back(p.covered);
	}
	return out;
}

} // namespace

TEST_CASE("CumulativeCoverageAccumulator - empty state yields no curve", "[cumulative_coverage]") {
	// Fails if Curve() fabricates a rank-0 point for a group nothing was added to.
	CumulativeCoverageAccumulator cc;
	REQUIRE(cc.Empty());
	REQUIRE(cc.Curve().empty());
}

TEST_CASE("CumulativeCoverageAccumulator - single rank single interval", "[cumulative_coverage]") {
	// Half-open width: [10,20) is 10 bases, not 11. Fails on a +1.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);

	auto curve = cc.Curve();
	REQUIRE(curve.size() == 1);
	REQUIRE(curve[0].rank == 0);
	REQUIRE(curve[0].covered == 10);
}

TEST_CASE("CumulativeCoverageAccumulator - disjoint ranks accumulate", "[cumulative_coverage]") {
	// Fails if later ranks replace rather than accumulate (would give 10, 10).
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(1, 100, 110);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {10, 20});
}

TEST_CASE("CumulativeCoverageAccumulator - a nested later interval adds nothing", "[cumulative_coverage]") {
	// rank 1 sits entirely inside rank 0, so it contributes no NEW base.
	// Fails if the union is computed as a sum of breadths (would give 20, 25).
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 30);
	cc.Add(1, 15, 20);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {20, 20});
}

TEST_CASE("CumulativeCoverageAccumulator - an enclosing later interval adds only the difference",
          "[cumulative_coverage]") {
	// The reverse nesting, which is where min-covering-rank earns its name: bases
	// 15..19 have minimum rank 0 even though rank 1 also covers them.
	// Fails if the sweep credits a segment to the LAST opener rather than the lowest.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 15, 20);
	cc.Add(1, 10, 30);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {5, 20});
}

TEST_CASE("CumulativeCoverageAccumulator - touching intervals neither overlap nor gap", "[cumulative_coverage]") {
	// [10,20) and [20,30) share no base and leave none uncovered: 10 then 20.
	// Fails at 19 (treating the boundary as shared) or 21 (inventing a base).
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(1, 20, 30);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {10, 20});
}

TEST_CASE("CumulativeCoverageAccumulator - identical intervals are not double counted", "[cumulative_coverage]") {
	// Two samples covering exactly the same span. Fails at {10, 20} if the second
	// rank's bases are added instead of unioned.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(1, 10, 20);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {10, 10});
}

TEST_CASE("CumulativeCoverageAccumulator - overlap WITHIN one rank counts once", "[cumulative_coverage]") {
	// One sample's own intervals overlapping each other. [10,20) u [15,25) = 15.
	// Fails at 20 if per-rank intervals are summed.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(0, 15, 25);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {15});
}

TEST_CASE("CumulativeCoverageAccumulator - a zero-coverage rank holds its position", "[cumulative_coverage]") {
	// The whole reason AddEmptyRank exists. A sample in the group with no coverage of
	// the target must still occupy a rank, or the group size is understated and curves
	// from groups with different detection rates stop being comparable.
	// Fails at {10} -- a one-point curve -- if the empty rank is dropped.
	CumulativeCoverageAccumulator cc;
	cc.AddEmptyRank(0);
	cc.Add(1, 10, 20);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {0, 10});
}

TEST_CASE("CumulativeCoverageAccumulator - a group with no coverage at all is a flat curve", "[cumulative_coverage]") {
	// Not an empty result: three samples were observed, all with zero coverage.
	// Fails if this collapses to an empty curve, which would silently drop the group.
	CumulativeCoverageAccumulator cc;
	cc.AddEmptyRank(0);
	cc.AddEmptyRank(1);
	cc.AddEmptyRank(2);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {0, 0, 0});
}

TEST_CASE("CumulativeCoverageAccumulator - degenerate intervals cover nothing but keep the rank",
          "[cumulative_coverage]") {
	// Zero-width and inverted intervals cover no base. Note this DIVERGES from
	// IntervalCompressor, which silently swaps an inverted pair -- there [20,10)
	// becomes 10 covered bases. Swapping invents coverage from transposed columns, so
	// the newer convention (region_presence, region_coverage) drops instead.
	// Fails at 10 if inverted intervals are swapped rather than dropped.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 15, 15);
	cc.Add(0, 20, 10);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {0});
}

// The three rejection cases below assert on the MESSAGE, not just the exception type,
// because they are three distinct causes with three distinct remedies and each has its
// own branch. Type-only assertions let any one branch absorb all three -- which is
// exactly what happened here: an earlier implementation had a pre-allocation size
// guard that fired first on every malformed input, leaving the gap and start-at-zero
// branches unreachable. Both tests still passed, because both still threw.
//
// Each fixture also carries enough well-formed rows that a "highest rank exceeds the
// row count" shortcut cannot be what rejects it. That is the input class those earlier
// tests were missing.

TEST_CASE("CumulativeCoverageAccumulator - a gap in the ranks is rejected", "[cumulative_coverage]") {
	// ROW_NUMBER() - 1 always produces 0..R-1. A gap means the caller built ranks some
	// other way, and the curve's x-axis would be a lie -- rank 1 would report what is
	// really rank 2's coverage. Four rows, highest rank 2, so row-count reasoning
	// alone cannot detect it.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(0, 25, 30);
	cc.Add(0, 35, 40);
	cc.Add(2, 50, 60);

	REQUIRE_THROWS_AS(cc.Curve(), miint::InvalidInputException);
	REQUIRE_THROWS_WITH(cc.Curve(), Catch::Matchers::ContainsSubstring("rank 1 is missing"));
}

TEST_CASE("CumulativeCoverageAccumulator - ranks must start at zero", "[cumulative_coverage]") {
	// The likely way to get here is ROW_NUMBER() without the - 1, so the message says
	// so. A curve whose first point is rank 1 has lost its base and every later point
	// is shifted. Four rows, highest rank 2.
	CumulativeCoverageAccumulator cc;
	cc.Add(1, 10, 20);
	cc.Add(1, 25, 30);
	cc.Add(1, 35, 40);
	cc.Add(2, 50, 60);

	REQUIRE_THROWS_WITH(cc.Curve(), Catch::Matchers::ContainsSubstring("ranks must start at 0"));
	REQUIRE_THROWS_WITH(cc.Curve(), Catch::Matchers::ContainsSubstring("note the - 1"));
}

TEST_CASE("CumulativeCoverageAccumulator - negative ranks are rejected", "[cumulative_coverage]") {
	// A negative rank would index the histogram out of bounds, so this must be caught
	// before anything is sized. Padded with valid ranks so the row count is ample.
	CumulativeCoverageAccumulator cc;
	cc.Add(-1, 10, 20);
	cc.Add(0, 25, 30);
	cc.Add(1, 35, 40);
	cc.Add(2, 50, 60);

	REQUIRE_THROWS_WITH(cc.Curve(), Catch::Matchers::ContainsSubstring("ranks must start at 0"));
}

TEST_CASE("CumulativeCoverageAccumulator - an implausibly high rank is a gap, not an allocation",
          "[cumulative_coverage]") {
	// Rejected as a gap before the histogram is sized, so this must not try to
	// allocate a billion buckets. Fails by exhausting memory rather than by returning
	// a wrong answer, which is why it is asserted explicitly.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(1'000'000'000, 30, 40);

	REQUIRE_THROWS_WITH(cc.Curve(), Catch::Matchers::ContainsSubstring("rank 1 is missing"));
}

TEST_CASE("CumulativeCoverageAccumulator - Absorb matches a single accumulated state", "[cumulative_coverage]") {
	// This is the DuckDB Combine path: parallel GROUP BY builds one state per thread.
	// Fails if Absorb forgets empty ranks -- a case a two-interval test would miss,
	// since the curve length would still look right.
	CumulativeCoverageAccumulator split_a;
	split_a.Add(0, 10, 20);
	split_a.AddEmptyRank(2);

	CumulativeCoverageAccumulator split_b;
	split_b.Add(1, 15, 30);
	split_b.Add(0, 50, 60);

	split_a.Absorb(split_b);

	CumulativeCoverageAccumulator single;
	single.Add(0, 10, 20);
	single.AddEmptyRank(2);
	single.Add(1, 15, 30);
	single.Add(0, 50, 60);

	REQUIRE(CoveredOf(split_a.Curve()) == CoveredOf(single.Curve()));
	REQUIRE(split_a.Curve().size() == 3);
}

TEST_CASE("CumulativeCoverageAccumulator - curve is independent of insertion order", "[cumulative_coverage]") {
	// A GROUP BY delivers rows in whatever order the scan produces. Fails if the sweep
	// leans on inputs arriving sorted by position or by rank.
	std::vector<std::tuple<int32_t, int64_t, int64_t>> obs = {{2, 40, 90},  {0, 10, 20}, {1, 15, 45},
	                                                          {0, 80, 100}, {2, 5, 12},  {1, 60, 70}};

	CumulativeCoverageAccumulator forward;
	for (const auto &o : obs) {
		forward.Add(std::get<0>(o), std::get<1>(o), std::get<2>(o));
	}

	CumulativeCoverageAccumulator reverse;
	for (auto it = obs.rbegin(); it != obs.rend(); ++it) {
		reverse.Add(std::get<0>(*it), std::get<1>(*it), std::get<2>(*it));
	}

	REQUIRE(CoveredOf(forward.Curve()) == CoveredOf(reverse.Curve()));
	REQUIRE(CoveredOf(forward.Curve()) == BruteForceCurve(obs, 3));
}

TEST_CASE("CumulativeCoverageAccumulator - matches a literal set-union oracle on random input",
          "[cumulative_coverage]") {
	// The assertion that would catch a logic error the hand-built cases above miss.
	// Small coordinate range on purpose: it forces dense overlap between ranks, which
	// is where a min-covering-rank bug shows up. Fixed seed so a failure reproduces.
	std::mt19937 rng(20260817);
	std::uniform_int_distribution<int32_t> n_ranks_dist(1, 8);
	std::uniform_int_distribution<int> n_obs_dist(0, 6);
	std::uniform_int_distribution<int64_t> pos_dist(1, 60);

	for (int trial = 0; trial < 400; trial++) {
		int32_t n_ranks = n_ranks_dist(rng);

		std::vector<std::tuple<int32_t, int64_t, int64_t>> obs;
		CumulativeCoverageAccumulator cc;
		for (int32_t r = 0; r < n_ranks; r++) {
			int n_obs = n_obs_dist(rng);
			// Every rank is registered even when it draws no intervals, so
			// zero-coverage samples are exercised throughout rather than only at the
			// edges.
			cc.AddEmptyRank(r);
			for (int j = 0; j < n_obs; j++) {
				int64_t a = pos_dist(rng);
				int64_t b = pos_dist(rng);
				if (a == b) {
					continue;
				}
				int64_t lo = std::min(a, b);
				int64_t hi = std::max(a, b);
				cc.Add(r, lo, hi);
				obs.emplace_back(r, lo, hi);
			}
		}

		auto expected = BruteForceCurve(obs, n_ranks);
		auto actual = CoveredOf(cc.Curve());
		REQUIRE(actual == expected);
	}
}

TEST_CASE("CumulativeCoverageAccumulator - sparse wide-coordinate curves match the oracle", "[cumulative_coverage]") {
	// Deliberately NOT a monotonicity check. The curve is a prefix sum of per-rank widths
	// that are non-negative by construction, so it is monotonic for ANY choice of
	// credited rank -- crediting every segment to the HIGHEST active rank would satisfy
	// it too. A monotonicity assertion here cannot fail, and the version of this case
	// that made one claimed to catch precisely the misattribution it is blind to.
	//
	// So this compares against the oracle instead, on a wider coordinate range than the
	// dense case above. Sparse intervals mean most ranks contribute disjoint runs, which
	// exercises the heap's lazy cleanup across long gaps rather than at coincident
	// endpoints. Curve length and rank labelling are still asserted, but they are cheap
	// structural checks, not the point of the case.
	std::mt19937 rng(20260818);
	std::uniform_int_distribution<int64_t> pos_dist(1, 500);

	for (int trial = 0; trial < 100; trial++) {
		CumulativeCoverageAccumulator cc;
		std::vector<std::tuple<int32_t, int64_t, int64_t>> obs;
		for (int32_t r = 0; r < 12; r++) {
			cc.AddEmptyRank(r);
			for (int j = 0; j < 4; j++) {
				int64_t a = pos_dist(rng);
				int64_t b = pos_dist(rng);
				if (a != b) {
					cc.Add(r, std::min(a, b), std::max(a, b));
					obs.emplace_back(r, std::min(a, b), std::max(a, b));
				}
			}
		}

		auto curve = cc.Curve();
		REQUIRE(curve.size() == 12);
		REQUIRE(CoveredOf(curve) == BruteForceCurve(obs, 12));
	}
}

TEST_CASE("CumulativeCoverageAccumulator - coordinates far beyond INT32 stay exact", "[cumulative_coverage]") {
	// Large POSITIONS. Note this alone does not catch 32-bit width arithmetic, because
	// the widths here are only a few thousand -- see the next case for that.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 1'000'000'000'000LL, 1'000'000'003'000LL);
	cc.Add(1, 1'000'000'002'000LL, 1'000'000'005'000LL);

	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {3000, 5000});
}

TEST_CASE("CumulativeCoverageAccumulator - a single segment wider than INT32 stays exact", "[cumulative_coverage]") {
	// The human genome is ~3.1 Gbp, so one elementary segment legitimately exceeds
	// 2^31 and a width accumulated in 32-bit arithmetic wraps NEGATIVE. This is the
	// input class that catches it; large coordinates with small widths do not.
	constexpr int64_t kWideA = 3'000'000'000LL; // > INT32_MAX
	constexpr int64_t kWideB = 2'500'000'000LL;

	CumulativeCoverageAccumulator cc;
	cc.Add(0, 1, 1 + kWideA);

	auto curve = cc.Curve();
	REQUIRE(curve.size() == 1);
	REQUIRE(curve[0].covered == kWideA);

	// And across two ranks -- abutting, so they share no base -- so the prefix sum is
	// exercised at that magnitude too. Widths are written as `start + width` rather
	// than as literal stops, because computing the stops by hand is how the first
	// version of this case ended up asserting a total that was off by one.
	CumulativeCoverageAccumulator two;
	two.Add(0, 1, 1 + kWideA);
	two.Add(1, 1 + kWideA, 1 + kWideA + kWideB);

	REQUIRE(CoveredOf(two.Curve()) == std::vector<int64_t> {kWideA, kWideA + kWideB});
}

// ---- Compact() ----
// Compaction exists to bound aggregate state: without it the accumulator holds every
// (rank, start, stop) triple until finalize, which for a deep metagenomic group is
// hundreds of MB per partial state, invisible to DuckDB's memory accounting so it OOMs
// rather than spills. compress_intervals compacts at a threshold and on every combine;
// this now does the same.
//
// It is only safe because merging a rank's intervals WITH EACH OTHER cannot change any
// base's minimum covering rank -- the base is still covered by that same rank. Merging
// across ranks would destroy the curve, and the tests below pin that it does not happen.

TEST_CASE("CumulativeCoverageAccumulator - Compact does not change the curve", "[cumulative_coverage]") {
	// The invariance that makes compaction legitimate. Randomized, with dense overlap so
	// there is plenty to merge. Fails if Compact ever merges across ranks.
	std::mt19937 rng(20260819);
	std::uniform_int_distribution<int32_t> n_ranks_dist(1, 8);
	std::uniform_int_distribution<int> n_obs_dist(0, 6);
	std::uniform_int_distribution<int64_t> pos_dist(1, 60);

	for (int trial = 0; trial < 300; trial++) {
		int32_t n_ranks = n_ranks_dist(rng);
		CumulativeCoverageAccumulator raw;
		CumulativeCoverageAccumulator compacted;
		for (int32_t r = 0; r < n_ranks; r++) {
			raw.AddEmptyRank(r);
			compacted.AddEmptyRank(r);
			int n_obs = n_obs_dist(rng);
			for (int j = 0; j < n_obs; j++) {
				int64_t a = pos_dist(rng);
				int64_t b = pos_dist(rng);
				if (a == b) {
					continue;
				}
				raw.Add(r, std::min(a, b), std::max(a, b));
				compacted.Add(r, std::min(a, b), std::max(a, b));
			}
		}

		compacted.Compact();
		REQUIRE(CoveredOf(compacted.Curve()) == CoveredOf(raw.Curve()));
	}
}

TEST_CASE("CumulativeCoverageAccumulator - Compact actually shrinks the state", "[cumulative_coverage]") {
	// If it did not shrink anything, it would be pure cost. Ten overlapping intervals on
	// one rank collapse to one; a second rank keeps its own separate run.
	CumulativeCoverageAccumulator cc;
	for (int i = 0; i < 10; i++) {
		cc.Add(0, 10 + i, 30 + i);
	}
	cc.Add(1, 500, 600);
	REQUIRE(cc.ObservationCount() == 11);

	cc.Compact();
	REQUIRE(cc.ObservationCount() == 2);
	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {29, 129});
}

TEST_CASE("CumulativeCoverageAccumulator - Compact keeps a zero-coverage rank", "[cumulative_coverage]") {
	// The case that would silently break the curve: a rank whose every observation is
	// degenerate has nothing to merge, so a naive compaction drops it, shortening the
	// curve and shifting every rank above it. Rank 1 must survive as a 0.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.AddEmptyRank(1);
	cc.Add(2, 30, 40);

	cc.Compact();
	REQUIRE(CoveredOf(cc.Curve()) == std::vector<int64_t> {10, 10, 20});

	// Same for a rank holding only inverted / zero-width intervals rather than the
	// AddEmptyRank marker.
	CumulativeCoverageAccumulator degen;
	degen.Add(0, 10, 20);
	degen.Add(1, 55, 55);
	degen.Add(1, 90, 70);
	degen.Add(2, 30, 40);

	degen.Compact();
	REQUIRE(CoveredOf(degen.Curve()) == std::vector<int64_t> {10, 10, 20});
}

TEST_CASE("CumulativeCoverageAccumulator - Compact is idempotent", "[cumulative_coverage]") {
	// Add() compacts past a threshold and Absorb() compacts every time, so a state can be
	// compacted many times over. A second pass must be a no-op in both size and answer.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(0, 15, 25);
	cc.AddEmptyRank(1);
	cc.Add(2, 30, 40);

	cc.Compact();
	const size_t after_first = cc.ObservationCount();
	auto curve_first = CoveredOf(cc.Curve());

	cc.Compact();
	REQUIRE(cc.ObservationCount() == after_first);
	REQUIRE(CoveredOf(cc.Curve()) == curve_first);
}

TEST_CASE("CumulativeCoverageAccumulator - Absorb merges correctly and does not compact a small state",
          "[cumulative_coverage]") {
	// The parallel-GROUP-BY path: combine partial states. The merged answer must equal
	// the single-state answer.
	//
	// Absorb must NOT compact unconditionally. Compacting on every Combine costs
	// O(m log m) per combine, so a parallel GROUP BY pays T sorts to save nothing when
	// the state is small -- measured 500k rows / 50 ranks going from 0.029 s at
	// threads=1 to 0.072 s at threads=16 (2.5x SLOWER), with the control case of
	// overlapping intervals, where Compact() genuinely shrinks, going 2x FASTER. The
	// threshold is what bounds the state; Combine does not need its own policy.
	CumulativeCoverageAccumulator a;
	CumulativeCoverageAccumulator b;
	for (int i = 0; i < 20; i++) {
		a.Add(0, 10 + i, 40 + i);
		b.Add(0, 20 + i, 50 + i);
	}
	a.AddEmptyRank(1);
	b.Add(1, 900, 950);

	a.Absorb(b);
	REQUIRE(a.CompactionCount() == 0);
	REQUIRE(a.ObservationCount() == 42);
	REQUIRE(CoveredOf(a.Curve()) == std::vector<int64_t> {59, 109});

	// An explicit Compact still shrinks it, and the answer is unchanged -- so the state
	// is bounded when it needs to be, just not on a schedule that costs more than it saves.
	a.Compact();
	REQUIRE(a.ObservationCount() < 42);
	REQUIRE(CoveredOf(a.Curve()) == std::vector<int64_t> {59, 109});
}

TEST_CASE("CumulativeCoverageAccumulator - compaction is amortized, not per-row", "[cumulative_coverage]") {
	// The input shape Compact() CANNOT shrink, which is also the documented shape: every
	// interval disjoint on one rank, because positions are normally compress_intervals
	// output. With a fixed threshold the first compaction fails to get below it and every
	// subsequent Add() re-sorts the whole state -- measured 17.5 ms per row past 1e6,
	// i.e. O(n^2 log n). 1.2M such rows did not finish in 90 s.
	//
	// Asserted as a BOUND ON COMPACTION COUNT rather than a wall-clock limit, so it is
	// deterministic. Just past the threshold, a growing floor compacts once and then
	// doubles out of the way; a fixed floor compacts on every one of the last 100 rows.
	const size_t n = CumulativeCoverageAccumulator::COMPACT_THRESHOLD + 100;

	CumulativeCoverageAccumulator cc;
	for (size_t i = 0; i < n; i++) {
		const int64_t start = static_cast<int64_t>(i) * 10 + 1;
		cc.Add(0, start, start + 5);
	}

	// Fixed threshold gives 101 here. Doubling gives 1.
	REQUIRE(cc.CompactionCount() <= 3);

	// The answer must be untouched by the policy: n disjoint 5-base intervals.
	const auto curve = CoveredOf(cc.Curve());
	REQUIRE(curve.size() == 1);
	REQUIRE(curve[0] == static_cast<int64_t>(n) * 5);
}

TEST_CASE("CumulativeCoverageAccumulator - self-Absorb is a no-op, not corruption", "[cumulative_coverage]") {
	// Absorb(*this) would read the source range while insert() reallocates it: undefined
	// behaviour. DuckDB never self-combines, so this is latent rather than live, but
	// Absorb is public and the header advertises it as the Combine path.
	CumulativeCoverageAccumulator cc;
	cc.Add(0, 10, 20);
	cc.Add(1, 15, 30);
	auto before = CoveredOf(cc.Curve());

	cc.Absorb(cc);
	REQUIRE(CoveredOf(cc.Curve()) == before);
}
