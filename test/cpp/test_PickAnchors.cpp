#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include "pick_anchors.hpp"

using Catch::Matchers::ContainsSubstring;

namespace {

// Ids are only hash input, so any stable naming works. Zero-padded so the
// lexicographic order the DuckDB reader imposes matches the numeric order used
// here -- otherwise a test reading "index 9" would not be reasoning about the
// same row the SQL layer would.
std::vector<std::string> Ids(uint32_t n) {
	std::vector<std::string> ids;
	ids.reserve(n);
	for (uint32_t i = 0; i < n; ++i) {
		std::string s = std::to_string(i);
		ids.push_back(std::string(4 - std::min<size_t>(4, s.size()), '0') + s);
	}
	return ids;
}

// ── FIXTURE 1: ten points evenly spaced on a line, coordinate i at index i.
//    Greedy farthest-point on this geometry is fully determined, so the exact
//    selection ORDER is assertable rather than just its shape. ──
std::vector<double> Line(uint32_t n) {
	std::vector<double> p(n);
	for (uint32_t i = 0; i < n; ++i) {
		p[i] = static_cast<double>(i);
	}
	return p;
}

// ── FIXTURE 2: three tight 2-D clusters, deliberately at UNEQUAL distances from
//    the cloud's centroid so that no step of the greedy sequence lands on a tie:
//      A = 4 points near (0,0), B = 4 near (100,0), C = 4 near (0,140).
//    centroid = (406/12, 566/12) = (33.83, 47.17); squared distance to it is
//    ~3.4e3 for A, ~6.6e3 for B, ~9.9e3 for C, so rank 0 must come from C.
//    Farthest from C is then B (~2.99e4 vs A's ~1.99e4), and the max-min step
//    after {C,B} lands in A (~1.02e4, against ~1e1 for the leftovers of B and C).
//    So the selection order is C, B, A -- one per cluster, most peripheral first. ──
const std::vector<double> CLUSTERS = {
    0,   0,   1,   0,   0,   1,   1,   1,   // A: indices 0-3
    100, 0,   101, 0,   100, 1,   101, 1,   // B: indices 4-7
    0,   140, 1,   140, 0,   141, 1,   141, // C: indices 8-11
};
char ClusterOf(uint32_t index) {
	return index < 4 ? 'A' : (index < 8 ? 'B' : 'C');
}

// ── FIXTURE 3: the proportional-allocation oracle. 100 points in 2-D whose
//    equal-frequency 2x2 grid has DELIBERATELY UNEQUAL cell sizes 40/10/10/40.
//
//    axis 0 = i, so its two equal-frequency bins are {0..49} and {50..99}.
//    axis 1 is small for {0..39} and {50..59} (50 samples) and large otherwise,
//    so its two bins are exactly those sets. The joint cells are therefore
//      (0,0) = {0..39}   40 samples      (0,1) = {40..49}   10 samples
//      (1,0) = {50..59}  10 samples      (1,1) = {60..99}   40 samples
//
//    This is the fixture that can tell proportional allocation apart from equal
//    allocation: at k=20 proportional draws 8/2/2/8 while one-per-stratum-equally
//    would draw 5/5/5/5. Nothing else about the function distinguishes them. ──
std::vector<double> Uneven() {
	std::vector<double> p;
	p.reserve(200);
	for (uint32_t i = 0; i < 100; ++i) {
		const bool small = i < 40 || (i >= 50 && i < 60);
		p.push_back(static_cast<double>(i));
		p.push_back(small ? static_cast<double>(i) : 1000.0 + static_cast<double>(i));
	}
	return p;
}
// The cell each index of Uneven() falls in, as derived in the comment above.
int UnevenCell(uint32_t i) {
	const bool low0 = i < 50;
	const bool low1 = i < 40 || (i >= 50 && i < 60);
	return (low0 ? 0 : 2) + (low1 ? 0 : 1);
}
// Cell -> number selected, for a selection over Uneven().
std::map<int, int> CellCounts(const std::vector<uint32_t> &picked) {
	std::map<int, int> counts;
	for (const auto i : picked) {
		counts[UnevenCell(i)]++;
	}
	return counts;
}

} // namespace

TEST_CASE("SelectFarthestPoint reproduces the hand-checkable order on a line", "[pick_anchors]") {
	// rank 0: farthest from the centroid (4.5) -- indices 0 and 9 tie at 4.5, so
	//         the tie-break to the lowest index decides -> 0
	// rank 1: farthest from {0}                                              -> 9
	// rank 2: max-min over {0,9} -- 4 and 5 both sit at 4, tie-break          -> 4
	// rank 3: max-min over {0,9,4} -- 2, 6 and 7 all sit at 2, tie-break      -> 2
	// Any change to the start rule, the max-min rule, or the tie-break shows up here.
	const auto picked = miint::SelectFarthestPoint(Line(10), 10, 1, 4);
	REQUIRE(picked == std::vector<uint32_t> {0, 9, 4, 2});
}

TEST_CASE("SelectFarthestPoint spreads across separated clusters", "[pick_anchors]") {
	// THE PROPERTY THAT JUSTIFIES THE RULE for diverse-subset use: a random draw
	// can miss a cluster entirely, leaving part of the space uncovered.
	// Farthest-point cannot -- asked for one anchor per cluster it must return
	// exactly one from each, in decreasing peripherality (see FIXTURE 2).
	const auto picked = miint::SelectFarthestPoint(CLUSTERS, 12, 2, 3);
	REQUIRE(picked.size() == 3);
	REQUIRE(ClusterOf(picked[0]) == 'C');
	REQUIRE(ClusterOf(picked[1]) == 'B');
	REQUIRE(ClusterOf(picked[2]) == 'A');
}

TEST_CASE("SelectFarthestPoint is prefix-stable and deterministic", "[pick_anchors]") {
	// Prefix stability is what lets a caller pick a generous k once and then trim;
	// it would break if the implementation ever re-optimized globally instead of
	// greedily. Determinism without a seed is what lets an anchor set be recorded
	// as a property of the data.
	const auto line = Line(10);
	const auto seven = miint::SelectFarthestPoint(line, 10, 1, 7);
	const auto four = miint::SelectFarthestPoint(line, 10, 1, 4);
	REQUIRE(std::vector<uint32_t>(seven.begin(), seven.begin() + 4) == four);
	REQUIRE(miint::SelectFarthestPoint(line, 10, 1, 7) == seven);

	// k == n is legal and degenerate: it just orders every sample.
	const auto all = miint::SelectFarthestPoint(line, 10, 1, 10);
	REQUIRE(std::set<uint32_t>(all.begin(), all.end()).size() == 10);
}

TEST_CASE("SelectStratified allocates proportionally to stratum size", "[pick_anchors]") {
	// THE CLAIM THE FUNCTION MAKES. On FIXTURE 3's 40/10/10/40 grid, k=20 must
	// come out 8/2/2/8 -- proportional. Equal-per-stratum allocation would give
	// 5/5/5/5, and simple random sampling would give something near 8/2/2/8 only
	// on average, not exactly. This assertion is the difference.
	const auto points = Uneven();
	const auto ids = Ids(100);
	const auto counts = CellCounts(miint::SelectStratified(points, ids, 100, 2, 20, 2, 42));
	REQUIRE(counts == std::map<int, int> {{0, 8}, {1, 2}, {2, 2}, {3, 8}});
}

TEST_CASE("SelectStratified re-draws within strata when the seed changes", "[pick_anchors]") {
	// The two halves of the rule, asserted together: the seed must change WHICH
	// samples are taken (unbiased within a stratum) without changing HOW MANY come
	// from each (coverage across strata is fixed by size alone). A rule that let
	// the seed move the counts would not be proportional allocation; one that let
	// it change nothing would not be a random draw.
	const auto points = Uneven();
	const auto ids = Ids(100);
	const auto a = miint::SelectStratified(points, ids, 100, 2, 20, 2, 1);
	const auto b = miint::SelectStratified(points, ids, 100, 2, 20, 2, 2);
	REQUIRE(CellCounts(a) == CellCounts(b));
	REQUIRE(std::set<uint32_t>(a.begin(), a.end()) != std::set<uint32_t>(b.begin(), b.end()));

	// Determinism at a fixed seed -- the counterpart of the above.
	REQUIRE(miint::SelectStratified(points, ids, 100, 2, 20, 2, 1) == a);
}

TEST_CASE("SelectStratified is unbiased within a stratum", "[pick_anchors]") {
	// THE MECHANISM that makes this rule beat random rather than joining the
	// losers. Leverage, medoid and farthest-point selection all lose because each
	// systematically prefers a KIND of point; stratified wins only because it has
	// no preference inside a stratum. So no sample may be structurally
	// unselectable: over 100 seeds every one of the 100 samples must be drawn at
	// least once. A hash keyed on position instead of identity, or any geometric
	// tie-break leaking into the within-stratum order, would leave the
	// never-selected set non-empty here.
	const auto points = Uneven();
	const auto ids = Ids(100);
	std::set<uint32_t> ever;
	for (int64_t seed = 0; seed < 100; ++seed) {
		const auto picked = miint::SelectStratified(points, ids, 100, 2, 20, 2, seed);
		ever.insert(picked.begin(), picked.end());
	}
	REQUIRE(ever.size() == 100);
}

TEST_CASE("SelectStratified is prefix-stable", "[pick_anchors]") {
	// Same contract as the farthest-point rule, and it holds for the same reason:
	// the global ordering is built without reference to k, so a k-selection's first
	// m entries are the m-selection. On FIXTURE 3, k=10 is the 4/1/1/4 allocation
	// nested inside k=20's 8/2/2/8.
	const auto points = Uneven();
	const auto ids = Ids(100);
	const auto twenty = miint::SelectStratified(points, ids, 100, 2, 20, 2, 42);
	const auto ten = miint::SelectStratified(points, ids, 100, 2, 10, 2, 42);
	REQUIRE(CellCounts(ten) == std::map<int, int> {{0, 4}, {1, 1}, {2, 1}, {3, 4}});
	REQUIRE(std::vector<uint32_t>(twenty.begin(), twenty.begin() + 10) == ten);
}

TEST_CASE("SelectStratified depends on sample identity, not row position", "[pick_anchors]") {
	// The hash is keyed on the id string so that a selection is a property of the
	// samples rather than of the order they happened to arrive in. Renaming the
	// samples must therefore change the draw even though the geometry is identical.
	const auto points = Uneven();
	const auto ids = Ids(100);
	std::vector<std::string> renamed;
	renamed.reserve(ids.size());
	for (const auto &id : ids) {
		renamed.push_back("s" + id);
	}
	const auto a = miint::SelectStratified(points, ids, 100, 2, 20, 2, 42);
	const auto b = miint::SelectStratified(points, renamed, 100, 2, 20, 2, 42);
	REQUIRE(CellCounts(a) == CellCounts(b));
	REQUIRE(std::set<uint32_t>(a.begin(), a.end()) != std::set<uint32_t>(b.begin(), b.end()));
}

TEST_CASE("SelectStratified degenerates gracefully when k <= 1 stratum", "[pick_anchors]") {
	// n_bins = 1 collapses the grid to a single stratum, which is exactly simple
	// random sampling -- the baseline the rule is measured against. It must stay
	// well-defined rather than dividing by an empty stratum count.
	const auto points = Uneven();
	const auto ids = Ids(100);
	const auto picked = miint::SelectStratified(points, ids, 100, 2, 20, 1, 42);
	REQUIRE(picked.size() == 20);
	REQUIRE(std::set<uint32_t>(picked.begin(), picked.end()).size() == 20);

	// k == n returns every sample regardless of how the strata fell out.
	const auto all = miint::SelectStratified(points, ids, 100, 2, 100, 2, 42);
	REQUIRE(std::set<uint32_t>(all.begin(), all.end()).size() == 100);
}

TEST_CASE("SelectStratified handles more bins than samples", "[pick_anchors]") {
	// n < n_bins is the other edge of the equal-frequency binning: every sample
	// gets its own bin, so the grid must not emit a bin index outside [0, n_bins)
	// and corrupt the stratum encoding.
	const auto picked = miint::SelectStratified(Line(10), Ids(10), 10, 1, 3, 16, 42);
	REQUIRE(picked.size() == 3);
	REQUIRE(std::set<uint32_t>(picked.begin(), picked.end()).size() == 3);
}

TEST_CASE("pick_anchors selectors reject impossible arguments", "[pick_anchors]") {
	const auto points = Uneven();
	const auto ids = Ids(100);
	// k must be a real request, and cannot exceed what exists.
	REQUIRE_THROWS_AS(miint::SelectStratified(points, ids, 100, 2, 0, 2, 42), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectStratified(points, ids, 100, 2, 101, 2, 42), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectFarthestPoint(points, 100, 2, 0), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectFarthestPoint(points, 100, 2, 101), std::invalid_argument);
	// A zero-dimensional or empty cloud has no geometry to select on.
	REQUIRE_THROWS_AS(miint::SelectStratified(points, ids, 100, 0, 5, 2, 42), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectFarthestPoint(points, 0, 2, 1), std::invalid_argument);
	// n_bins < 1 would make equal-frequency binning meaningless.
	REQUIRE_THROWS_AS(miint::SelectStratified(points, ids, 100, 2, 5, 0, 42), std::invalid_argument);
	// Shape mismatches are caller bugs and must not be read past the end.
	REQUIRE_THROWS_AS(miint::SelectStratified(points, ids, 100, 3, 5, 2, 42), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectStratified(points, Ids(99), 100, 2, 5, 2, 42), std::invalid_argument);
	REQUIRE_THROWS_AS(miint::SelectFarthestPoint(points, 100, 3, 5), std::invalid_argument);
	// A grid too large to index must be refused rather than silently aliased into
	// the wrong strata. 3000000^3 overflows a 64-bit cell index; two axes cannot
	// overflow it at all, since UINT32_MAX^2 still fits.
	const std::vector<double> tiny3d = {0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3};
	REQUIRE_THROWS_WITH(miint::SelectStratified(tiny3d, Ids(4), 4, 3, 2, 3000000u, 42), ContainsSubstring("strata"));
}
