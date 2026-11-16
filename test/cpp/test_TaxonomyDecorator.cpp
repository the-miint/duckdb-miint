#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "taxonomy_decorator_fixtures.hpp"

#include <compare>
#include <cstdint>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace miint;
using namespace miint::test;
using Catch::Matchers::ContainsSubstring;
using Catch::Matchers::WithinRel;

namespace {
enum Rank : size_t { kDomain, kPhylum, kClass, kOrder };

// A hand-built tree, numbered out of preorder so that a node id read as a preorder position,
// or the reverse, gives wrong counts:
//
//   root(3) +-- n1(7) +-- t4(1)  A: D1; P1; C1
//           |         +-- t5(6)  B: D1; P1; C2; O1
//           +-- n2(0) --- n3(8) +-- t6(2)  C: D1; P2; C3
//               (one child)     +-- t7(4)  D: D1 only, so cleaning wipes it (no phylum)
//                               +-- t8(5)  no row
//   Row E (D1; P3) is on no tip.
constexpr uint32_t kRoot = 3, kN1 = 7, kN2 = 0, kN3 = 8, kT4 = 1, kT5 = 6, kT6 = 2, kT7 = 4, kT8 = 5;
constexpr uint32_t kNodes = 9;
constexpr uint32_t NP = TreeIndex::NO_PARENT;

TaxonomyRow Row(std::string id, std::initializer_list<std::string> names) {
	TaxonomyRow row {std::move(id), {}};
	size_t i = 0;
	for (const auto &name : names) {
		row.ranks[i++] = name;
	}
	return row;
}

TaxonomyStore HandStore() {
	return TaxonomyStore({Row("A", {"D1", "P1", "C1"}), Row("B", {"D1", "P1", "C2", "O1"}),
	                      Row("C", {"D1", "P2", "C3"}), Row("D", {"D1"}), Row("E", {"D1", "P3"})});
}

TreeIndex HandTree() {
	std::vector<uint32_t> parent(kNodes);
	parent[kRoot] = NP;
	parent[kN1] = kRoot;
	parent[kN2] = kRoot;
	parent[kN3] = kN2;
	parent[kT4] = kN1;
	parent[kT5] = kN1;
	parent[kT6] = kN3;
	parent[kT7] = kN3;
	parent[kT8] = kN3;
	return TreeIndex::Build(parent, std::vector<double>(kNodes, 1.0));
}

std::vector<std::optional<size_t>> HandRows() {
	std::vector<std::optional<size_t>> rows(kNodes);
	rows[kT4] = 0;
	rows[kT5] = 1;
	rows[kT6] = 2;
	rows[kT7] = 3;
	return rows;
}

// Per-node values, from (node, value) pairs; unlisted nodes are 0.
std::vector<uint32_t> ByNode(std::initializer_list<std::pair<uint32_t, uint32_t>> values) {
	std::vector<uint32_t> out(kNodes, 0);
	for (const auto &[node, value] : values) {
		out[node] = value;
	}
	return out;
}

template <class F>
std::vector<uint32_t> Observed(F f) {
	std::vector<uint32_t> out(kNodes);
	for (uint32_t v = 0; v < kNodes; v++) {
		out[v] = f(v);
	}
	return out;
}

std::string ErrorOf(const std::vector<std::optional<size_t>> &rows) {
	auto store = HandStore();
	try {
		CladeCountIndex index(HandTree(), store, rows);
	} catch (const std::invalid_argument &e) {
		return e.what();
	}
	return "";
}

int Sign(std::strong_ordering o) {
	return o == std::strong_ordering::less ? -1 : o == std::strong_ordering::greater ? 1 : 0;
}
} // namespace

// Precision is c / NumTips and the fold breaks ties by LeafCount, so the two must differ in
// exactly one way: a tip with no row, or one whose lineage cleaning wiped, is a leaf but not an
// informative tip. Counts must also follow the tree through a one-child node and a
// three-child node.
TEST_CASE("CladeCountIndex counts a hand-built tree", "[taxonomy_decorator]") {
	auto store = HandStore();
	REQUIRE_FALSE(store.IsInformative(3));
	const NameId d1 = store.Name(0, kDomain), p1 = store.Name(0, kPhylum), p2 = store.Name(2, kPhylum);
	const NameId c1 = store.Name(0, kClass), c2 = store.Name(1, kClass), c3 = store.Name(2, kClass);
	const NameId o1 = store.Name(1, kOrder);
	REQUIRE(store.NameText(kPhylum, p2) == "P2");
	REQUIRE(store.NameText(kClass, c3) == "C3");

	CladeCountIndex index(HandTree(), store, HandRows());
	REQUIRE(index.Tree().size() == kNodes);

	CHECK(Observed([&](uint32_t v) { return index.LeafCount(v); }) ==
	      ByNode({{kRoot, 5}, {kN1, 2}, {kN2, 3}, {kN3, 3}, {kT4, 1}, {kT5, 1}, {kT6, 1}, {kT7, 1}, {kT8, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.NumTips(v); }) ==
	      ByNode({{kRoot, 3}, {kN1, 2}, {kN2, 1}, {kN3, 1}, {kT4, 1}, {kT5, 1}, {kT6, 1}}));

	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kDomain, d1); }) ==
	      ByNode({{kRoot, 3}, {kN1, 2}, {kN2, 1}, {kN3, 1}, {kT4, 1}, {kT5, 1}, {kT6, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kPhylum, p1); }) ==
	      ByNode({{kRoot, 2}, {kN1, 2}, {kT4, 1}, {kT5, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kPhylum, p2); }) ==
	      ByNode({{kRoot, 1}, {kN2, 1}, {kN3, 1}, {kT6, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kClass, c1); }) ==
	      ByNode({{kRoot, 1}, {kN1, 1}, {kT4, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kClass, c2); }) ==
	      ByNode({{kRoot, 1}, {kN1, 1}, {kT5, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kClass, c3); }) ==
	      ByNode({{kRoot, 1}, {kN2, 1}, {kN3, 1}, {kT6, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kOrder, o1); }) ==
	      ByNode({{kRoot, 1}, {kN1, 1}, {kT5, 1}}));

	CHECK(index.Total(kPhylum, p1) == 2);
	CHECK(index.Total(kPhylum, p2) == 1);
	CHECK(index.Total(kClass, c1) == 1);
	CHECK(index.Total(kOrder, o1) == 1);
}

// Total is the denominator of recall, so it counts the tree's informative tips, not the
// store's rows: D1 is on five rows, but one sits on no tip and one was wiped, so its total is
// 3. A name only on rows that no tip carries (P3) counts nowhere.
TEST_CASE("CladeCountIndex totals count the tree's tips, not the taxonomy's rows", "[taxonomy_decorator]") {
	auto store = HandStore();
	const NameId d1 = store.Name(0, kDomain), p3 = store.Name(4, kPhylum);
	REQUIRE(store.NameText(kPhylum, p3) == "P3");
	CladeCountIndex index(HandTree(), store, HandRows());

	CHECK(index.Total(kDomain, d1) == 3);
	CHECK(index.Total(kPhylum, p3) == 0);
	CHECK(Observed([&](uint32_t v) { return index.CladeCount(v, kPhylum, p3); }) == ByNode({}));
}

// A row is a tip's taxonomy: one on an internal node, or one that does not exist, is a caller
// bug that would otherwise be counted silently or read out of bounds.
TEST_CASE("CladeCountIndex rejects a bad row mapping", "[taxonomy_decorator]") {
	auto rows = HandRows();
	rows.pop_back();
	auto message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("8 entries"));
	CHECK_THAT(message, ContainsSubstring("9 nodes"));

	rows = HandRows();
	rows[kN1] = 0;
	CHECK_THAT(ErrorOf(rows), ContainsSubstring("node 7 is not a tip"));

	rows = HandRows();
	rows[kT8] = 99;
	message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("row 99"));
	CHECK_THAT(message, ContainsSubstring("5 rows"));
}

// A node is a candidate when its count reaches max(min_count, ceil(total/2)): recall >= 1/2,
// tested in integers, so a clade holding exactly half of an even total qualifies. min_count
// only matters for small totals, and 0 would admit nodes with no informative tips.
TEST_CASE("F1Score::Threshold is max(min_count, ceil(total/2))", "[taxonomy_decorator]") {
	struct Case {
		uint32_t total, min_count, expected;
	};
	const Case cases[] = {{0, 1, 1}, {0, 2, 2}, {1, 1, 1}, {1, 2, 2}, {2, 1, 1}, {2, 2, 2},  {3, 1, 2},
	                      {3, 2, 2}, {4, 1, 2}, {4, 2, 2}, {5, 1, 3}, {5, 2, 3}, {10, 2, 5}, {11, 2, 6}};
	for (const auto &c : cases) {
		INFO("total " << c.total << ", min_count " << c.min_count);
		CHECK(F1Score::Threshold(c.total, c.min_count) == c.expected);
	}
	// The largest total a tree can have; ceil must not overflow.
	CHECK(F1Score::Threshold(F1Score::kMaxTips, 2) == F1Score::kMaxTips / 2);
	CHECK_THROWS_AS(F1Score::Threshold(4, 0), std::invalid_argument);
}

// The fold keeps the best-scoring candidate for a name, so Compare must order by
// F1 = 2c/(total+m) exactly: equal rationals tie (the fold's tie-break then decides), and
// values doubles cannot tell apart are still ordered.
TEST_CASE("F1Score::Compare orders candidates by F1, exactly", "[taxonomy_decorator]") {
	// Same c, fewer tips below: higher precision wins.
	CHECK(Sign(F1Score::Compare(3, 4, 3, 6, 4)) == 1);
	CHECK(Sign(F1Score::Compare(3, 6, 3, 4, 4)) == -1);
	// More of the name but much less precise: 8/12 against 6/7.
	CHECK(Sign(F1Score::Compare(4, 8, 3, 3, 4)) == -1);
	// Different (c, m), equal F1: 4/6 and 6/9.
	CHECK(Sign(F1Score::Compare(2, 2, 3, 5, 4)) == 0);
	CHECK(Sign(F1Score::Compare(3, 5, 2, 2, 4)) == 0);

	// At the tip cap the cross products reach 2^63 and must not overflow.
	constexpr uint32_t kCap = F1Score::kMaxTips;
	CHECK(Sign(F1Score::Compare(kCap, kCap, kCap, kCap, kCap)) == 0);
	CHECK(Sign(F1Score::Compare(kCap, kCap, kCap - 1, kCap, kCap)) == 1);
	// Cross products about 7.8e18 that differ by exactly 1, so the two F1 values are the same
	// double: only exact integer arithmetic orders them.
	const uint32_t c1 = 1839998538, m1 = 2139392963, c2 = 1818229075, m2 = 2088673900;
	REQUIRE(F1Score::ScoreForReport(c1, m1, kCap) == F1Score::ScoreForReport(c2, m2, kCap));
	CHECK(Sign(F1Score::Compare(c1, m1, c2, m2, kCap)) == -1);
	CHECK(Sign(F1Score::Compare(c2, m2, c1, m1, kCap)) == 1);
}

// ScoreForReport is what a user reads beside a placement, so it must be F1 itself.
TEST_CASE("F1Score::ScoreForReport is 2c/(total+m)", "[taxonomy_decorator]") {
	CHECK_THAT(F1Score::ScoreForReport(3, 4, 4), WithinRel(0.75, 1e-15));
	CHECK_THAT(F1Score::ScoreForReport(2, 2, 4), WithinRel(2.0 / 3.0, 1e-15));
	CHECK_THAT(F1Score::ScoreForReport(1, 1, 1), WithinRel(1.0, 1e-15));
}

// Every later stage reads these counts, so they must equal a count made the slow way, on trees
// with chains, multifurcations, unmapped and wiped tips, and names down to strain.
TEST_CASE("CladeCountIndex matches counts made by walking up from every tip", "[taxonomy_decorator]") {
	std::mt19937 rng(20261116);
	for (int trial = 0; trial < 300; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		const auto index = build_index(c, store);
		const uint32_t root = root_of(c.tree);
		INFO("trial " << trial << ", " << n << " nodes");
		for (uint32_t v = 0; v < n; v++) {
			INFO("node " << v);
			REQUIRE(index.LeafCount(v) == naive.leaf_count[v]);
			REQUIRE(index.NumTips(v) == naive.num_tips[v]);
			if (v != root) {
				REQUIRE(index.NumTips(c.tree.parent[v]) >= index.NumTips(v));
			}
			for (size_t r = 0; r < kRankCount; r++) {
				for (NameId name = 0; name < store.NameCount(r); name++) {
					REQUIRE(index.CladeCount(v, r, name) == naive.clade[r][name][v]);
				}
			}
		}
		for (size_t r = 0; r < kRankCount; r++) {
			for (NameId name = 0; name < store.NameCount(r); name++) {
				REQUIRE(index.Total(r, name) == naive.Total(r, name, root));
			}
		}
	}
}

// Placement searches for a name's lowest candidates along at most two root paths. That is
// sound only if candidates are closed upward and there are at most two lowest ones, and two
// only when two disjoint clades each hold exactly half of an even total.
TEST_CASE("CladeCountIndex candidates are closed upward with at most two lowest", "[taxonomy_decorator]") {
	std::mt19937 rng(20261117);
	bool saw_two = false;
	for (int trial = 0; trial < 300; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		INFO("trial " << trial << ", " << n << " nodes");
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			for (size_t r = 0; r < kRankCount; r++) {
				for (NameId name = 0; name < store.NameCount(r); name++) {
					const uint32_t total = index.Total(r, name);
					const uint32_t threshold = F1Score::Threshold(total, min_count);
					auto eligible = [&](uint32_t v) {
						return index.CladeCount(v, r, name) >= threshold;
					};
					// A lowest candidate: eligible, with no eligible child.
					std::vector<bool> has_eligible_child(n, false);
					for (uint32_t v = 0; v < n; v++) {
						if (eligible(v) && c.tree.parent[v] != kNoParent) {
							REQUIRE(eligible(c.tree.parent[v]));
							has_eligible_child[c.tree.parent[v]] = true;
						}
					}
					std::vector<uint32_t> lowest;
					for (uint32_t v = 0; v < n; v++) {
						if (eligible(v) && !has_eligible_child[v]) {
							lowest.push_back(v);
						}
					}
					INFO("min_count " << min_count << ", rank " << r << ", name " << name << ", total " << total);
					REQUIRE(lowest.size() <= 2);
					if (lowest.size() == 2) {
						saw_two = true;
						REQUIRE(total % 2 == 0);
						CHECK(index.CladeCount(lowest[0], r, name) == total / 2);
						CHECK(index.CladeCount(lowest[1], r, name) == total / 2);
					}
				}
			}
		}
	}
	// The 50/50 case must actually occur, or this test proves nothing about it.
	CHECK(saw_two);
}

// Compare's integer formula, against F1 reduced to lowest terms: equal reduced fractions must
// tie, and every other pair must be ordered as the fractions are.
TEST_CASE("F1Score::Compare agrees with exact fractions", "[taxonomy_decorator]") {
	std::mt19937 rng(20261118);
	bool saw_tie = false, saw_less = false, saw_greater = false;
	auto random_candidate = [&](uint32_t total) {
		const uint32_t c = std::uniform_int_distribution<uint32_t>(1, total)(rng);
		const uint32_t m = std::uniform_int_distribution<uint32_t>(c, 40)(rng);
		return std::pair<uint32_t, uint32_t>(c, m);
	};
	for (int trial = 0; trial < 20000; trial++) {
		const uint32_t total = std::uniform_int_distribution<uint32_t>(1, 30)(rng);
		const auto [c1, m1] = random_candidate(total);
		const auto [c2, m2] = random_candidate(total);
		const uint64_t n1 = 2 * c1, d1 = total + m1, g1 = std::gcd(n1, d1);
		const uint64_t n2 = 2 * c2, d2 = total + m2, g2 = std::gcd(n2, d2);
		int expected;
		if (n1 / g1 == n2 / g2 && d1 / g1 == d2 / g2) {
			expected = 0;
		} else {
			const long double f1 = static_cast<long double>(n1) / d1, f2 = static_cast<long double>(n2) / d2;
			expected = f1 < f2 ? -1 : 1;
		}
		saw_tie |= expected == 0 && (c1 != c2 || m1 != m2);
		saw_less |= expected < 0;
		saw_greater |= expected > 0;
		INFO("(" << c1 << ", " << m1 << ") vs (" << c2 << ", " << m2 << "), total " << total);
		REQUIRE(Sign(F1Score::Compare(c1, m1, c2, m2, total)) == expected);
	}
	CHECK(saw_tie);
	CHECK(saw_less);
	CHECK(saw_greater);
}
