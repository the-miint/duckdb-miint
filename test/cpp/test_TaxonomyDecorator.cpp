#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "taxonomy_decorator_fixtures.hpp"

#include <algorithm>
#include <array>
#include <compare>
#include <cstdint>
#include <functional>
#include <numeric>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <string_view>
#include <tuple>
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
	parent[kRoot] = kNoParent;
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

// Precision is c / InformativeTipCount and score ties are broken by TipCount, so the two must differ in
// exactly one way: a tip with no row, or one whose lineage cleaning wiped, is a tip but not an
// informative one. Counts must also follow the tree through a one-child node and a
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

	CHECK(Observed([&](uint32_t v) { return index.TipCount(v); }) ==
	      ByNode({{kRoot, 5}, {kN1, 2}, {kN2, 3}, {kN3, 3}, {kT4, 1}, {kT5, 1}, {kT6, 1}, {kT7, 1}, {kT8, 1}}));
	CHECK(Observed([&](uint32_t v) { return index.InformativeTipCount(v); }) ==
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

	// The boundaries: one entry too many, and the first row past the end.
	rows = HandRows();
	rows.push_back(std::nullopt);
	message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("10 entries"));
	CHECK_THAT(message, ContainsSubstring("9 nodes"));

	rows = HandRows();
	rows[kT8] = 5;
	message = ErrorOf(rows);
	CHECK_THAT(message, ContainsSubstring("row 5,"));
	CHECK_THAT(message, ContainsSubstring("5 rows"));
}

// Pinned here (the rule itself is documented at F1Score::Threshold): a clade holding exactly
// half of an even total is eligible, min_count matters only for small totals, and 0 is
// rejected.
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

// A name is assigned to its best-scoring node, so Compare must order exactly: equal rationals
// tie, and values doubles cannot tell apart are still ordered.
TEST_CASE("F1Score::Compare orders scores by F1, exactly", "[taxonomy_decorator]") {
	// Same c, fewer tips below: higher precision wins.
	CHECK(Sign(F1Score::Compare(3, 4, 3, 6, 4)) == 1);
	CHECK(Sign(F1Score::Compare(3, 6, 3, 4, 4)) == -1);
	// More of the name but much less precise: 8/12 against 6/7.
	CHECK(Sign(F1Score::Compare(4, 8, 3, 3, 4)) == -1);
	// Different (c, m), equal F1: 4/6 and 6/9. A divergence from tax2tree, whose doubles rank
	// the first higher (see F1Score).
	CHECK(Sign(F1Score::Compare(2, 2, 3, 5, 4)) == 0);
	CHECK(Sign(F1Score::Compare(3, 5, 2, 2, 4)) == 0);

	// At the tip cap the cross products reach 2^63 and must not overflow.
	constexpr uint32_t kCap = F1Score::kMaxTips;
	CHECK(Sign(F1Score::Compare(kCap, kCap, kCap, kCap, kCap)) == 0);
	CHECK(Sign(F1Score::Compare(kCap, kCap, kCap - 1, kCap, kCap)) == 1);
	// Cross products c1(total+m2) = 7794523695057664824 and c2(total+m1) = 7794523695057664825:
	// (c1, m1)'s cross product is lower by one, so its F1 is lower, yet the two F1 values are the
	// same double: only exact integer arithmetic orders them.
	const uint32_t c1 = 1839998538, m1 = 2139392963, c2 = 1818229075, m2 = 2088673900;
	REQUIRE(F1Score::ScoreForReport(c1, m1, kCap) == F1Score::ScoreForReport(c2, m2, kCap));
	CHECK(Sign(F1Score::Compare(c1, m1, c2, m2, kCap)) == -1);
	CHECK(Sign(F1Score::Compare(c2, m2, c1, m1, kCap)) == 1);
}

// ScoreForReport is what a user reads beside an assigned name, so it must be F1 itself.
TEST_CASE("F1Score::ScoreForReport is 2c/(total+m)", "[taxonomy_decorator]") {
	CHECK_THAT(F1Score::ScoreForReport(3, 4, 4), WithinRel(0.75, 1e-15));
	CHECK_THAT(F1Score::ScoreForReport(2, 2, 4), WithinRel(2.0 / 3.0, 1e-15));
	CHECK_THAT(F1Score::ScoreForReport(1, 1, 1), WithinRel(1.0, 1e-15));
}

// Decoration reads every count from this index, so each must equal the count found by walking
// up from every tip, on trees with chains, multifurcations, tips with no row or a wiped one, and
// names down to strain.
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
			REQUIRE(index.TipCount(v) == naive.tip_count[v]);
			REQUIRE(index.InformativeTipCount(v) == naive.informative_tip_count[v]);
			if (v != root) {
				REQUIRE(index.InformativeTipCount(c.tree.parent[v]) >= index.InformativeTipCount(v));
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

// Decoration searches at most two root paths for a name's lowest eligible nodes, which part 1
// of F1Score's contract makes sound. This checks that contract on the index's counts.
TEST_CASE("CladeCountIndex eligible nodes are closed upward with at most two lowest", "[taxonomy_decorator]") {
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
					// A lowest eligible node: one with no eligible child.
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
	auto random_score = [&](uint32_t total) {
		const uint32_t c = std::uniform_int_distribution<uint32_t>(1, total)(rng);
		const uint32_t m = std::uniform_int_distribution<uint32_t>(c, 40)(rng);
		return std::pair<uint32_t, uint32_t>(c, m);
	};
	for (int trial = 0; trial < 20000; trial++) {
		const uint32_t total = std::uniform_int_distribution<uint32_t>(1, 30)(rng);
		const auto [c1, m1] = random_score(total);
		const auto [c2, m2] = random_score(total);
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

// LowestEligibleNodes walks up from a name's middle tip in preorder, so it needs each name's
// tips by preorder position and the node at a position.
TEST_CASE("CladeCountIndex lists a hand-built tree's names in preorder", "[taxonomy_decorator]") {
	auto store = HandStore();
	const NameId d1 = store.Name(0, kDomain), p1 = store.Name(0, kPhylum), p2 = store.Name(2, kPhylum);
	const NameId p3 = store.Name(4, kPhylum);
	CladeCountIndex index(HandTree(), store, HandRows());

	// Children are visited in ascending index, so n2(0) comes before n1(7).
	std::vector<uint32_t> node_at(kNodes);
	for (uint32_t p = 0; p < kNodes; p++) {
		node_at[p] = index.NodeAtPreorder(p);
	}
	CHECK(node_at == std::vector<uint32_t> {kRoot, kN2, kN3, kT6, kT7, kT8, kN1, kT4, kT5});

	auto occurrences = [&](size_t rank, NameId name) {
		const auto occ = index.Occurrences(rank, name);
		return std::vector<uint32_t>(occ.begin(), occ.end());
	};
	// t6, t4 and t5 are at positions 3, 7 and 8. The wiped t7 and the rowless t8 carry nothing.
	CHECK(occurrences(kDomain, d1) == std::vector<uint32_t> {3, 7, 8});
	CHECK(occurrences(kPhylum, p1) == std::vector<uint32_t> {7, 8});
	CHECK(occurrences(kPhylum, p2) == std::vector<uint32_t> {3});
	CHECK(occurrences(kPhylum, p3).empty());

	// Every name in the store has an id, including P3, which is on no tip.
	CHECK(index.NameCount(kDomain) == 1);
	CHECK(index.NameCount(kPhylum) == 3);
	CHECK(index.NameCount(kClass) == 3);
	CHECK(index.NameCount(kOrder) == 1);
	CHECK(index.NameCount(kRankCount - 1) == 0);
}

// The same on random trees: each name's list must be exactly the informative tips carrying it,
// strictly ascending.
TEST_CASE("CladeCountIndex lists each name's tips in preorder", "[taxonomy_decorator]") {
	std::mt19937 rng(20261119);
	for (int trial = 0; trial < 300; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		const auto index = build_index(c, store);
		const uint32_t root = root_of(c.tree);
		INFO("trial " << trial << ", " << n << " nodes");
		for (uint32_t v = 0; v < n; v++) {
			REQUIRE(index.NodeAtPreorder(index.Tree().pre(v)) == v);
		}
		for (size_t r = 0; r < kRankCount; r++) {
			REQUIRE(index.NameCount(r) == store.NameCount(r));
			for (NameId name = 0; name < store.NameCount(r); name++) {
				INFO("rank " << r << ", name " << name);
				const auto occ = index.Occurrences(r, name);
				REQUIRE(occ.size() == naive.Total(r, name, root));
				REQUIRE(std::adjacent_find(occ.begin(), occ.end(), std::greater_equal<>()) == occ.end());
				for (uint32_t p : occ) {
					const uint32_t v = index.NodeAtPreorder(p);
					REQUIRE(index.Tree().is_tip(v));
					REQUIRE(naive.clade[r][name][v] == 1);
				}
			}
		}
	}
}

namespace {
// An index over a tree given by parent pointers, where each listed tip carries its lineage
// (domain first) and every other tip carries no row.
struct HandCase {
	TaxonomyStore store;
	CladeCountIndex index;

	NameId Name(size_t rank, std::string_view text) const {
		for (NameId k = 0; k < store.NameCount(rank); k++) {
			if (store.NameText(rank, k) == text) {
				return k;
			}
		}
		throw std::logic_error("HandCase: no name " + std::string(text));
	}
};

using Lineages = std::vector<std::pair<uint32_t, std::vector<std::string>>>;

HandCase MakeCase(const std::vector<uint32_t> &parent, const Lineages &lineage_of) {
	std::vector<TaxonomyRow> rows;
	std::vector<std::optional<size_t>> row_of_node(parent.size());
	for (const auto &[tip, lineage] : lineage_of) {
		TaxonomyRow row {"T" + std::to_string(tip), {}};
		std::copy(lineage.begin(), lineage.end(), row.ranks.begin());
		row_of_node[tip] = rows.size();
		rows.push_back(std::move(row));
	}
	TaxonomyStore store(rows);
	CladeCountIndex index(TreeIndex::Build(parent, std::vector<double>(parent.size(), 1.0)), store, row_of_node);
	return {std::move(store), std::move(index)};
}

// Each of `tips` carries the lineage (D1; N), so N is a phylum.
Lineages WithN(std::initializer_list<uint32_t> tips) {
	Lineages out;
	for (uint32_t tip : tips) {
		out.push_back({tip, {"D1", "N"}});
	}
	return out;
}

std::vector<uint32_t> LowestOf(const HandCase &c, std::string_view phylum, uint32_t min_count) {
	const NameId name = c.Name(kPhylum, phylum);
	return LowestEligibleNodes(c.index, kPhylum, name, F1Score::Threshold(c.index.Total(kPhylum, name), min_count));
}
} // namespace

// Only a middle tip of a name, in preorder, is sure to lie below a lowest eligible node; every
// eligible ancestor of another tip can be above one. Here a walk up from the name's first tip,
// a, stops at the root, yet B holds 9 of the 10 and is eligible below it. A name found only at
// the root would be assigned to the wrong node. The mirror puts a last, against a walk up from
// the last tip.
TEST_CASE("LowestEligibleNodes finds a clade below the first tip's eligible ancestor", "[taxonomy_decorator]") {
	// root(0) +-- a(1)
	//         +-- B(2) --- tips 3..11
	std::vector<uint32_t> parent(12, 2);
	parent[0] = kNoParent;
	parent[1] = 0;
	parent[2] = 0;
	auto first = MakeCase(parent, WithN({1, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
	REQUIRE(first.index.Total(kPhylum, first.Name(kPhylum, "N")) == 10);
	CHECK(LowestOf(first, "N", 2) == std::vector<uint32_t> {2});

	// root(0) +-- B(1) --- tips 2..10
	//         +-- a(11)
	std::vector<uint32_t> mirror(12, 1);
	mirror[0] = kNoParent;
	mirror[1] = 0;
	mirror[11] = 0;
	auto last = MakeCase(mirror, WithN({2, 3, 4, 5, 6, 7, 8, 9, 10, 11}));
	CHECK(LowestOf(last, "N", 2) == std::vector<uint32_t> {1});
}

// A name split evenly between two disjoint clades has both as lowest eligible nodes: neither
// holds more of it. They come in preorder, which here is not index order, and the one-child p
// above x is eligible too but holds x.
TEST_CASE("LowestEligibleNodes returns both halves of a 50/50 split, in preorder", "[taxonomy_decorator]") {
	// root(0) +-- p(1) --- x(5) +-- t6(6)  N
	//         |                 +-- t7(7)  N
	//         +-- y(2) +-- t3(3)  N
	//                  +-- t4(4)  N
	const std::vector<uint32_t> parent {kNoParent, 0, 0, 2, 2, 1, 5, 5};
	auto c = MakeCase(parent, WithN({3, 4, 6, 7}));
	CHECK(LowestOf(c, "N", 2) == std::vector<uint32_t> {5, 2});
	CHECK(LowestOf(c, "N", 1) == std::vector<uint32_t> {5, 2});
}

// With an even total the middle of the name's tips is two tips. When the walks up from them end
// at nested nodes, only the deeper is lowest: the other holds it. Both orders are checked, so
// neither walk's result is kept just for being first or second.
TEST_CASE("LowestEligibleNodes keeps the deeper of two nested results", "[taxonomy_decorator]") {
	// root(0) +-- q(1) +-- r(2) +-- t4(4)  N
	//         |        |        +-- t5(5)  N
	//         |        +-- t3(3)  N
	//         +-- t6(6)  N
	// The middle tips are t5, whose walk ends at r, and t3, whose walk ends at q.
	auto deeper_first = MakeCase({kNoParent, 0, 1, 1, 2, 2, 0}, WithN({3, 4, 5, 6}));
	CHECK(LowestOf(deeper_first, "N", 2) == std::vector<uint32_t> {2});

	// root(0) +-- t1(1)  N
	//         +-- q(2) +-- t3(3)  N
	//                  +-- r(4) +-- t5(5)  N
	//                           +-- t6(6)  N
	// The middle tips are t3, whose walk ends at q, and t5, whose walk ends at r.
	auto deeper_second = MakeCase({kNoParent, 0, 0, 2, 2, 4, 4}, WithN({1, 3, 5, 6}));
	CHECK(LowestOf(deeper_second, "N", 2) == std::vector<uint32_t> {4});
}

// An odd total has one lowest node. A name on fewer tips than min_count, or on no tip, is
// eligible nowhere. With min_count 1 a name on one tip is eligible at that tip, which no
// internal node can hold alone.
TEST_CASE("LowestEligibleNodes with an odd total, too few tips, or none", "[taxonomy_decorator]") {
	// root(0) +-- s(1) +-- t3(3)  N
	//         |        +-- t4(4)  N
	//         |        +-- t5(5)  (D1; O)
	//         +-- t2(2)  N
	//         +-- h(6) --- t7(7)  (D1; M)
	Lineages lineages = WithN({2, 3, 4});
	lineages.push_back({5, {"D1", "O"}});
	lineages.push_back({7, {"D1", "M"}});
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 1, 0, 6}, lineages);
	CHECK(LowestOf(c, "N", 2) == std::vector<uint32_t> {1});
	CHECK(LowestOf(c, "M", 2).empty());
	CHECK(LowestOf(c, "M", 1) == std::vector<uint32_t> {7});

	// P3 sits on a row no tip carries, so its total is 0.
	auto store = HandStore();
	const NameId p3 = store.Name(4, kPhylum);
	CladeCountIndex index(HandTree(), store, HandRows());
	REQUIRE(index.Total(kPhylum, p3) == 0);
	CHECK(LowestEligibleNodes(index, kPhylum, p3, F1Score::Threshold(0, 1)).empty());
}

// Against the brute-force set, over counts made by walking up from every tip, on random trees
// with chains, multifurcations and 50/50 splits, and min_count 1 to 3.
TEST_CASE("LowestEligibleNodes matches the brute-force lowest eligible nodes", "[taxonomy_decorator]") {
	std::mt19937 rng(20261120);
	bool saw_two = false, saw_tip = false, saw_none = false;
	for (int trial = 0; trial < 300; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		const auto index = build_index(c, store);
		const uint32_t root = root_of(c.tree);
		const auto internal = has_child(c.tree);
		auto by_preorder = [&](uint32_t a, uint32_t b) {
			return index.Tree().pre(a) < index.Tree().pre(b);
		};
		INFO("trial " << trial << ", " << n << " nodes");
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			for (size_t r = 0; r < kRankCount; r++) {
				for (NameId name = 0; name < store.NameCount(r); name++) {
					const uint32_t total = naive.Total(r, name, root);
					const uint32_t threshold = F1Score::Threshold(total, min_count);
					auto expected = naive_lowest_eligible(c, naive, r, name, threshold);
					std::sort(expected.begin(), expected.end(), by_preorder);
					INFO("min_count " << min_count << ", rank " << r << ", name " << name << ", total " << total);
					REQUIRE(LowestEligibleNodes(index, r, name, threshold) == expected);
					saw_two |= expected.size() == 2;
					saw_tip |= expected.size() == 1 && !internal[expected[0]];
					saw_none |= expected.empty();
				}
			}
		}
	}
	// Each kind of result must occur, or this test proves nothing about it.
	CHECK(saw_two);
	CHECK(saw_tip);
	CHECK(saw_none);
}

// tax2tree applies min_count at every rank. Here strain uses 1: a strain is about one genome, and
// under a threshold of 2 it could never be assigned.
TEST_CASE("MinCountAt is min_count, except 1 at strain", "[taxonomy_decorator]") {
	const size_t strain = kRankCount - 1;
	for (size_t r = 0; r < strain; r++) {
		INFO("rank " << r);
		CHECK(MinCountAt(r, 1) == 1);
		CHECK(MinCountAt(r, 2) == 2);
		CHECK(MinCountAt(r, 5) == 5);
	}
	CHECK(MinCountAt(strain, 1) == 1);
	CHECK(MinCountAt(strain, 2) == 1);
	CHECK(MinCountAt(strain, 5) == 1);
}

namespace {
// Masks are compared as unsigned, which Catch prints as numbers; it prints uint8_t as characters.
std::vector<unsigned> Widen(const std::vector<uint8_t> &masks) {
	return std::vector<unsigned>(masks.begin(), masks.end());
}

std::vector<unsigned> SafeMasks(const CladeCountIndex &index, uint32_t min_count) {
	return Widen(RankSafeMasks(index, min_count));
}

// Per-node masks, from (node, mask) pairs; unlisted nodes are 0.
std::vector<unsigned> MasksOf(size_t n, std::initializer_list<std::pair<uint32_t, unsigned>> masks) {
	std::vector<unsigned> out(n, 0);
	for (const auto &[node, mask] : masks) {
		out[node] = mask;
	}
	return out;
}

constexpr unsigned kDomainBit = 1u << kDomain, kPhylumBit = 1u << kPhylum;
} // namespace

// Two names are eligible at v and at c1, so the phylum is unsafe there. Counting each name's
// lowest eligible nodes up the tree gives v 3 (X's two halves and Y), less 1 at X's LCA, so 2.
// A count capped at 2 on the way up would give 2 - 1 = 1 and call v safe.
TEST_CASE("RankSafeMasks counts a 50/50 name and another name exactly", "[taxonomy_decorator]") {
	// v(0) +-- c1(1) +-- xa(3) +-- t6, t7   X
	//      |         +-- ya(4) +-- t8, t9   Y
	//      |         +-- t5                 Y
	//      +-- c2(2) +-- t10, t11           X
	// Phylum: X's lowest nodes are xa and c2 (2 of 4 each), Y's is ya (2 of 3).
	// Domain: D1 is on all 7 tips, so it needs 4, which v and c1 hold.
	Lineages lineages;
	for (uint32_t tip : {6, 7, 10, 11}) {
		lineages.push_back({tip, {"D1", "X"}});
	}
	for (uint32_t tip : {5, 8, 9}) {
		lineages.push_back({tip, {"D1", "Y"}});
	}
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 1, 3, 3, 4, 4, 2, 2}, lineages);
	CHECK(SafeMasks(c.index, 2) ==
	      MasksOf(12, {{0, kDomainBit}, {1, kDomainBit}, {2, kPhylumBit}, {3, kPhylumBit}, {4, kPhylumBit}}));
}

// A name split 50/50 is still one name at, and above, the LCA of its halves: without the
// correction there, root and v would count X twice and call the phylum unsafe. D1 splits 4 of 5,
// so its lowest node is v, and Z, on one tip, is eligible nowhere under min_count 2.
TEST_CASE("RankSafeMasks counts a 50/50 name once at the LCA of its halves", "[taxonomy_decorator]") {
	// root(0) +-- v(1) +-- c1(3) +-- t5, t6   X
	//         |        +-- c2(4) +-- t7, t8   X
	//         +-- t2   Z
	Lineages lineages = {{2, {"D1", "Z"}}};
	for (uint32_t tip : {5, 6, 7, 8}) {
		lineages.push_back({tip, {"D1", "X"}});
	}
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 3, 3, 4, 4}, lineages);
	const unsigned both = kDomainBit | kPhylumBit;
	CHECK(SafeMasks(c.index, 2) == MasksOf(9, {{0, both}, {1, both}, {3, kPhylumBit}, {4, kPhylumBit}}));
}

// Each strain is on one genome, under its own one-child holder. Strain uses a min_count of 1, so
// each holder is safe at strain, while the root holds both strains and is not. Every shallower
// rank needs 2, which only the root holds. A tip, eligible for its strain, is never safe.
TEST_CASE("RankSafeMasks makes a singleton strain safe at its holder", "[taxonomy_decorator]") {
	// root(0) +-- h1(1) --- t3   ... S1; T1
	//         +-- h2(2) --- t4   ... S1; T2
	const std::vector<std::string> shared {"D1", "P1", "C1", "O1", "F1", "G1", "S1"};
	auto lineage = [&](std::string strain) {
		auto out = shared;
		out.push_back(std::move(strain));
		return out;
	};
	auto c = MakeCase({kNoParent, 0, 0, 1, 2}, {{3, lineage("T1")}, {4, lineage("T2")}});
	const unsigned strain_bit = 1u << (kRankCount - 1);
	const unsigned shallower = strain_bit - 1;
	CHECK(SafeMasks(c.index, 2) == MasksOf(5, {{0, shallower}, {1, strain_bit}, {2, strain_bit}}));
	// Under min_count 1 every shared name needs only ceil(2/2) = 1, so each holder is safe at
	// every rank: why min_count stays 2 above strain.
	CHECK(SafeMasks(c.index, 1) == MasksOf(5, {{0, shallower}, {1, 0xFF}, {2, 0xFF}}));
}

// min_count 0 would let a node with no informative tips be eligible, with precision 0/0.
TEST_CASE("RankSafeMasks rejects a min_count of 0", "[taxonomy_decorator]") {
	auto c = MakeCase({kNoParent, 0, 0}, WithN({1, 2}));
	CHECK_THROWS_AS(RankSafeMasks(c.index, 0), std::invalid_argument);
}

// Against brute force: at each internal node, count the names whose naive count reaches their
// threshold, with strain's min_count of 1 written out here rather than taken from MinCountAt.
TEST_CASE("RankSafeMasks matches brute-force rank safety", "[taxonomy_decorator]") {
	std::mt19937 rng(20261121);
	bool saw_safe = false, saw_unsafe = false, saw_split_safe = false;
	for (int trial = 0; trial < 300; trial++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		const auto index = build_index(c, store);
		const uint32_t root = root_of(c.tree);
		const auto internal = has_child(c.tree);
		INFO("trial " << trial << ", " << n << " nodes");
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			std::vector<unsigned> expected(n, 0);
			for (size_t r = 0; r < kRankCount; r++) {
				const uint32_t rank_min_count = r == kRankCount - 1 ? 1 : min_count;
				std::vector<uint32_t> eligible_names(n, 0);
				std::vector<std::vector<uint32_t>> split_lowest;
				for (NameId name = 0; name < store.NameCount(r); name++) {
					const uint32_t threshold = F1Score::Threshold(naive.Total(r, name, root), rank_min_count);
					for (uint32_t v = 0; v < n; v++) {
						eligible_names[v] += naive.clade[r][name][v] >= threshold;
					}
					auto lowest = naive_lowest_eligible(c, naive, r, name, threshold);
					if (lowest.size() == 2) {
						split_lowest.push_back(std::move(lowest));
					}
				}
				for (uint32_t v = 0; v < n; v++) {
					if (internal[v] && eligible_names[v] == 1) {
						expected[v] |= 1u << r;
					}
					saw_safe |= internal[v] && eligible_names[v] == 1;
					saw_unsafe |= internal[v] && eligible_names[v] >= 2;
				}
				// The LCA correction must matter somewhere: a 50/50 name alone at its halves' LCA.
				for (const auto &lowest : split_lowest) {
					const uint32_t lca = naive_lca(c.tree, lowest[0], lowest[1]);
					saw_split_safe |= eligible_names[lca] == 1;
				}
			}
			INFO("min_count " << min_count);
			REQUIRE(SafeMasks(index, min_count) == expected);
		}
	}
	CHECK(saw_safe);
	CHECK(saw_unsafe);
	CHECK(saw_split_safe);
}

// tax2tree's pick_names walks the ranks from domain down, assigning at each safe rank and
// stopping at the first unsafe rank after an assignment; unsafe ranks before the first safe one
// are skipped.
TEST_CASE("PickedRanks takes the run of safe ranks from the shallowest", "[taxonomy_decorator]") {
	auto picked = [](unsigned mask) {
		return unsigned(PickedRanks(static_cast<uint8_t>(mask)));
	};
	CHECK(picked(0) == 0);
	CHECK(picked(0b0000'0110) == 0b0000'0110);
	CHECK(picked(0b0000'1011) == 0b0000'0011);
	CHECK(picked(0b1010'1010) == 0b0000'0010);
	CHECK(picked(0b1111'0000) == 0b1111'0000);
	CHECK(picked(0b1111'1111) == 0b1111'1111);
	CHECK(picked(0b1000'0000) == 0b1000'0000);
	CHECK(picked(0b0000'0001) == 0b0000'0001);

	// Every mask, against pick_names' loop written out.
	for (unsigned mask = 0; mask < 256; mask++) {
		unsigned expected = 0;
		bool assigned = false;
		for (size_t r = 0; r < kRankCount; r++) {
			if (mask >> r & 1) {
				expected |= 1u << r;
				assigned = true;
			} else if (assigned) {
				break;
			}
		}
		INFO("mask " << mask);
		REQUIRE(picked(mask) == expected);
	}
}

namespace {
std::vector<NameScore> NamesAt(const NameAssignment &assignment, uint32_t node) {
	const auto names = assignment.NamesAt(node);
	return std::vector<NameScore>(names.begin(), names.end());
}

// Per-node names, from (node, names) pairs; unlisted nodes have none.
std::vector<std::vector<NameScore>> NamesOf(size_t n,
                                            std::initializer_list<std::pair<uint32_t, std::vector<NameScore>>> names) {
	std::vector<std::vector<NameScore>> out(n);
	for (const auto &[node, scores] : names) {
		out[node] = scores;
	}
	return out;
}

std::vector<std::vector<NameScore>> AllNames(const NameAssignment &assignment, size_t n) {
	std::vector<std::vector<NameScore>> out(n);
	for (uint32_t v = 0; v < n; v++) {
		out[v] = NamesAt(assignment, v);
	}
	return out;
}

std::vector<NotPlacedName> ByRankAndName(std::vector<NotPlacedName> names) {
	std::sort(names.begin(), names.end(), [](const NotPlacedName &a, const NotPlacedName &b) {
		return std::tie(a.rank, a.name) < std::tie(b.rank, b.name);
	});
	return names;
}
} // namespace

// N is 2 of 2 tips at both u and its one child w, so they tie on every count: the shallower u,
// later in postorder, keeps it. The root also holds N but scores lower (2 of 3 tips), while for
// D1 the root scores higher than u or w (3 of 3 against 2 of 2 of a total of 3): a strict win
// each way. Z, on one tip, is eligible nowhere under min_count 2.
TEST_CASE("AssignNames breaks a unary-chain tie toward the shallower node", "[taxonomy_decorator]") {
	// root(0) +-- u(1) --- w(2) +-- t3  N
	//         |                 +-- t4  N
	//         +-- t5  Z
	Lineages lineages = WithN({3, 4});
	lineages.push_back({5, {"D1", "Z"}});
	auto c = MakeCase({kNoParent, 0, 1, 2, 2, 0}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), n = c.Name(kPhylum, "N"), z = c.Name(kPhylum, "Z");
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 6) == NamesOf(6, {{0, {{kDomain, d1, 0, 3, 3, 3}}}, {1, {{kPhylum, n, 1, 2, 2, 2}}}}));
	CHECK(assignment.NotPlaced() == std::vector<NotPlacedName> {{kPhylum, z, 1, NotPlacedReason::kNeverEligible}});
}

// X splits 2 and 2 between a and b, and Y makes v unsafe at the phylum, so X's candidates are a,
// its one-child parent p, and b, all at 2 of 2 tips: the latest in postorder, b, keeps it. b's
// index is below a's, so this is not index order.
TEST_CASE("AssignNames breaks a 50/50 tie toward the node later in postorder", "[taxonomy_decorator]") {
	// v(0) +-- p(1) --- a(9) +-- t10, t11   X
	//      +-- b(2) +-- t4, t5              X
	//      +-- c(3) +-- t6, t7, t8          Y
	Lineages lineages;
	for (uint32_t tip : {4, 5, 10, 11}) {
		lineages.push_back({tip, {"D1", "X"}});
	}
	for (uint32_t tip : {6, 7, 8}) {
		lineages.push_back({tip, {"D1", "Y"}});
	}
	auto c = MakeCase({kNoParent, 0, 0, 0, 2, 2, 3, 3, 3, 1, 9, 9}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), x = c.Name(kPhylum, "X"), y = c.Name(kPhylum, "Y");
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 12) == NamesOf(12, {{0, {{kDomain, d1, 0, 7, 7, 7}}},
	                                               {2, {{kPhylum, x, 2, 2, 2, 4}}},
	                                               {3, {{kPhylum, y, 3, 3, 3, 3}}}}));
	CHECK(assignment.NotPlaced().empty());
}

// The tie-breaks apply only to equal scores. X splits 2 and 2 between a and b, but b also holds
// an informative Z tip, so a's F1 (4/6) beats b's (4/7) though b is later in postorder.
TEST_CASE("AssignNames prefers a higher F1 to a later node", "[taxonomy_decorator]") {
	// v(0) +-- a(1) +-- t4, t5             X
	//      +-- b(2) +-- t6, t7  X;  t8  Z
	//      +-- c(3) +-- t9, t10, t11       Y
	Lineages lineages;
	for (uint32_t tip : {4, 5, 6, 7}) {
		lineages.push_back({tip, {"D1", "X"}});
	}
	lineages.push_back({8, {"D1", "Z"}});
	for (uint32_t tip : {9, 10, 11}) {
		lineages.push_back({tip, {"D1", "Y"}});
	}
	auto c = MakeCase({kNoParent, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), x = c.Name(kPhylum, "X"), y = c.Name(kPhylum, "Y");
	const NameId z = c.Name(kPhylum, "Z");
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 12) == NamesOf(12, {{0, {{kDomain, d1, 0, 8, 8, 8}}},
	                                               {1, {{kPhylum, x, 1, 2, 2, 4}}},
	                                               {3, {{kPhylum, y, 3, 3, 3, 3}}}}));
	CHECK(assignment.NotPlaced() == std::vector<NotPlacedName> {{kPhylum, z, 1, NotPlacedReason::kNeverEligible}});
}

// With equal F1, the node with fewer tips keeps the name, counting every tip: b holds the same 2
// of 2 informative tips as a, plus a tip with no row, so a keeps X though b is later in postorder.
TEST_CASE("AssignNames breaks an F1 tie toward fewer tips", "[taxonomy_decorator]") {
	// v(0) +-- a(1) +-- t4, t5             X
	//      +-- b(2) +-- t6, t7  X;  t8  no row
	//      +-- c(3) +-- t9, t10, t11       Y
	Lineages lineages;
	for (uint32_t tip : {4, 5, 6, 7}) {
		lineages.push_back({tip, {"D1", "X"}});
	}
	for (uint32_t tip : {9, 10, 11}) {
		lineages.push_back({tip, {"D1", "Y"}});
	}
	auto c = MakeCase({kNoParent, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3, 3}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), x = c.Name(kPhylum, "X"), y = c.Name(kPhylum, "Y");
	REQUIRE(c.index.TipCount(2) == 3);
	REQUIRE(c.index.InformativeTipCount(2) == 2);
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 12) == NamesOf(12, {{0, {{kDomain, d1, 0, 7, 7, 7}}},
	                                               {1, {{kPhylum, x, 1, 2, 2, 4}}},
	                                               {3, {{kPhylum, y, 3, 3, 3, 3}}}}));
}

// A node's preliminary name is its deepest surviving one, not its deepest picked one: the root
// picks domain, phylum and class, but C1 scores higher at q (2 of 2) than at the root (2 of 3),
// so the root's preliminary name is P1.
TEST_CASE("AssignNames gives the deepest surviving name as preliminary", "[taxonomy_decorator]") {
	// root(0) +-- q(1) +-- t3, t4  (D1; P1; C1)
	//         +-- t2   (D1; P1; C2)
	auto c =
	    MakeCase({kNoParent, 0, 0, 1, 1}, {{2, {"D1", "P1", "C2"}}, {3, {"D1", "P1", "C1"}}, {4, {"D1", "P1", "C1"}}});
	const NameId d1 = c.Name(kDomain, "D1"), p1 = c.Name(kPhylum, "P1");
	const NameId c1 = c.Name(kClass, "C1"), c2 = c.Name(kClass, "C2");
	const auto assignment = AssignNames(c.index, 2);
	const NameScore root_p1 {kPhylum, p1, 0, 3, 3, 3}, q_c1 {kClass, c1, 1, 2, 2, 2};
	CHECK(AllNames(assignment, 5) == NamesOf(5, {{0, {{kDomain, d1, 0, 3, 3, 3}, root_p1}}, {1, {q_c1}}}));
	REQUIRE(assignment.Preliminary(0) != nullptr);
	CHECK(*assignment.Preliminary(0) == root_p1);
	REQUIRE(assignment.Preliminary(1) != nullptr);
	CHECK(*assignment.Preliminary(1) == q_c1);
	CHECK(assignment.Preliminary(2) == nullptr);
	CHECK(assignment.NotPlaced() == std::vector<NotPlacedName> {{kClass, c2, 1, NotPlacedReason::kNeverEligible}});
}

// A name with tips that goes nowhere is reported with why. On this star, PA and PB are both
// eligible only at the root, so neither is ever the one eligible phylum; PC is on one tip; and C1
// is safe at the root, but the root's run of ranks stops at the unsafe phylum.
TEST_CASE("AssignNames reports why each unassigned name went nowhere", "[taxonomy_decorator]") {
	// root(0) --- t1, t2, t3 (D1; PA; C1);  t4, t5 (D1; PB);  t6 (D1; PC)
	Lineages lineages;
	for (uint32_t tip : {1, 2, 3}) {
		lineages.push_back({tip, {"D1", "PA", "C1"}});
	}
	lineages.push_back({4, {"D1", "PB"}});
	lineages.push_back({5, {"D1", "PB"}});
	lineages.push_back({6, {"D1", "PC"}});
	auto c = MakeCase({kNoParent, 0, 0, 0, 0, 0, 0}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), c1 = c.Name(kClass, "C1");
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 7) == NamesOf(7, {{0, {{kDomain, d1, 0, 6, 6, 6}}}}));
	CHECK(assignment.NotPlaced() == ByRankAndName({{kPhylum, c.Name(kPhylum, "PA"), 3, NotPlacedReason::kNeverRankSafe},
	                                               {kPhylum, c.Name(kPhylum, "PB"), 2, NotPlacedReason::kNeverRankSafe},
	                                               {kPhylum, c.Name(kPhylum, "PC"), 1, NotPlacedReason::kNeverEligible},
	                                               {kClass, c1, 3, NotPlacedReason::kBlockedByContiguity}}));
}

// min_count 0 would let a node with no informative tips be eligible, with precision 0/0.
TEST_CASE("AssignNames rejects a min_count of 0", "[taxonomy_decorator]") {
	auto c = MakeCase({kNoParent, 0, 0}, WithN({1, 2}));
	CHECK_THROWS_AS(AssignNames(c.index, 0), std::invalid_argument);
}

// Two genomes of species S1, each under a one-child holder, beside two of S2. Strain uses a
// min_count of 1, so each strain lands on its holder. Under min_count 2 the species needs both
// genomes, and the root holds both species, so S1 goes nowhere for either genome: both are
// treated the same. Under min_count 1 each holder is a candidate for S1 at an equal score, and
// the tie gives S1 to h2 alone, so h1's genome loses its species.
TEST_CASE("AssignNames puts singleton strains on their holders, the same for both genomes", "[taxonomy_decorator]") {
	// root(0) +-- h1(1) --- t3  ...; G1; S1; T1
	//         +-- h2(2) --- t4  ...; G1; S1; T2
	//         +-- t5, t6        ...; G1; S2
	const std::vector<std::string> shared {"D1", "P1", "C1", "O1", "F1", "G1"};
	auto lineage = [&](std::vector<std::string> deeper) {
		auto out = shared;
		out.insert(out.end(), deeper.begin(), deeper.end());
		return out;
	};
	auto c =
	    MakeCase({kNoParent, 0, 0, 1, 2, 0, 0},
	             {{3, lineage({"S1", "T1"})}, {4, lineage({"S1", "T2"})}, {5, lineage({"S2"})}, {6, lineage({"S2"})}});
	const size_t species = kRankCount - 2, strain = kRankCount - 1;
	const NameId s1 = c.Name(species, "S1"), s2 = c.Name(species, "S2");
	const NameId t1 = c.Name(strain, "T1"), t2 = c.Name(strain, "T2");
	std::vector<NameScore> root_names;
	for (size_t r = 0; r < shared.size(); r++) {
		root_names.push_back({r, c.Name(r, shared[r]), 0, 4, 4, 4});
	}

	const auto at_two = AssignNames(c.index, 2);
	CHECK(AllNames(at_two, 7) ==
	      NamesOf(7, {{0, root_names}, {1, {{strain, t1, 1, 1, 1, 1}}}, {2, {{strain, t2, 2, 1, 1, 1}}}}));
	CHECK(at_two.NotPlaced() == ByRankAndName({{species, s1, 2, NotPlacedReason::kNeverRankSafe},
	                                           {species, s2, 2, NotPlacedReason::kNeverRankSafe}}));

	const auto at_one = AssignNames(c.index, 1);
	CHECK(AllNames(at_one, 7) == NamesOf(7, {{0, root_names},
	                                         {1, {{strain, t1, 1, 1, 1, 1}}},
	                                         {2, {{species, s1, 2, 1, 1, 2}, {strain, t2, 2, 1, 1, 1}}}}));
	CHECK(at_one.NotPlaced() == std::vector<NotPlacedName> {{species, s2, 2, NotPlacedReason::kNeverRankSafe}});
}

// Against brute force (naive_assign): eligibility, safety, contiguity and the fold applied at
// every node directly, on random trees with chains, multifurcations, 50/50 splits and strains.
// The brute force runs first over every case, so the coverage it requires is checked even
// while AssignNames fails.
TEST_CASE("AssignNames matches the brute-force assignment", "[taxonomy_decorator]") {
	std::mt19937 rng(20261122);
	struct Trial {
		TaxonomyCase c;
		uint32_t min_count;
		NaiveAssignment expected;
	};
	std::vector<Trial> trials;
	size_t won_on_tip_count = 0, won_on_postorder = 0;
	std::array<size_t, 3> reasons {};
	for (int k = 0; k < 300; k++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			auto expected = naive_assign(c, store, naive, min_count);
			won_on_tip_count += expected.won_on_tip_count;
			won_on_postorder += expected.won_on_postorder;
			for (const auto &name : expected.not_placed) {
				reasons[static_cast<size_t>(name.reason)]++;
			}
			trials.push_back({c, min_count, std::move(expected)});
		}
	}
	// Each tie-break and each reason must occur, or this test proves nothing about it.
	CHECK(won_on_tip_count > 0);
	CHECK(won_on_postorder > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kNeverEligible)] > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kNeverRankSafe)] > 0);
	CHECK(reasons[static_cast<size_t>(NotPlacedReason::kBlockedByContiguity)] > 0);

	for (size_t k = 0; k < trials.size(); k++) {
		const auto &[c, min_count, expected] = trials[k];
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		const uint32_t n = static_cast<uint32_t>(c.tree.parent.size());
		INFO("trial " << k / 3 << ", " << n << " nodes, min_count " << min_count);
		const auto assignment = AssignNames(index, min_count);
		REQUIRE(assignment.Scores() == expected.scores);
		REQUIRE(assignment.NotPlaced() == expected.not_placed);
		std::vector<std::vector<NameScore>> by_node(n);
		for (const auto &score : expected.scores) {
			by_node[score.node].push_back(score);
		}
		for (uint32_t v = 0; v < n; v++) {
			INFO("node " << v);
			std::sort(by_node[v].begin(), by_node[v].end(),
			          [](const NameScore &a, const NameScore &b) { return a.rank < b.rank; });
			REQUIRE(NamesAt(assignment, v) == by_node[v]);
			const NameScore *preliminary = assignment.Preliminary(v);
			if (by_node[v].empty()) {
				REQUIRE(preliminary == nullptr);
			} else {
				REQUIRE(preliminary != nullptr);
				REQUIRE(*preliminary == by_node[v].back());
			}
		}
	}
}
