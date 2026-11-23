#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "taxonomy_decorator_fixtures.hpp"

#include <algorithm>
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
enum Rank : size_t { kDomain, kPhylum, kClass, kOrder, kFamily, kGenus, kSpecies, kStrain };

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
	CHECK(index.NameCount(kStrainRank) == 0);
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

// Gives each of `tips` the lineage `lineage`.
void Add(Lineages &lineages, const std::vector<std::string> &lineage, std::initializer_list<uint32_t> tips) {
	for (uint32_t tip : tips) {
		lineages.push_back({tip, lineage});
	}
}

// Each of `tips` carries the lineage (D1; N), so N is a phylum.
Lineages WithN(std::initializer_list<uint32_t> tips) {
	Lineages out;
	Add(out, {"D1", "N"}, tips);
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

	// An odd total has one middle tip and one walk, so a first walk from a fails here even if an
	// even total's second walk, from the upper middle tip, would have found B.
	// root(0) +-- a(1)
	//         +-- B(2) --- tips 3..12
	std::vector<uint32_t> odd(13, 2);
	odd[0] = kNoParent;
	odd[1] = 0;
	odd[2] = 0;
	auto odd_case = MakeCase(odd, WithN({1, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12}));
	REQUIRE(odd_case.index.Total(kPhylum, odd_case.Name(kPhylum, "N")) == 11);
	CHECK(LowestOf(odd_case, "N", 2) == std::vector<uint32_t> {2});
}

// The search walks up from the middle tips only, which finds every lowest eligible node only
// when the threshold is at least 1 and at least half the name's total. Below that it would
// return wrong nodes, so it throws instead.
TEST_CASE("LowestEligibleNodes rejects a threshold below 1 or below half the total", "[taxonomy_decorator]") {
	// root(0) +-- x(1) +-- t3, t4  N
	//         +-- y(2) +-- t5, t6  N
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 2, 2}, WithN({3, 4, 5, 6}));
	const NameId n = c.Name(kPhylum, "N");
	REQUIRE(c.index.Total(kPhylum, n) == 4);
	CHECK_THROWS_AS(LowestEligibleNodes(c.index, kPhylum, n, 0), std::invalid_argument);
	CHECK_THROWS_AS(LowestEligibleNodes(c.index, kPhylum, n, 1), std::invalid_argument);
	CHECK(LowestEligibleNodes(c.index, kPhylum, n, 2) == std::vector<uint32_t> {1, 2});
	// With an odd total, half rounds up: 5 tips need 3.
	// root(0) +-- x(1) +-- t2, t3, t4  N
	//         +-- t5, t6  N
	auto odd = MakeCase({kNoParent, 0, 1, 1, 1, 0, 0}, WithN({2, 3, 4, 5, 6}));
	const NameId odd_n = odd.Name(kPhylum, "N");
	REQUIRE(odd.index.Total(kPhylum, odd_n) == 5);
	CHECK_THROWS_AS(LowestEligibleNodes(odd.index, kPhylum, odd_n, 2), std::invalid_argument);
	CHECK(LowestEligibleNodes(odd.index, kPhylum, odd_n, 3) == std::vector<uint32_t> {1});
	// A name on no tip still needs a threshold of at least 1.
	auto store = HandStore();
	CladeCountIndex index(HandTree(), store, HandRows());
	CHECK_THROWS_AS(LowestEligibleNodes(index, kPhylum, store.Name(4, kPhylum), 0), std::invalid_argument);
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
// eligible nowhere. With min_count 1 a name on one tip is eligible at that tip, and the tip is
// its lowest node, not its one-child parent h, which holds it too.
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

// tax2tree applies min_count at every rank. Here strain uses 1 (MinCountAt says why).
TEST_CASE("MinCountAt is min_count, except 1 at strain", "[taxonomy_decorator]") {
	const size_t strain = kStrainRank;
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

// Per-node values, from (node, value) pairs; unlisted nodes get T's default.
template <class T>
std::vector<T> PerNode(size_t n, std::initializer_list<std::pair<uint32_t, T>> values) {
	std::vector<T> out(n);
	for (const auto &[node, value] : values) {
		out[node] = value;
	}
	return out;
}

// Per-node masks; unlisted nodes are 0.
std::vector<unsigned> MasksOf(size_t n, std::initializer_list<std::pair<uint32_t, unsigned>> masks) {
	return PerNode(n, masks);
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
	Add(lineages, {"D1", "X"}, {6, 7, 10, 11});
	Add(lineages, {"D1", "Y"}, {5, 8, 9});
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
	Add(lineages, {"D1", "X"}, {5, 6, 7, 8});
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
	const unsigned strain_bit = 1u << kStrainRank;
	const unsigned shallower = strain_bit - 1;
	CHECK(SafeMasks(c.index, 2) == MasksOf(5, {{0, shallower}, {1, strain_bit}, {2, strain_bit}}));
	// Under min_count 1 every shared name needs only ceil(2/2) = 1, so each holder is safe at
	// every rank: why only strain uses 1.
	CHECK(SafeMasks(c.index, 1) == MasksOf(5, {{0, shallower}, {1, 0xFF}, {2, 0xFF}}));
}

// F1Score::Threshold says why min_count 0 is rejected.
TEST_CASE("RankSafeMasks rejects a min_count of 0", "[taxonomy_decorator]") {
	auto c = MakeCase({kNoParent, 0, 0}, WithN({1, 2}));
	CHECK_THROWS_AS(RankSafeMasks(c.index, 0), std::invalid_argument);
}

// Against brute force (naive_eligible): at each internal node, count the names whose naive count
// reaches their threshold, with strain's min_count of 1 written out rather than taken from
// MinCountAt.
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
			const auto eligible = naive_eligible(c, store, naive, min_count);
			std::vector<unsigned> expected(n, 0);
			for (size_t r = 0; r < kRankCount; r++) {
				const auto &eligible_names = eligible.count[r];
				std::vector<std::vector<uint32_t>> split_lowest;
				for (NameId name = 0; name < store.NameCount(r); name++) {
					const uint32_t threshold = naive_threshold(r, naive.Total(r, name, root), min_count);
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

	// Every mask, against pick_names' loop written out (naive_picked).
	for (unsigned mask = 0; mask < 256; mask++) {
		INFO("mask " << mask);
		REQUIRE(picked(mask) == naive_picked(mask));
	}
}

namespace {
std::vector<NameScore> NamesAt(const NameAssignment &assignment, uint32_t node) {
	const auto names = assignment.NamesAt(node);
	return std::vector<NameScore>(names.begin(), names.end());
}

// Per-node names; unlisted nodes have none.
std::vector<std::vector<NameScore>> NamesOf(size_t n,
                                            std::initializer_list<std::pair<uint32_t, std::vector<NameScore>>> names) {
	return PerNode(n, names);
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
	Add(lineages, {"D1", "X"}, {4, 5, 10, 11});
	Add(lineages, {"D1", "Y"}, {6, 7, 8});
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
	Add(lineages, {"D1", "X"}, {4, 5, 6, 7});
	lineages.push_back({8, {"D1", "Z"}});
	Add(lineages, {"D1", "Y"}, {9, 10, 11});
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
	Add(lineages, {"D1", "X"}, {4, 5, 6, 7});
	Add(lineages, {"D1", "Y"}, {9, 10, 11});
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
	Add(lineages, {"D1", "PA", "C1"}, {1, 2, 3});
	Add(lineages, {"D1", "PB"}, {4, 5});
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

// A node where a name's rank is safe but not picked is not the end of the name's walk. At v the
// class is safe (C1 alone is eligible) but v picks only the domain: no phylum is eligible there
// yet, since P1 has 2 of its 5 tips at v, and Pa and Pb one each. At u, P1 is eligible, so u
// picks domain, phylum and class, and C1 is assigned there.
TEST_CASE("AssignNames looks past a node that skips the name's rank", "[taxonomy_decorator]") {
	// u(0) +-- v(1) +-- t2, t3  (D1; P1; C1)
	//      |        +-- t4  (D1; Pa);  t5  (D1; Pb)
	//      +-- t6, t7, t8  (D1; P1)
	Lineages lineages;
	Add(lineages, {"D1", "P1", "C1"}, {2, 3});
	Add(lineages, {"D1", "Pa"}, {4});
	Add(lineages, {"D1", "Pb"}, {5});
	Add(lineages, {"D1", "P1"}, {6, 7, 8});
	auto c = MakeCase({kNoParent, 0, 1, 1, 1, 1, 0, 0, 0}, lineages);
	const NameId d1 = c.Name(kDomain, "D1"), p1 = c.Name(kPhylum, "P1"), c1 = c.Name(kClass, "C1");
	REQUIRE(SafeMasks(c.index, 2)[1] == (kDomainBit | 1u << kClass));
	const auto assignment = AssignNames(c.index, 2);
	CHECK(AllNames(assignment, 9) ==
	      NamesOf(9, {{0, {{kDomain, d1, 0, 7, 7, 7}, {kPhylum, p1, 0, 5, 7, 5}, {kClass, c1, 0, 2, 7, 2}}}}));
	CHECK(assignment.NotPlaced() ==
	      ByRankAndName({{kPhylum, c.Name(kPhylum, "Pa"), 1, NotPlacedReason::kNeverEligible},
	                     {kPhylum, c.Name(kPhylum, "Pb"), 1, NotPlacedReason::kNeverEligible}}));
}

// F1Score::Threshold says why min_count 0 is rejected.
TEST_CASE("AssignNames rejects a min_count of 0", "[taxonomy_decorator]") {
	auto c = MakeCase({kNoParent, 0, 0}, WithN({1, 2}));
	CHECK_THROWS_AS(AssignNames(c.index, 0), std::invalid_argument);
}

// Two genomes of species S1, each under a one-child holder, beside two of S2. Strain uses a
// min_count of 1, so each strain lands on its holder. Under min_count 2 the species needs both
// genomes, and the root holds both species, so S1 goes nowhere for either genome: both are
// treated the same. Under min_count 1 each holder is a candidate for S1 at an equal score, and
// the tie gives S1 to h2 alone, so h1's genome loses its species.
TEST_CASE("AssignNames puts singleton strains on their holders, and under min_count 2 treats both genomes alike",
          "[taxonomy_decorator]") {
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
	const size_t species = kSpeciesRank, strain = kStrainRank;
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
	NaiveCoverage coverage;
	for (int k = 0; k < 300; k++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto naive = naive_counts(c, store);
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			auto expected = naive_assign(c, store, naive, min_count);
			coverage.Add(expected);
			trials.push_back({c, min_count, std::move(expected)});
		}
	}
	CHECK(coverage.Missing().empty());

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

namespace {
// A label as text, "p__P1; c__C1", in the order it lists its names. The stage tests read labels
// through this and LabelFrom rather than FormatLabel and ParseNodeLabel, so they don't depend on
// the code the label tests check.
std::string LabelText(const NodeLabel &label, const TaxonomyStore &store) {
	std::string out;
	for (const auto &[rank, name] : label) {
		if (!out.empty()) {
			out += "; ";
		}
		out += kLineageRanks[rank].letter;
		out += "__" + store.NameText(rank, name);
	}
	return out;
}

std::vector<std::string> LabelTexts(const std::vector<NodeLabel> &labels, const TaxonomyStore &store) {
	std::vector<std::string> out;
	for (const auto &label : labels) {
		out.push_back(LabelText(label, store));
	}
	return out;
}

std::vector<std::string> LabelTexts(const Backfill &backfill, const TaxonomyStore &store) {
	return LabelTexts(backfill.labels, store);
}

// Per-node label text; unlisted nodes have an empty label.
std::vector<std::string> TextsOf(size_t n, std::initializer_list<std::pair<uint32_t, std::string>> texts) {
	return PerNode(n, texts);
}

// The names the fold gave a node, as label text.
std::string FoldText(const NameAssignment &assignment, uint32_t node, const TaxonomyStore &store) {
	NodeLabel label;
	for (const auto &score : assignment.NamesAt(node)) {
		label.push_back({score.rank, score.name});
	}
	return LabelText(label, store);
}
} // namespace

// q's fold left it only O1: its phylum and class are unsafe there, since P1 and C1 have most of
// their tips under r. Its nearest named ancestor, the root, is named at the domain, so the
// phylum and class between them come from the hierarchy, O1's parent C1 and C1's parent P1. s
// is named one rank below the root, so there is no gap and it keeps P2 alone. Tips, and
// internal nodes without a name, get empty labels.
TEST_CASE("BackfillNames fills the ranks between a node and its named ancestor from the hierarchy",
          "[taxonomy_decorator]") {
	// root(0) +-- q(1) +-- t4, t5      D1; P1; C1; O1
	//         +-- r(2) +-- t6, t7, t8  D1; P1; C1; O2
	//         +-- s(3) +-- t9, t10     D1; P2
	Lineages lineages;
	Add(lineages, {"D1", "P1", "C1", "O1"}, {4, 5});
	Add(lineages, {"D1", "P1", "C1", "O2"}, {6, 7, 8});
	Add(lineages, {"D1", "P2"}, {9, 10});
	auto c = MakeCase({kNoParent, 0, 0, 0, 1, 1, 2, 2, 2, 3, 3}, lineages);
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store) == "d__D1");
	REQUIRE(FoldText(assignment, 1, c.store) == "o__O1");
	REQUIRE(FoldText(assignment, 2, c.store) == "p__P1; c__C1; o__O2");
	REQUIRE(FoldText(assignment, 3, c.store) == "p__P2");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) ==
	      TextsOf(11, {{0, "d__D1"}, {1, "p__P1; c__C1; o__O1"}, {2, "p__P1; c__C1; o__O2"}, {3, "p__P2"}}));
	CHECK(backfill.inverted.empty());
}

// tax2tree keeps only a node's deepest name and rebuilds the shallower ones from the hierarchy,
// so a fold name that disagrees with the hierarchy is replaced. n holds Py, the one phylum
// eligible there, and Cx, whose parent is Pz: Pz is not eligible at n, having most of its tips
// under m. n's label is Cx's lineage, Pz then Cx.
TEST_CASE("BackfillNames replaces the fold's shallower names with the hierarchy's", "[taxonomy_decorator]") {
	// root(0) +-- n(1) +-- t3, t4            D1; Pz; Cx
	//         |        +-- t5, t6, t7        D1; Py
	//         +-- m(2) +-- t8, t9, t10, t11  D1; Pz; C2
	Lineages lineages;
	Add(lineages, {"D1", "Pz", "Cx"}, {3, 4});
	Add(lineages, {"D1", "Py"}, {5, 6, 7});
	Add(lineages, {"D1", "Pz", "C2"}, {8, 9, 10, 11});
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2}, lineages);
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store) == "d__D1");
	REQUIRE(FoldText(assignment, 1, c.store) == "p__Py; c__Cx");
	REQUIRE(FoldText(assignment, 2, c.store) == "p__Pz; c__C2");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) == TextsOf(12, {{0, "d__D1"}, {1, "p__Pz; c__Cx"}, {2, "p__Pz; c__C2"}}));
	CHECK(backfill.inverted.empty());
}

// The root has no named ancestor, so its fill reaches the domain. tax2tree's backfill skips the
// root, so the root keeps its deepest name alone: D1, which the fold also gave it, is dropped
// and not rebuilt, and no lineage on this tree would carry the domain. q is named one rank below
// the root's P1, so it has no gap to fill.
TEST_CASE("BackfillNames fills the root up to the domain", "[taxonomy_decorator]") {
	// root(0) +-- q(1) +-- t3, t4  (D1; P1; C1)
	//         +-- t2   (D1; P1; C2)
	auto c =
	    MakeCase({kNoParent, 0, 0, 1, 1}, {{2, {"D1", "P1", "C2"}}, {3, {"D1", "P1", "C1"}}, {4, {"D1", "P1", "C1"}}});
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store) == "d__D1; p__P1");
	REQUIRE(FoldText(assignment, 1, c.store) == "c__C1");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) == TextsOf(5, {{0, "d__D1; p__P1"}, {1, "c__C1"}}));
	CHECK(backfill.inverted.empty());
}

// On a two-domain tree the root is unnamed, every rank having two eligible names there, so y,
// which the fold named at every rank from A to S3, has no named ancestor. tax2tree skips such a
// node and keeps only S3, losing the six names above it. Here the fill reaches the domain. x,
// named at the domain alone, has nothing to fill; ab and cd fill up to x's domain.
TEST_CASE("BackfillNames fills a node with no named ancestor up to the domain", "[taxonomy_decorator]") {
	// root(0) +-- x(1) +-- ab(3) +-- t5, t6  B; P1; C1; O1; F1; G1; S1
	//         |        +-- cd(4) +-- t7, t8  B; P2; C2; O2; F2; G2; S2
	//         +-- y(2) +-- t9, t10           A; P3; C3; O3; F3; G3; S3
	auto lineage = [](std::string domain, std::string k) {
		return std::vector<std::string> {domain, "P" + k, "C" + k, "O" + k, "F" + k, "G" + k, "S" + k};
	};
	Lineages lineages;
	Add(lineages, lineage("B", "1"), {5, 6});
	Add(lineages, lineage("B", "2"), {7, 8});
	Add(lineages, lineage("A", "3"), {9, 10});
	auto c = MakeCase({kNoParent, 0, 0, 1, 1, 3, 3, 4, 4, 2, 2}, lineages);
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store).empty());
	REQUIRE(FoldText(assignment, 1, c.store) == "d__B");
	REQUIRE(FoldText(assignment, 2, c.store) == "d__A; p__P3; c__C3; o__O3; f__F3; g__G3; s__S3");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) == TextsOf(11, {{1, "d__B"},
	                                                    {2, "d__A; p__P3; c__C3; o__O3; f__F3; g__G3; s__S3"},
	                                                    {3, "p__P1; c__C1; o__O1; f__F1; g__G1; s__S1"},
	                                                    {4, "p__P2; c__C2; o__O2; f__F2; g__G2; s__S2"}}));
	CHECK(backfill.inverted.empty());
}

// n is named at the family, below a, named at the genus: n won F2 from a, which holds G's tips
// too, while a kept G, eligible nowhere else. With its ancestor's rank deeper than its own,
// tax2tree empties n's label, losing F2; here n keeps F2 and the inversion is recorded. a fills
// the family between it and the root's O1 from the hierarchy: G's parent F1, not the F2 below it.
TEST_CASE("BackfillNames keeps an inverted-rank node's name and records it", "[taxonomy_decorator]") {
	// root(0) +-- a(1) +-- n(2) +-- t4, t5  D1; P1; C1; O1; F2
	//         |        +-- t6, t7           D1; P1; C1; O1; F1; G
	//         +-- m(3) +-- t8 .. t12        D1; P1; C1; O1; F1
	Lineages lineages;
	Add(lineages, {"D1", "P1", "C1", "O1", "F2"}, {4, 5});
	Add(lineages, {"D1", "P1", "C1", "O1", "F1", "G"}, {6, 7});
	Add(lineages, {"D1", "P1", "C1", "O1", "F1"}, {8, 9, 10, 11, 12});
	auto c = MakeCase({kNoParent, 0, 1, 0, 2, 2, 1, 1, 3, 3, 3, 3, 3}, lineages);
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store) == "d__D1; p__P1; c__C1; o__O1");
	REQUIRE(FoldText(assignment, 1, c.store) == "g__G");
	REQUIRE(FoldText(assignment, 2, c.store) == "f__F2");
	REQUIRE(FoldText(assignment, 3, c.store) == "f__F1");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) ==
	      TextsOf(13, {{0, "d__D1; p__P1; c__C1; o__O1"}, {1, "f__F1; g__G"}, {2, "f__F2"}, {3, "f__F1"}}));
	CHECK(backfill.inverted == std::vector<InvertedBackfill> {{2, kFamily, c.Name(kFamily, "F2"), 1, kGenus}});
}

// An inversion at the domain is recorded too. tax2tree never compares a domain's rank with its
// ancestor's, so it keeps the name without a word; the label is the same here, but a domain below
// a phylum is as inverted as a family below a genus. Two domains are eligible at the root, so it
// is named at the phylum, P; n holds D's two tips and wins D, and its phylum is unsafe, Pa and Pb
// being on one tip each.
TEST_CASE("BackfillNames records an inverted node named at the domain", "[taxonomy_decorator]") {
	// root(0) +-- n(1) +-- t2  D; Pa
	//         |        +-- t3  D; Pb
	//         +-- t4, t5       D2; P
	Lineages lineages {{2, {"D", "Pa"}}, {3, {"D", "Pb"}}};
	Add(lineages, {"D2", "P"}, {4, 5});
	auto c = MakeCase({kNoParent, 0, 1, 1, 0, 0}, lineages);
	const auto assignment = AssignNames(c.index, 2);
	REQUIRE(FoldText(assignment, 0, c.store) == "p__P");
	REQUIRE(FoldText(assignment, 1, c.store) == "d__D");
	const auto backfill = BackfillNames(assignment, c.store, c.index.Tree());
	CHECK(LabelTexts(backfill, c.store) == TextsOf(6, {{0, "d__D2; p__P"}, {1, "d__D"}}));
	CHECK(backfill.inverted == std::vector<InvertedBackfill> {{1, kDomain, c.Name(kDomain, "D"), 0, kPhylum}});
}

// Against brute force (naive_backfill): nearest named ancestors found by walking parents, fills
// by walking Parent(), on random trees through AssignNames at min_count 1 to 3. The brute force
// runs first over every case, so the coverage it requires is checked even while BackfillNames
// fails.
TEST_CASE("BackfillNames matches the brute-force backfill", "[taxonomy_decorator]") {
	std::mt19937 rng(20261120);
	struct Trial {
		TaxonomyCase c;
		uint32_t min_count;
		NaiveBackfill expected;
	};
	std::vector<Trial> trials;
	BackfillCoverage coverage;
	for (int k = 0; k < 300; k++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			auto expected = naive_backfill(c, store, AssignNames(index, min_count));
			coverage.Add(expected);
			trials.push_back({c, min_count, std::move(expected)});
		}
	}
	CHECK(coverage.Missing().empty());

	for (size_t k = 0; k < trials.size(); k++) {
		const auto &[c, min_count, expected] = trials[k];
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		INFO("trial " << k / 3 << ", " << c.tree.parent.size() << " nodes, min_count " << min_count);
		const auto backfill = BackfillNames(AssignNames(index, min_count), store, index.Tree());
		REQUIRE(backfill.labels == expected.labels);
		REQUIRE(backfill.inverted == expected.inverted);
	}
}

namespace {
// A store holding each lineage (domain first) as a row, so labels can name its names.
TaxonomyStore StoreOf(const std::vector<std::vector<std::string>> &lineages) {
	std::vector<TaxonomyRow> rows;
	for (const auto &lineage : lineages) {
		TaxonomyRow row {"R" + std::to_string(rows.size()), {}};
		std::copy(lineage.begin(), lineage.end(), row.ranks.begin());
		rows.push_back(std::move(row));
	}
	return TaxonomyStore(rows);
}

TreeIndex TreeOf(const std::vector<uint32_t> &parent) {
	return TreeIndex::Build(parent, std::vector<double>(parent.size(), 1.0));
}

// A label from text, "g__G1; s__S1", each name found at its rank in `store`.
NodeLabel LabelFrom(const TaxonomyStore &store, std::string_view text) {
	NodeLabel label;
	while (!text.empty()) {
		const size_t end = std::min(text.find("; "), text.size());
		const std::string_view part = text.substr(0, end);
		text.remove_prefix(std::min(end + 2, text.size()));
		const auto rank = std::find_if(kLineageRanks.begin(), kLineageRanks.end(),
		                               [&](const LineageRank &r) { return r.letter == part[0]; }) -
		                  kLineageRanks.begin();
		NameId name = 0;
		while (name < store.NameCount(rank) && store.NameText(rank, name) != part.substr(3)) {
			name++;
		}
		if (name == store.NameCount(rank) || part.substr(1, 2) != "__") {
			throw std::logic_error("LabelFrom: no name " + std::string(part));
		}
		label.push_back({static_cast<size_t>(rank), name});
	}
	return label;
}

// A lineage D1; P1; C1; O1; F1, then `deeper`.
std::vector<std::string> UnderF1(std::vector<std::string> deeper) {
	std::vector<std::string> out {"D1", "P1", "C1", "O1", "F1"};
	out.insert(out.end(), deeper.begin(), deeper.end());
	return out;
}

// The labels of n nodes, from text; unlisted nodes have empty labels.
std::vector<NodeLabel> LabelsOf(const TaxonomyStore &store, size_t n,
                                std::initializer_list<std::pair<uint32_t, std::string>> labels) {
	std::vector<NodeLabel> out(n);
	for (const auto &[node, text] : labels) {
		out[node] = LabelFrom(store, text);
	}
	return out;
}

// u(0) +-- v(1) +-- a(3) +-- t5, t6
//      |        +-- b(4) +-- t7, t8
//      +-- c(2) +-- t9, t10
// v is unnamed, so its labeled descendants are a and b. c's label, unlike a's and b's, keeps u's
// own promotion from taking their names.
const std::vector<uint32_t> kUvTree {kNoParent, 0, 0, 1, 1, 3, 3, 4, 4, 2, 2};
constexpr uint32_t kU = 0, kV = 1, kC = 2, kA = 3, kB = 4;
} // namespace

// A name moves up to the topmost node whose only labeled descendant holds it, not only to the
// labeled node's parent: tax2tree's behavior, kept. At R, H1 and Z disagree; at Q, H1 is the only
// labeled descendant, P being unnamed and W's subtree holding no names, so S lands on Q, two
// levels above H1, and spans W's tips, which S does not cover.
TEST_CASE("PromoteCommonNames hoists a lone descendant's name to the topmost node", "[taxonomy_decorator]") {
	// R(0) +-- Q(1) +-- P(3) +-- H1(5) --- t7
	//      |        |        +-- H2(6) --- t8
	//      |        +-- W(4) +-- t9, t10
	//      +-- Z(2) +-- t11, t12
	const auto store = StoreOf({{"D1", "P1", "C1", "O1", "F1", "G1", "S"}, {"D1", "P1", "C1", "O1", "F1", "Other"}});
	const auto tree = TreeOf({kNoParent, 0, 0, 1, 1, 3, 3, 5, 6, 4, 4, 2, 2});
	const auto labels = PromoteCommonNames(LabelsOf(store, 13, {{5, "s__S"}, {2, "g__Other"}}), store, tree);
	CHECK(LabelTexts(labels, store) == TextsOf(13, {{1, "s__S"}, {2, "g__Other"}}));
}

// A descendant's last name can be taken, as in tax2tree. a holds only F, which b also starts with, so F
// goes up to v, which spans them both, and a is left unnamed. Refusing to take a last name would
// leave F on a and on v, two nodes for one family.
TEST_CASE("PromoteCommonNames takes a descendant's last name", "[taxonomy_decorator]") {
	const auto store = StoreOf({{"D1", "P1", "C1", "O1", "F", "G"}, {"D1", "P1", "C1", "O1", "F9"}});
	const auto tree = TreeOf(kUvTree);
	const auto labels = PromoteCommonNames(
	    LabelsOf(store, 11, {{kU, "d__D1; p__P1; c__C1; o__O1"}, {kA, "f__F"}, {kB, "f__F; g__G"}, {kC, "f__F9"}}),
	    store, tree);
	CHECK(LabelTexts(labels, store) ==
	      TextsOf(11, {{kU, "d__D1; p__P1; c__C1; o__O1"}, {kV, "f__F"}, {kB, "g__G"}, {kC, "f__F9"}}));
}

// Promotion goes round by round while the first names agree, and stops after a round leaves a
// descendant with one name: here after C1, though a and b would agree on O1 next.
TEST_CASE("PromoteCommonNames stops after a round leaves a descendant with one name", "[taxonomy_decorator]") {
	const auto store = StoreOf({{"D1", "P1", "C1", "O1", "F2"}, {"D1", "P1", "C9"}});
	const auto tree = TreeOf(kUvTree);
	const auto labels = PromoteCommonNames(
	    LabelsOf(store, 11, {{kU, "d__D1; p__P1"}, {kA, "c__C1; o__O1"}, {kB, "c__C1; o__O1; f__F2"}, {kC, "c__C9"}}),
	    store, tree);
	CHECK(LabelTexts(labels, store) ==
	      TextsOf(11, {{kU, "d__D1; p__P1"}, {kV, "c__C1"}, {kA, "o__O1"}, {kB, "o__O1; f__F2"}, {kC, "c__C9"}}));
}

// Nodes are visited in preorder, so u, seeing a, b and c agree on D1 to O1, takes all four before
// v is visited. Visited bottom-up, v would take them first, and u would then take only D1 to C1
// from v and c, stopping when v is down to O1.
TEST_CASE("PromoteCommonNames visits ancestors before descendants", "[taxonomy_decorator]") {
	const auto store =
	    StoreOf({{"D1", "P1", "C1", "O1", "F1"}, {"D1", "P1", "C1", "O1", "F2"}, {"D1", "P1", "C1", "O1", "F3"}});
	const auto tree = TreeOf(kUvTree);
	const auto labels = PromoteCommonNames(LabelsOf(store, 11,
	                                                {{kA, "d__D1; p__P1; c__C1; o__O1; f__F1"},
	                                                 {kB, "d__D1; p__P1; c__C1; o__O1; f__F2"},
	                                                 {kC, "d__D1; p__P1; c__C1; o__O1; f__F3"}}),
	                                       store, tree);
	CHECK(LabelTexts(labels, store) ==
	      TextsOf(11, {{kU, "d__D1; p__P1; c__C1; o__O1"}, {kA, "f__F1"}, {kB, "f__F2"}, {kC, "f__F3"}}));
}

// Two names agree when one is the other's polyphyletic base, and the suffixed one goes up, as
// tax2tree merges them to keep GTDB's names. Nothing else agrees. tax2tree matches by substring, so there
// Baci and Bacillus_A agree, and so do Candidatus and Candidatus_Foo.
TEST_CASE("PromoteCommonNames merges a name with its polyphyletic variant, and nothing else", "[taxonomy_decorator]") {
	struct Pair {
		std::string a, b, promoted;
	};
	const Pair pairs[] = {{"Bacillus", "Bacillus_A", "Bacillus_A"},
	                      {"Bacillus_A", "Bacillus", "Bacillus_A"},
	                      {"Baci", "Bacillus_A", ""},
	                      {"Foo", "Foo_1", ""},
	                      {"Bacillus_A", "Bacillus_B", ""},
	                      {"Candidatus", "Candidatus_Foo", ""}};
	for (const auto &[a, b, promoted] : pairs) {
		INFO(a << " and " << b);
		const auto store = StoreOf({UnderF1({a, "S1"}), UnderF1({b, "S2"}), UnderF1({"Other"})});
		const auto tree = TreeOf(kUvTree);
		const std::string u = "d__D1; p__P1; c__C1; o__O1; f__F1";
		const auto labels = PromoteCommonNames(
		    LabelsOf(store, 11, {{kU, u}, {kA, "g__" + a + "; s__S1"}, {kB, "g__" + b + "; s__S2"}, {kC, "g__Other"}}),
		    store, tree);
		if (promoted.empty()) {
			CHECK(LabelTexts(labels, store) ==
			      TextsOf(11, {{kU, u}, {kA, "g__" + a + "; s__S1"}, {kB, "g__" + b + "; s__S2"}, {kC, "g__Other"}}));
		} else {
			CHECK(LabelTexts(labels, store) ==
			      TextsOf(11, {{kU, u}, {kV, "g__" + promoted}, {kA, "s__S1"}, {kB, "s__S2"}, {kC, "g__Other"}}));
		}
	}
}

// A name and its polyphyletic variant agree only at one rank: at two, they are two taxa. Here a
// holds the family Foo_A and b, inverted below u's order, the class Foo; merging them would put
// a family on v and take b's class away.
TEST_CASE("PromoteCommonNames never merges names at two ranks", "[taxonomy_decorator]") {
	const auto store =
	    StoreOf({{"D1", "P1", "C1", "O1", "Foo_A", "G1"}, {"D1", "P1", "Foo"}, {"D1", "P1", "C1", "O1", "Other"}});
	const auto tree = TreeOf(kUvTree);
	const std::string u = "d__D1; p__P1; c__C1; o__O1";
	const auto labels = PromoteCommonNames(
	    LabelsOf(store, 11, {{kU, u}, {kA, "f__Foo_A; g__G1"}, {kB, "c__Foo"}, {kC, "f__Other"}}), store, tree);
	CHECK(LabelTexts(labels, store) ==
	      TextsOf(11, {{kU, u}, {kA, "f__Foo_A; g__G1"}, {kB, "c__Foo"}, {kC, "f__Other"}}));
}

// Agreement counts distinct names: three never agree, while a name held twice is one name.
TEST_CASE("PromoteCommonNames never merges three distinct names", "[taxonomy_decorator]") {
	// u(0) +-- v(1) +-- a(3) +-- t6, t7
	//      |        +-- b(4) +-- t8, t9
	//      |        +-- e(5) +-- t10, t11
	//      +-- c(2) +-- t12, t13
	const auto store = StoreOf({UnderF1({"Bacillus", "S1"}), UnderF1({"Bacillus_A", "S2"}),
	                            UnderF1({"Bacillus_B", "S3"}), UnderF1({"Bacillus", "S4"}), UnderF1({"Other"})});
	const auto tree = TreeOf({kNoParent, 0, 0, 1, 1, 1, 3, 3, 4, 4, 5, 5, 2, 2});
	const std::string u = "d__D1; p__P1; c__C1; o__O1; f__F1";

	const auto three = PromoteCommonNames(LabelsOf(store, 14,
	                                               {{0, u},
	                                                {3, "g__Bacillus; s__S1"},
	                                                {4, "g__Bacillus_A; s__S2"},
	                                                {5, "g__Bacillus_B; s__S3"},
	                                                {2, "g__Other"}}),
	                                      store, tree);
	CHECK(LabelTexts(three, store) == TextsOf(14, {{0, u},
	                                               {3, "g__Bacillus; s__S1"},
	                                               {4, "g__Bacillus_A; s__S2"},
	                                               {5, "g__Bacillus_B; s__S3"},
	                                               {2, "g__Other"}}));

	const auto two = PromoteCommonNames(LabelsOf(store, 14,
	                                             {{0, u},
	                                              {3, "g__Bacillus; s__S1"},
	                                              {4, "g__Bacillus_A; s__S2"},
	                                              {5, "g__Bacillus; s__S4"},
	                                              {2, "g__Other"}}),
	                                    store, tree);
	CHECK(LabelTexts(two, store) ==
	      TextsOf(14, {{0, u}, {1, "g__Bacillus_A"}, {3, "s__S1"}, {4, "s__S2"}, {5, "s__S4"}, {2, "g__Other"}}));
}

// The suffix is the one tax2tree's promotion looks for: an underscore and capital letters, here at
// the end, after a non-empty base. tax2tree's pattern is anchored only at the start, so it finds a
// suffix in Candidatus_Foo and Foo_A1.
TEST_CASE("PolyphyleticBase strips a trailing capital-letter suffix", "[taxonomy_decorator]") {
	const std::pair<std::string_view, std::string_view> cases[] = {
	    {"Bacillus_A", "Bacillus"},
	    {"Bacillus_AB", "Bacillus"},
	    {"g__Bacillus_A", "g__Bacillus"},
	    {"A_B_C", "A_B"},
	    {"Bacillus_A subtilis_B", "Bacillus_A subtilis"},
	    {"Bacillus", "Bacillus"},
	    {"Foo_1", "Foo_1"},
	    {"Foo_A1", "Foo_A1"},
	    {"Foo_a", "Foo_a"},
	    {"Foo_", "Foo_"},
	    {"_A", "_A"},
	    {"Candidatus_Foo", "Candidatus_Foo"},
	    {"[Eubacterium]", "[Eubacterium]"},
	};
	for (const auto &[name, base] : cases) {
		INFO(name);
		CHECK(PolyphyleticBase(name) == base);
	}
}

// A promoted name says the node is that taxon, so it must fit below what the node already says.
// v is named at the species Z, and n, its one labeled descendant, at the genus Y, inverted below
// it: Y stays on n, the clade it was placed on. tax2tree moves it to v, writing "s__Z; g__Y", a
// genus over a species of another genus.
TEST_CASE("PromoteCommonNames never moves a name onto a node named deeper", "[taxonomy_decorator]") {
	// u(0) +-- v(1) +-- n(3) +-- t4, t5
	//      |        +-- t6
	//      +-- c(2) +-- t7, t8
	const auto store = StoreOf({UnderF1({"G1", "Z"}), UnderF1({"G1", "Other"}), UnderF1({"Y"})});
	const auto tree = TreeOf({kNoParent, 0, 0, 1, 3, 3, 1, 2, 2});
	const std::string u = "d__D1; p__P1; c__C1; o__O1; f__F1; g__G1";
	const auto labels =
	    PromoteCommonNames(LabelsOf(store, 9, {{0, u}, {1, "s__Z"}, {3, "g__Y"}, {2, "s__Other"}}), store, tree);
	CHECK(LabelTexts(labels, store) == TextsOf(9, {{0, u}, {1, "s__Z"}, {3, "g__Y"}, {2, "s__Other"}}));
}

// A label holds one name per rank. v holds F1 and G, and n, inverted below it, holds F2 alone:
// tax2tree moves F2 to v, writing "f__F1; g__G; f__F2", two families on one node. Here F2 stays on n.
TEST_CASE("PromoteCommonNames never gives a node a second name at one rank", "[taxonomy_decorator]") {
	// u(0) +-- v(1) +-- n(3) +-- t4, t5
	//      |        +-- t6
	//      +-- c(2) +-- t7, t8
	const auto store =
	    StoreOf({{"D1", "P1", "C1", "O1", "F1", "G"}, {"D1", "P1", "C1", "O1", "F2"}, {"D1", "P1", "C1", "O1", "F3"}});
	const auto tree = TreeOf({kNoParent, 0, 0, 1, 3, 3, 1, 2, 2});
	const std::string u = "d__D1; p__P1; c__C1; o__O1";
	const auto labels =
	    PromoteCommonNames(LabelsOf(store, 9, {{0, u}, {1, "f__F1; g__G"}, {3, "f__F2"}, {2, "f__F3"}}), store, tree);
	CHECK(LabelTexts(labels, store) == TextsOf(9, {{0, u}, {1, "f__F1; g__G"}, {3, "f__F2"}, {2, "f__F3"}}));
}

// A promoted name must be deeper than every rank the node's path already holds, not only its own
// label's: an empty-labeled node already says what its ancestors say. So an inverted name stays
// on the node it was placed on, where the inversion record finds it, and an inverted node takes
// no name at a rank its path holds above it.
TEST_CASE("PromoteCommonNames leaves an inverted name where it was placed", "[taxonomy_decorator]") {
	SECTION("an empty-labeled node takes no inverted name") {
		// u [d..o] +-- v [] +-- a [c__C2] +-- t5, t6
		//          |        +-- b []      +-- t7, t8
		//          +-- c [f__F9]          +-- t9, t10
		const auto store = StoreOf({{"D1", "P1", "C1", "O1", "F9"}, {"D1", "P1", "C2"}});
		const auto tree = TreeOf(kUvTree);
		const std::string u = "d__D1; p__P1; c__C1; o__O1";
		const auto labels =
		    PromoteCommonNames(LabelsOf(store, 11, {{kU, u}, {kA, "c__C2"}, {kC, "f__F9"}}), store, tree);
		CHECK(LabelTexts(labels, store) == TextsOf(11, {{kU, u}, {kA, "c__C2"}, {kC, "f__F9"}}));
	}
	SECTION("an inverted node takes no name at a rank its path holds") {
		// u(0) [d..o] +-- d(1) [c__C2] +-- e(2) [o__Ox; f__Fx] +-- t3, t4
		//             +-- c(5) [f__F9] +-- t6, t7
		const auto store = StoreOf({{"D1", "P1", "C1", "O1", "F9"}, {"D1", "P1", "C2", "Ox", "Fx"}});
		const auto tree = TreeOf({kNoParent, 0, 1, 2, 2, 0, 5, 5});
		const std::string u = "d__D1; p__P1; c__C1; o__O1";
		const auto labels = PromoteCommonNames(
		    LabelsOf(store, 8, {{0, u}, {1, "c__C2"}, {2, "o__Ox; f__Fx"}, {5, "f__F9"}}), store, tree);
		CHECK(LabelTexts(labels, store) == TextsOf(8, {{0, u}, {1, "c__C2"}, {2, "o__Ox; f__Fx"}, {5, "f__F9"}}));
	}
}

namespace {
std::vector<LabelPart> Parts(std::initializer_list<std::pair<size_t, std::string>> parts) {
	std::vector<LabelPart> out;
	for (const auto &[rank, name] : parts) {
		out.push_back({rank, name});
	}
	return out;
}

// Each tip's lineage as "tip: name|name|...", eight fields, empty where no label gives the rank.
std::vector<std::string> LineageTexts(const EmittedLineages &emitted,
                                      const std::vector<std::vector<LabelPart>> &labels) {
	std::vector<std::string> out;
	for (size_t i = 0; i < emitted.tips.size() && i < emitted.holders.size(); i++) {
		std::string line = std::to_string(emitted.tips[i]) + ": ";
		for (size_t r = 0; r < kLineageRanks.size(); r++) {
			if (r > 0) {
				line += "|";
			}
			line += HolderName(labels, emitted.holders[i][r], r);
		}
		out.push_back(line);
	}
	return out;
}

std::string ParseError(std::string_view label) {
	try {
		ParseNodeLabel(label);
	} catch (const std::invalid_argument &e) {
		return e.what();
	}
	return "";
}
} // namespace

// Labels come from FormatLabel, from other tools, and from hand editing, so spacing after ';'
// varies and ranks need not come in order: each component is placed by its letter. k is the
// domain's letter in trees labeled with kingdoms.
TEST_CASE("ParseNodeLabel reads components by their rank letter", "[taxonomy_decorator]") {
	CHECK(ParseNodeLabel("f__A; g__B") == Parts({{kFamily, "A"}, {kGenus, "B"}}));
	CHECK(ParseNodeLabel("f__A;g__B") == Parts({{kFamily, "A"}, {kGenus, "B"}}));
	CHECK(ParseNodeLabel(" f__A ;\tg__B ") == Parts({{kFamily, "A"}, {kGenus, "B"}}));
	CHECK(ParseNodeLabel("f__A\f;\vg__B") == Parts({{kFamily, "A"}, {kGenus, "B"}}));
	CHECK(ParseNodeLabel("g__B; f__A") == Parts({{kFamily, "A"}, {kGenus, "B"}}));
	CHECK(ParseNodeLabel("k__Bacteria") == Parts({{kDomain, "Bacteria"}}));
	CHECK(ParseNodeLabel("g__Escherichia; s__Escherichia coli; t__K-12") ==
	      Parts({{kGenus, "Escherichia"}, {kSpecies, "Escherichia coli"}, {kStrain, "K-12"}}));
	CHECK(ParseNodeLabel("g__[Eubacterium]") == Parts({{kGenus, "[Eubacterium]"}}));
	CHECK(ParseNodeLabel("").empty());
}

// An empty name says the label has no name at that rank, as in Greengenes' "g__; s__", so long as
// the label names no deeper rank; above a named rank it is a gap, and the label is malformed. The
// error names the shallowest empty component, the one the deeper names leave a gap below.
TEST_CASE("ParseNodeLabel reads an empty name below every named rank as no name", "[taxonomy_decorator]") {
	CHECK(ParseNodeLabel("f__").empty());
	CHECK(ParseNodeLabel("g__; s__").empty());
	CHECK(ParseNodeLabel("p__B; c__") == Parts({{kPhylum, "B"}}));
	CHECK(ParseNodeLabel("s__; g__G") == Parts({{kGenus, "G"}}));
	for (const std::string_view gap : {"p__B; c__; o__D", "o__D; c__", "c__; f__; g__G", "s__; c__; o__D"}) {
		INFO(gap);
		const std::string error = ParseError(gap);
		CHECK_THAT(error, ContainsSubstring(std::string(gap)));
		CHECK_THAT(error, ContainsSubstring("\"c__\"") && ContainsSubstring("gap"));
	}
}

// The ranks a label names are contiguous: decoration fills every rank between a node's name and its
// named ancestor's from the taxonomy, so a label that leaves one out is malformed, like an empty
// name above a named rank. The error names the component below the rank left out.
TEST_CASE("ParseNodeLabel rejects a label that leaves out a rank between two it names", "[taxonomy_decorator]") {
	const std::pair<std::string_view, std::string_view> cases[] = {
	    {"p__B; o__D", "o__D"},
	    {"o__D; p__B", "o__D"},
	    {"d__Bacteria; s__Escherichia coli; t__K-12", "s__Escherichia coli"},
	};
	for (const auto &[label, component] : cases) {
		INFO(label);
		const std::string error = ParseError(label);
		CHECK_THAT(error, ContainsSubstring(std::string(label)));
		CHECK_THAT(error, ContainsSubstring("\"" + std::string(component) + "\"") && ContainsSubstring("gap"));
	}
}

// A label that cannot be read is an error naming it and the component, never a guess: tax2tree
// indexes a rank table with the first character and fails with a KeyError, or reads "x__" as a
// name. A name may not start with whitespace, as no taxonomy name does, or with a rank prefix of
// its own; whitespace at either end of a component is trimmed.
TEST_CASE("ParseNodeLabel rejects a malformed label, naming it", "[taxonomy_decorator]") {
	const std::pair<std::string_view, std::string_view> cases[] = {
	    {"x__A", "x__A"},       {"D__A", "D__A"},     {"fA", "fA"},           {"f_A", "f_A"},
	    {"f__ A", "f__ A"},     {"f__A;", ""},        {"f__A; f__B", "f__B"}, {"c__; c__C", "c__C"},
	    {"d__A; k__B", "k__B"}, {"f__\vA", "f__\vA"}, {"c__c__A", "c__c__A"},
	};
	for (const auto &[label, component] : cases) {
		INFO(label);
		const std::string error = ParseError(label);
		CHECK_THAT(error, ContainsSubstring(std::string(label)));
		CHECK_THAT(error, ContainsSubstring("\"" + std::string(component) + "\""));
	}
}

// The text form of a label is the one tax2tree writes: names in rank order, each after its
// letter, joined by "; ".
TEST_CASE("FormatLabel writes names in rank order with their letters", "[taxonomy_decorator]") {
	const auto store = StoreOf({{"Bacteria", "Pseudomonadota", "Gammaproteobacteria", "Enterobacterales",
	                             "Enterobacteriaceae", "Escherichia", "Escherichia coli"}});
	CHECK(FormatLabel(LabelFrom(store, "d__Bacteria; p__Pseudomonadota"), store) == "d__Bacteria; p__Pseudomonadota");
	CHECK(FormatLabel(LabelFrom(store, "g__Escherichia; s__Escherichia coli"), store) ==
	      "g__Escherichia; s__Escherichia coli");
	CHECK(FormatLabel(LabelFrom(store, "d__Bacteria"), store) == "d__Bacteria");
	CHECK(FormatLabel({}, store).empty());
}

// A tip's lineage has every rank, NULL where no label on its path gives one. The root's label
// counts, a node with no label adds nothing, and a tip's own label is ignored: a tip's name is
// its id.
TEST_CASE("EmitLineages gives each tip eight ranks from the labels on its path", "[taxonomy_decorator]") {
	// root(0) [d__D1; p__P1] +-- v(1) [c__C1] +-- t3, t4
	//                        +-- w(2)         +-- t5
	//                        +-- t6 [g__X]
	const auto tree = TreeOf({kNoParent, 0, 0, 1, 1, 2, 0});
	std::vector<std::vector<LabelPart>> labels(7);
	labels[0] = Parts({{kDomain, "D1"}, {kPhylum, "P1"}});
	labels[1] = Parts({{kClass, "C1"}});
	labels[6] = Parts({{kGenus, "X"}});
	const auto emitted = EmitLineages(tree, labels);
	CHECK(emitted.tips == std::vector<uint32_t> {3, 4, 5, 6});
	CHECK(LineageTexts(emitted, labels) ==
	      std::vector<std::string> {"3: D1|P1|C1|||||", "4: D1|P1|C1|||||", "5: D1|P1||||||", "6: D1|P1||||||"});
	CHECK(emitted.repeats.empty());
}

// A rank held twice on a path is recorded, even with the same name both times, and the shallower
// name is the one emitted, as tax2tree does. The caller decides whether a repeat is an error.
TEST_CASE("EmitLineages records a rank held twice on a path", "[taxonomy_decorator]") {
	// root(0) [d__D1; p__P1] +-- v(1) [p__P2; c__C1] +-- x(2) [c__C1] +-- t3, t4
	//                        |                       +-- t5
	//                        +-- t6
	const auto tree = TreeOf({kNoParent, 0, 1, 2, 2, 1, 0});
	std::vector<std::vector<LabelPart>> labels(7);
	labels[0] = Parts({{kDomain, "D1"}, {kPhylum, "P1"}});
	labels[1] = Parts({{kPhylum, "P2"}, {kClass, "C1"}});
	labels[2] = Parts({{kClass, "C1"}});
	const auto emitted = EmitLineages(tree, labels);
	CHECK(LineageTexts(emitted, labels) ==
	      std::vector<std::string> {"3: D1|P1|C1|||||", "4: D1|P1|C1|||||", "5: D1|P1|C1|||||", "6: D1|P1||||||"});
	CHECK(emitted.repeats == std::vector<RankRepeat> {{1, kPhylum, 0}, {2, kClass, 1}});
}

// Repeats come by node, then rank, whatever order the walk meets them in. Here node 2 is above
// node 1, so the walk meets 2's repeat first.
TEST_CASE("EmitLineages lists repeats by node, then rank", "[taxonomy_decorator]") {
	// root(0) [d__D1; p__P1] +-- w(2) [p__P2; c__C2] +-- v(1) [p__P3; c__C3] +-- t3, t4
	//                                                +-- t5
	const auto tree = TreeOf({kNoParent, 2, 0, 1, 1, 2});
	std::vector<std::vector<LabelPart>> labels(6);
	labels[0] = Parts({{kDomain, "D1"}, {kPhylum, "P1"}});
	labels[2] = Parts({{kPhylum, "P2"}, {kClass, "C2"}});
	labels[1] = Parts({{kPhylum, "P3"}, {kClass, "C3"}});
	const auto emitted = EmitLineages(tree, labels);
	CHECK(LineageTexts(emitted, labels) ==
	      std::vector<std::string> {"3: D1|P1|C2|||||", "4: D1|P1|C2|||||", "5: D1|P1|C2|||||"});
	CHECK(emitted.repeats == std::vector<RankRepeat> {{1, kPhylum, 0}, {1, kClass, 2}, {2, kPhylum, 0}});
}

TEST_CASE("EmitLineages on a one-node tree, and with labels of the wrong size", "[taxonomy_decorator]") {
	const auto tree = TreeOf({kNoParent});
	const std::vector<std::vector<LabelPart>> labels(1);
	const auto emitted = EmitLineages(tree, labels);
	CHECK(LineageTexts(emitted, labels) == std::vector<std::string> {"0: |||||||"});
	CHECK_THROWS_AS(EmitLineages(TreeOf({kNoParent, 0, 0}), labels), std::invalid_argument);
}

// ParseNodeLabel never gives a label that holds a rank twice, so one breaks EmitLineages'
// precondition: an error naming the node and the rank, never a repeat whose ancestor is the node.
TEST_CASE("EmitLineages rejects a label that holds a rank twice", "[taxonomy_decorator]") {
	const auto tree = TreeOf({kNoParent, 0, 1, 1});
	std::vector<std::vector<LabelPart>> labels(4);
	labels[1] = Parts({{kFamily, "A"}, {kFamily, "B"}});
	CHECK_THROWS_AS(EmitLineages(tree, labels), std::invalid_argument);
	CHECK_THROWS_WITH(EmitLineages(tree, labels), ContainsSubstring("node 1") && ContainsSubstring("family"));
	// With the rank already held above, as a repeat would be.
	labels[0] = Parts({{kFamily, "R"}});
	CHECK_THROWS_WITH(EmitLineages(tree, labels), ContainsSubstring("node 1") && ContainsSubstring("family"));
	// A rank past the last is no rank at all.
	labels[1] = Parts({{kLineageRanks.size(), "A"}});
	CHECK_THROWS_WITH(EmitLineages(tree, labels), ContainsSubstring("node 1") && ContainsSubstring("rank 8"));
}

// The stages in order. F's tips are split between v's two clades and w, so F is eligible at w
// alone, and a and b are named at the genus. Backfill gives a and b their family from the
// hierarchy, and promotion then moves F up to v, which spans them both.
TEST_CASE("Decorate places, backfills and promotes", "[taxonomy_decorator]") {
	// root(0) +-- v(1) +-- a(4) +-- t6, t7  D1; P1; C1; O1; F; G1
	//         |        +-- b(5) +-- t8, t9  D1; P1; C1; O1; F; G2
	//         +-- w(2) +-- t10 .. t14       D1; P1; C1; O1; F
	//         +-- c(3) +-- t15, t16         D1; P1; C1; O1; F9
	Lineages lineages;
	Add(lineages, {"D1", "P1", "C1", "O1", "F", "G1"}, {6, 7});
	Add(lineages, {"D1", "P1", "C1", "O1", "F", "G2"}, {8, 9});
	Add(lineages, {"D1", "P1", "C1", "O1", "F"}, {10, 11, 12, 13, 14});
	Add(lineages, {"D1", "P1", "C1", "O1", "F9"}, {15, 16});
	auto c = MakeCase({kNoParent, 0, 0, 0, 1, 1, 4, 4, 5, 5, 2, 2, 2, 2, 2, 3, 3}, lineages);
	const auto result = Decorate(c.index, c.store, 2);
	CHECK(result.assignment.Scores() == AssignNames(c.index, 2).Scores());
	CHECK(LabelTexts(result.labels, c.store) ==
	      TextsOf(
	          17,
	          {{0, "d__D1; p__P1; c__C1; o__O1"}, {1, "f__F"}, {2, "f__F"}, {3, "f__F9"}, {4, "g__G1"}, {5, "g__G2"}}));
	CHECK(result.inverted.empty());
}

// tax2tree's validate_all_paths, on our own output: on random trees, every label survives
// FormatLabel and ParseNodeLabel, and wherever no inverted-rank backfill was recorded, each tip's
// path holds every rank at most once, in order from the root, so EmitLineages records no repeat.
// An inversion is the one source of a repeat: an inverted node keeps its name, every repeat is at
// an inverted node or below one, and the paths to the other tips stay in rank order. Rank safety
// keeps two names at one rank off a path: a node holds a name only where at least half of that
// name's tips are, so every ancestor has as many, and an ancestor holding another name at that rank
// would have two eligible names there. The trials must include both kinds, and some promotion.
TEST_CASE("Decorate's labels keep every path in rank order unless a backfill was inverted", "[taxonomy_decorator]") {
	std::mt19937 rng(20261121);
	size_t with_inversion = 0, without_inversion = 0, promoted = 0;
	for (int k = 0; k < 300; k++) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 40)(rng);
		const auto c = random_taxonomy_case(rng, n);
		TaxonomyStore store(c.rows);
		const auto index = build_index(c, store);
		for (uint32_t min_count = 1; min_count <= 3; min_count++) {
			INFO("trial " << k << ", " << n << " nodes, min_count " << min_count);
			const auto result = Decorate(index, store, min_count);
			REQUIRE(result.labels.size() == n);
			const auto backfill = BackfillNames(result.assignment, store, index.Tree());
			std::vector<std::vector<LabelPart>> parsed(n);
			for (uint32_t v = 0; v < n; v++) {
				INFO("node " << v);
				promoted += result.labels[v].size() > backfill.labels[v].size();
				parsed[v] = ParseNodeLabel(FormatLabel(result.labels[v], store));
				REQUIRE(parsed[v].size() == result.labels[v].size());
				for (size_t i = 0; i < parsed[v].size(); i++) {
					REQUIRE(parsed[v][i].rank == result.labels[v][i].rank);
					REQUIRE(parsed[v][i].name == store.NameText(result.labels[v][i].rank, result.labels[v][i].name));
				}
			}
			const auto emitted = EmitLineages(index.Tree(), parsed);
			REQUIRE(emitted.tips.size() == emitted.holders.size());
			// Whether v is an inverted node or below one.
			const auto below_inversion = [&](uint32_t v) {
				for (; v != kNoParent; v = c.tree.parent[v]) {
					for (const auto &inverted : result.inverted) {
						if (inverted.node == v) {
							return true;
						}
					}
				}
				return false;
			};
			for (const auto &inverted : result.inverted) {
				INFO("inverted node " << inverted.node);
				const auto &label = result.labels[inverted.node];
				CHECK(std::find(label.begin(), label.end(), RankName {inverted.rank, inverted.name}) != label.end());
			}
			for (const auto &repeat : emitted.repeats) {
				INFO("repeat at node " << repeat.node);
				CHECK(below_inversion(repeat.node));
			}
			if (result.inverted.empty()) {
				without_inversion++;
				CHECK(emitted.repeats.empty());
			} else {
				with_inversion++;
			}
			for (const uint32_t tip : emitted.tips) {
				if (below_inversion(tip)) {
					continue;
				}
				// The ranks of the labels on the tip's path, from the root down.
				std::vector<size_t> ranks;
				for (uint32_t v = c.tree.parent[tip]; v != kNoParent; v = c.tree.parent[v]) {
					for (auto it = result.labels[v].rbegin(); it != result.labels[v].rend(); ++it) {
						ranks.push_back(it->rank);
					}
				}
				std::reverse(ranks.begin(), ranks.end());
				INFO("tip " << tip);
				CHECK(std::adjacent_find(ranks.begin(), ranks.end(), std::greater_equal<size_t>()) == ranks.end());
			}
		}
	}
	CHECK(with_inversion > 0);
	CHECK(without_inversion > 0);
	CHECK(promoted > 0);
}
