#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "NearestTip.hpp"
#include "NewickTree.hpp"
#include "TreeIndex.hpp"
#include "tree_test_util.hpp"

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using miint::NearestTipOptions;
using miint::NearestTipResult;
using miint::NearestTipScope;
using miint::NearestTipTargets;
using miint::TreeIndex;
using namespace miint::test;

namespace {

constexpr uint32_t NO_TARGET = NearestTipResult::NO_TARGET;

// A Newick tree searched by name. Target ranks follow name order, as the SQL function
// ranks them.
struct NamedTree {
	miint::NewickTree tree;
	TreeIndex index;

	explicit NamedTree(const char *newick) : tree(miint::NewickTree::parse(newick)), index(index_of(tree)) {
	}

	uint32_t node(const std::string &name) const {
		auto found = tree.find_node_by_name(name);
		REQUIRE(found.has_value());
		return found.value();
	}

	std::string name(uint32_t node) const {
		return node == NO_TARGET ? std::string("<none>") : tree.name(node);
	}

	NearestTipTargets targets(std::vector<std::string> names) const {
		std::sort(names.begin(), names.end());
		std::vector<uint32_t> nodes;
		std::vector<uint32_t> ranks;
		for (uint32_t r = 0; r < names.size(); ++r) {
			nodes.push_back(node(names[r]));
			ranks.push_back(r);
		}
		return NearestTipTargets(index, nodes, ranks);
	}
};

struct Found {
	std::string target;
	double distance;
	uint64_t n_tied;
	std::string lca;
};

Found search(const NamedTree &t, const NearestTipTargets &targets, const std::string &query, bool exclude_self = true,
             NearestTipScope scope = NearestTipScope::Path) {
	NearestTipOptions options;
	options.exclude_self = exclude_self;
	options.scope = scope;
	auto r = targets.Search(t.node(query), options);
	return {t.name(r.target), r.distance, r.n_tied, t.name(r.lca)};
}

// The brute-force answer: every eligible target scored with the same distance formula,
// from the naive model's LCA and root distances. For the clade scope, only the targets
// in the lowest ancestor-or-self of the query that holds any eligible target compete.
NearestTipResult oracle(const TreeArrays &t, const std::vector<uint32_t> &nodes, const std::vector<uint32_t> &ranks,
                        uint32_t query, bool exclude_self, NearestTipScope scope = NearestTipScope::Path) {
	auto eligible = [&](uint32_t target) {
		return !(exclude_self && target == query);
	};
	uint32_t clade = kNoParent;
	if (scope == NearestTipScope::Clade) {
		for (uint32_t a = query; a != kNoParent && clade == kNoParent; a = t.parent[a]) {
			for (uint32_t target : nodes) {
				if (eligible(target) && naive_is_ancestor_or_self(t, a, target)) {
					clade = a;
					break;
				}
			}
		}
	}

	NearestTipResult best;
	uint32_t best_rank = UINT32_MAX;
	const double dq = naive_root_distance(t, query);
	for (size_t i = 0; i < nodes.size(); ++i) {
		const uint32_t target = nodes[i];
		if (!eligible(target)) {
			continue;
		}
		if (scope == NearestTipScope::Clade && !naive_is_ancestor_or_self(t, clade, target)) {
			continue;
		}
		const uint32_t l = naive_lca(t, query, target);
		const double dl = naive_root_distance(t, l);
		const double d = (dq - dl) + (naive_root_distance(t, target) - dl);
		if (best.target == NO_TARGET || d < best.distance) {
			best = NearestTipResult {target, d, 1, l};
			best_rank = ranks[i];
		} else if (d == best.distance) {
			++best.n_tied;
			if (ranks[i] < best_rank) {
				best.target = target;
				best.lca = l;
				best_rank = ranks[i];
			}
		}
	}
	return best;
}

} // namespace

TEST_CASE("NearestTip: a worked example, every edge on the path counted", "[NearestTip]") {
	NamedTree t("(((A:0.1,B:0.2)C:0.3,(D:0.4)E:0.5)x:0.1,(F:0.6,G:0.7)H:0.8)root;");
	auto targets = t.targets({"A", "D"});

	auto b = search(t, targets, "B");
	REQUIRE(b.target == "A");
	REQUIRE(b.distance == Approx(0.3).margin(1e-12));
	REQUIRE(b.n_tied == 1);
	REQUIRE(b.lca == "C");

	auto f = search(t, targets, "F");
	REQUIRE(f.target == "A");
	REQUIRE(f.distance == Approx(1.9).margin(1e-12));
	REQUIRE(f.lca == "root");

	auto g = search(t, targets, "G");
	REQUIRE(g.target == "A");
	REQUIRE(g.distance == Approx(2.0).margin(1e-12));
	REQUIRE(g.lca == "root");
}

TEST_CASE("NearestTip: a tie is counted and goes to the smallest name", "[NearestTip]") {
	NamedTree t("(Q:0.1,(T1:0.2,T2:0.2)n:0.1)r;");
	auto q = search(t, t.targets({"T1", "T2"}), "Q");
	REQUIRE(q.target == "T1");
	REQUIRE(q.distance == Approx(0.4).margin(1e-12));
	REQUIRE(q.n_tied == 2);
	REQUIRE(q.lca == "r");
}

TEST_CASE("NearestTip: a tie across levels goes to the smallest name, not the first found", "[NearestTip]") {
	// Z (under a) and B (under the root) are both at distance 2 from Q. The search reaches
	// Z first; B wins because its name is smaller, not because it was found first.
	NamedTree t("((Q:1,Z:1)a:1,B:0)r;");
	auto q = search(t, t.targets({"B", "Z"}), "Q");
	REQUIRE(q.target == "B");
	REQUIRE(q.distance == 2.0);
	REQUIRE(q.n_tied == 2);
	REQUIRE(q.lca == "r");
}

TEST_CASE("NearestTip: the nearest target counts the edge down to it", "[NearestTip]") {
	// From Q, the clade under `far` is entered by a long edge and the one under `near` by a
	// short one; the nearest tip is under `near`, though `far`'s own tip edge is shorter.
	// Comparing children without the edge down to each would pick `far`.
	NamedTree t("(Q:1,(F:0.1)far:5,(N:2)near:1)r;");
	auto q = search(t, t.targets({"F", "N"}), "Q");
	REQUIRE(q.target == "N");
	REQUIRE(q.distance == 4.0);
}

TEST_CASE("NearestTip: a query that is a target", "[NearestTip]") {
	NamedTree t("((T1:0.1,T2:0.2)a:0.1,T3:0.5)r;");
	auto targets = t.targets({"T1", "T2", "T3"});

	auto excluded = search(t, targets, "T1", true);
	REQUIRE(excluded.target == "T2");
	REQUIRE(excluded.distance == Approx(0.3).margin(1e-12));
	REQUIRE(excluded.n_tied == 1);
	REQUIRE(excluded.lca == "a");

	auto included = search(t, targets, "T1", false);
	REQUIRE(included.target == "T1");
	REQUIRE(included.distance == 0.0);
	REQUIRE(included.n_tied == 1);
	REQUIRE(included.lca == "T1");
}

TEST_CASE("NearestTip: no eligible target", "[NearestTip]") {
	NamedTree t("((A:1,B:1)c:1,D:1)r;");
	auto none = search(t, t.targets({}), "A");
	REQUIRE(none.target == "<none>");
	REQUIRE(none.n_tied == 0);
	REQUIRE(none.lca == "<none>");

	// The only target is the query itself, excluded.
	auto self_only = search(t, t.targets({"A"}), "A", true);
	REQUIRE(self_only.target == "<none>");
	REQUIRE(self_only.n_tied == 0);
}

TEST_CASE("NearestTip: answers do not depend on the order of children", "[NearestTip]") {
	// The same tree written with every node's children reversed. Re-serializing a tree
	// must not change any answer, ties included.
	NamedTree forward("(((A:1,B:1)c:0,(Q:1,Z:1)a:1)m:0,(Y:2,(X:0,W:0)n:2)k:0)r;");
	NamedTree reversed("(((W:0,X:0)n:2,Y:2)k:0,((Z:1,Q:1)a:1,(B:1,A:1)c:0)m:0)r;");
	const std::vector<std::string> target_names = {"A", "B", "W", "X", "Y", "Z"};
	auto forward_targets = forward.targets(target_names);
	auto reversed_targets = reversed.targets(target_names);

	for (const char *query : {"A", "B", "Q", "W", "X", "Y", "Z"}) {
		for (bool exclude_self : {true, false}) {
			CAPTURE(query, exclude_self);
			auto f = search(forward, forward_targets, query, exclude_self);
			auto r = search(reversed, reversed_targets, query, exclude_self);
			REQUIRE(f.target == r.target);
			REQUIRE(f.distance == r.distance);
			REQUIRE(f.n_tied == r.n_tied);
			REQUIRE(f.lca == r.lca);
		}
	}
}

TEST_CASE("NearestTip: a deep caterpillar", "[NearestTip]") {
	// Tip k (k < levels) hangs off internal node k at root distance k + 1. From tip k, a
	// target tip m < k is reached through internal node m, at distance (k + 1 - m) + 1.
	const uint32_t levels = 200000;
	auto t = caterpillar(levels);
	auto index = TreeIndex::Build(t.parent, t.length);
	auto tip = [&](uint32_t k) {
		return levels + k;
	};
	NearestTipTargets targets(index, {tip(0), tip(levels / 2)}, {0, 1});

	const uint32_t k = levels - 1;
	auto r = targets.Search(tip(k));
	REQUIRE(r.target == tip(levels / 2));
	REQUIRE(r.distance == static_cast<double>(k - levels / 2 + 2));
	REQUIRE(r.n_tied == 1);
	REQUIRE(r.lca == levels / 2);
}

TEST_CASE("NearestTip: path and clade scope disagree", "[NearestTip]") {
	// T2 is nearer by path, but T1 is in the smallest clade around Q that holds a target.
	NamedTree t("((Q:0.1,T1:5.0)n1:0.1,T2:0.1)r;");
	auto targets = t.targets({"T1", "T2"});

	auto path = search(t, targets, "Q", true, NearestTipScope::Path);
	REQUIRE(path.target == "T2");
	REQUIRE(path.distance == Approx(0.3).margin(1e-12));
	REQUIRE(path.lca == "r");

	auto clade = search(t, targets, "Q", true, NearestTipScope::Clade);
	REQUIRE(clade.target == "T1");
	REQUIRE(clade.distance == Approx(5.1).margin(1e-12));
	REQUIRE(clade.n_tied == 1);
	REQUIRE(clade.lca == "n1");
}

TEST_CASE("NearestTip: the clade scope counts ties inside the clade", "[NearestTip]") {
	// The smallest clade around Q holding a target is the root's; T1 and T2 tie in it.
	NamedTree t("(Q:1,(T1:1,T2:1)n:1)r;");
	auto q = search(t, t.targets({"T1", "T2"}), "Q", true, NearestTipScope::Clade);
	REQUIRE(q.target == "T1");
	REQUIRE(q.distance == 3.0);
	REQUIRE(q.n_tied == 2);
	REQUIRE(q.lca == "r");
}

TEST_CASE("NearestTip: a query that may be its own target, by scope", "[NearestTip]") {
	// Q's sibling T is at distance 0. Under 'path' it ties with Q itself; under 'clade' the
	// smallest clade holding an eligible target is Q's own (a tip), so T is not in it.
	NamedTree t("((Q:0,T:0)a:1,U:1)r;");
	auto targets = t.targets({"Q", "T", "U"});

	auto path = search(t, targets, "Q", false, NearestTipScope::Path);
	REQUIRE(path.target == "Q");
	REQUIRE(path.distance == 0.0);
	REQUIRE(path.n_tied == 2);

	auto clade = search(t, targets, "Q", false, NearestTipScope::Clade);
	REQUIRE(clade.target == "Q");
	REQUIRE(clade.distance == 0.0);
	REQUIRE(clade.n_tied == 1);
	REQUIRE(clade.lca == "Q");

	// Excluding itself, Q's smallest clade with a target is a's, holding T.
	auto excluded = search(t, targets, "Q", true, NearestTipScope::Clade);
	REQUIRE(excluded.target == "T");
	REQUIRE(excluded.lca == "a");
}

TEST_CASE("NearestTip: random trees and groups agree with a brute-force search", "[NearestTip]") {
	// Several target sets (groups) over one tree index, as the SQL function builds them:
	// a name rank per tip shared by every group, and groups that overlap. Every tip is a
	// query against every group, in both scopes, with and without itself.
	std::mt19937 rng(20261106);
	const double target_fractions[] = {0.0, 0.1, 0.5, 1.0};
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	for (int trial = 0; trial < 150; ++trial) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 120)(rng);
		const bool dyadic = trial % 3 != 0;
		auto t = random_tree(rng, n, dyadic);
		auto index = TreeIndex::Build(t.parent, t.length);

		std::vector<uint32_t> tips;
		for (uint32_t v = 0; v < n; ++v) {
			if (naive_is_tip(t, v)) {
				tips.push_back(v);
			}
		}
		std::vector<uint32_t> tip_rank(n, 0);
		{
			std::vector<uint32_t> order(tips.size());
			std::iota(order.begin(), order.end(), 0u);
			std::shuffle(order.begin(), order.end(), rng);
			for (size_t i = 0; i < tips.size(); ++i) {
				tip_rank[tips[i]] = order[i];
			}
		}

		const int n_groups = std::uniform_int_distribution<int>(1, 4)(rng);
		std::vector<std::vector<uint32_t>> group_nodes(n_groups);
		std::vector<std::vector<uint32_t>> group_ranks(n_groups);
		std::vector<NearestTipTargets> groups;
		for (int g = 0; g < n_groups; ++g) {
			const double fraction = target_fractions[(trial + g) % 4];
			for (uint32_t v : tips) {
				if (unit(rng) < fraction) {
					group_nodes[g].push_back(v);
					group_ranks[g].push_back(tip_rank[v]);
				}
			}
			groups.emplace_back(index, group_nodes[g], group_ranks[g]);
		}

		for (uint32_t q : tips) {
			for (int g = 0; g < n_groups; ++g) {
				for (auto scope : {NearestTipScope::Path, NearestTipScope::Clade}) {
					for (bool exclude_self : {true, false}) {
						CAPTURE(trial, n, dyadic, g, q, scope == NearestTipScope::Clade, exclude_self);
						auto expected = oracle(t, group_nodes[g], group_ranks[g], q, exclude_self, scope);
						NearestTipOptions options;
						options.scope = scope;
						options.exclude_self = exclude_self;
						auto got = groups[g].Search(q, options);
						REQUIRE(got.target == expected.target);
						REQUIRE(got.n_tied == expected.n_tied);
						if (expected.target != NO_TARGET) {
							REQUIRE(got.distance == expected.distance);
							REQUIRE(got.lca == expected.lca);
						}
					}
				}
			}
		}
	}
}

TEST_CASE("NearestTip: rejects malformed target sets", "[NearestTip]") {
	NamedTree t("((A:1,B:1)c:1,D:1)r;");
	const uint32_t a = t.node("A");
	const uint32_t b = t.node("B");
	REQUIRE_THROWS_WITH(NearestTipTargets(t.index, {a, b}, {0}), ContainsSubstring("size"));
	REQUIRE_THROWS_WITH(NearestTipTargets(t.index, {a, a}, {0, 1}), ContainsSubstring("twice"));
	REQUIRE_THROWS_WITH(NearestTipTargets(t.index, {a, 99}, {0, 1}), ContainsSubstring("out of range"));
}
