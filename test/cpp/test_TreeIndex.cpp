#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "NewickTree.hpp"
#include "TreeIndex.hpp"
#include "tree_test_util.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <random>
#include <string>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::TreeIndex;
using namespace miint::test;

namespace {

constexpr uint32_t NP = TreeIndex::NO_PARENT;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

uint32_t node_named(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

void require_matches_naive(const TreeIndex &idx, const TreeArrays &t) {
	const uint32_t n = static_cast<uint32_t>(t.parent.size());
	REQUIRE(idx.size() == n);

	// Preorder positions are a permutation of 0..n-1.
	std::vector<bool> seen(n, false);
	for (uint32_t v = 0; v < n; ++v) {
		REQUIRE(idx.pre(v) < n);
		REQUIRE_FALSE(seen[idx.pre(v)]);
		seen[idx.pre(v)] = true;
	}

	for (uint32_t v = 0; v < n; ++v) {
		REQUIRE(idx.parent(v) == t.parent[v]);
		REQUIRE(idx.depth(v) == naive_depth(t, v));
		REQUIRE(idx.root_distance(v) == naive_root_distance(t, v));
	}
	for (uint32_t a = 0; a < n; ++a) {
		for (uint32_t v = 0; v < n; ++v) {
			// Every subtree is exactly one preorder interval.
			REQUIRE(idx.is_ancestor_or_self(a, v) == naive_is_ancestor_or_self(t, a, v));
			REQUIRE(idx.lca(a, v) == naive_lca(t, a, v));
		}
	}
}

} // namespace

TEST_CASE("TreeIndex: root distance and depth on a hand-written tree", "[TreeIndex]") {
	// Dyadic lengths, so every expected sum is exact. G's own length (3) is ignored.
	auto tree = miint::NewickTree::parse("((A:0.5,B:0.25)C:1,(D:0.125,E:2)F:0.75)G:3;");
	auto idx = TreeIndex::FromNewickTree(tree);

	REQUIRE(idx.root() == node_named(tree, "G"));
	REQUIRE(idx.root_distance(node_named(tree, "G")) == 0.0);
	REQUIRE(idx.root_distance(node_named(tree, "C")) == 1.0);
	REQUIRE(idx.root_distance(node_named(tree, "A")) == 1.5);
	REQUIRE(idx.root_distance(node_named(tree, "B")) == 1.25);
	REQUIRE(idx.root_distance(node_named(tree, "F")) == 0.75);
	REQUIRE(idx.root_distance(node_named(tree, "D")) == 0.875);
	REQUIRE(idx.root_distance(node_named(tree, "E")) == 2.75);

	REQUIRE(idx.depth(node_named(tree, "G")) == 0);
	REQUIRE(idx.depth(node_named(tree, "C")) == 1);
	REQUIRE(idx.depth(node_named(tree, "E")) == 2);

	REQUIRE(idx.lca(node_named(tree, "A"), node_named(tree, "B")) == node_named(tree, "C"));
	REQUIRE(idx.lca(node_named(tree, "A"), node_named(tree, "E")) == node_named(tree, "G"));
	REQUIRE(idx.lca(node_named(tree, "D"), node_named(tree, "F")) == node_named(tree, "F"));
	REQUIRE(idx.lca(node_named(tree, "B"), node_named(tree, "B")) == node_named(tree, "B"));
}

TEST_CASE("TreeIndex: lowest_ancestor_containing climbs to the first subtree holding a position", "[TreeIndex]") {
	auto tree = miint::NewickTree::parse("((A:1,B:1)C:1,(D:1,(E:1,F:1)G:1)H:1)R;");
	auto idx = TreeIndex::FromNewickTree(tree);
	auto at = [&](const char *name) {
		return node_named(tree, name);
	};

	// From E: its own subtree holds only E; G adds F; H adds D; R adds A, B, C.
	REQUIRE(idx.lowest_ancestor_containing(at("E"), idx.pre(at("E"))) == at("E"));
	REQUIRE(idx.lowest_ancestor_containing(at("E"), idx.pre(at("F"))) == at("G"));
	REQUIRE(idx.lowest_ancestor_containing(at("E"), idx.pre(at("D"))) == at("H"));
	REQUIRE(idx.lowest_ancestor_containing(at("E"), idx.pre(at("A"))) == at("R"));
	// Starting part-way up is the same climb.
	REQUIRE(idx.lowest_ancestor_containing(at("G"), idx.pre(at("D"))) == at("H"));
	REQUIRE(idx.lowest_ancestor_containing(at("H"), idx.pre(at("F"))) == at("H"));
}

TEST_CASE("TreeIndex: random trees in arbitrary numbering agree with a naive walk", "[TreeIndex]") {
	std::mt19937 rng(20261105);
	for (int trial = 0; trial < 60; ++trial) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 120)(rng);
		const bool dyadic = trial % 2 == 0;
		auto t = random_tree(rng, n, dyadic);
		CAPTURE(trial, n, dyadic);
		auto idx = TreeIndex::Build(t.parent, t.length);
		require_matches_naive(idx, t);
	}
}

TEST_CASE("TreeIndex: a single node is a tree", "[TreeIndex]") {
	auto idx = TreeIndex::Build({NP}, {kNaN});
	REQUIRE(idx.size() == 1);
	REQUIRE(idx.root() == 0);
	REQUIRE(idx.root_distance(0) == 0.0);
	REQUIRE(idx.lca(0, 0) == 0);
}

TEST_CASE("TreeIndex: a 200,000-level caterpillar builds and answers without recursion", "[TreeIndex]") {
	// Placement trees nest this deep. A recursive build would overflow the stack, and a
	// parent-by-parent climb would be slow; the answers here are known in closed form.
	const uint32_t levels = 200000;
	auto t = caterpillar(levels);
	auto idx = TreeIndex::Build(t.parent, t.length);

	auto tip = [&](uint32_t k) {
		return levels + k;
	};
	REQUIRE(idx.depth(tip(levels)) == levels);
	REQUIRE(idx.root_distance(tip(levels)) == static_cast<double>(levels));
	REQUIRE(idx.depth(tip(0)) == 1);

	std::mt19937 rng(7);
	std::uniform_int_distribution<uint32_t> pick(0, levels);
	for (int i = 0; i < 2000; ++i) {
		uint32_t a = pick(rng);
		uint32_t b = pick(rng);
		// Tips k and m (k < m) meet at internal node k; the last tip shares node levels-1.
		uint32_t expected = std::min(std::min(a, b), levels - 1);
		CAPTURE(a, b);
		REQUIRE(idx.lca(tip(a), tip(b)) == expected);
	}
}

TEST_CASE("TreeIndex: the root's branch length is ignored", "[TreeIndex]") {
	for (double root_length : {kNaN, -5.0, std::numeric_limits<double>::infinity(), 7.0}) {
		auto idx = TreeIndex::Build({NP, 0, 0}, {root_length, 1.0, 2.0});
		REQUIRE(idx.root_distance(0) == 0.0);
		REQUIRE(idx.root_distance(2) == 2.0);
	}
}

TEST_CASE("TreeIndex: rejects input that is not a rooted tree", "[TreeIndex]") {
	SECTION("empty") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({}, {}), ContainsSubstring("empty"));
	}
	SECTION("array sizes differ") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0}, {1.0}), ContainsSubstring("size"));
	}
	SECTION("two roots") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, NP}, {1.0, 1.0}), ContainsSubstring("root"));
	}
	SECTION("no root") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({1, 0}, {1.0, 1.0}), ContainsSubstring("root"));
	}
	SECTION("a parent out of range") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 5}, {1.0, 1.0}), ContainsSubstring("out of range"));
	}
	SECTION("a cycle beside the root") {
		// 0 is the root; 1 and 2 point at each other and never reach it.
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 2, 1}, {1.0, 1.0, 1.0}), ContainsSubstring("cycle"));
	}
	SECTION("a node that is its own parent") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 1}, {1.0, 1.0}), ContainsSubstring("cycle"));
	}
}

TEST_CASE("TreeIndex: rejects a NULL, infinite or negative branch length below the root", "[TreeIndex]") {
	// Root distances must grow away from the root, or a search that prunes by distance
	// gives wrong answers; a missing length has no distance at all.
	std::vector<int64_t> ids = {1000, 1001, 1002};
	SECTION("NULL") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0, 0}, {1.0, 1.0, kNaN}, &ids),
		                    ContainsSubstring("node_index 1002") && ContainsSubstring("NULL"));
	}
	SECTION("negative") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0, 0}, {1.0, -0.5, 1.0}, &ids),
		                    ContainsSubstring("node_index 1001") && ContainsSubstring("negative"));
	}
	SECTION("infinite") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0, 0}, {1.0, std::numeric_limits<double>::infinity(), 1.0}, &ids),
		                    ContainsSubstring("node_index 1001") && ContainsSubstring("finite"));
	}
	SECTION("without caller ids the dense index is named") {
		REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0, 0}, {1.0, 1.0, kNaN}), ContainsSubstring("node 2 "));
	}
}

TEST_CASE("TreeIndex: FromNewickTree applies the same validation", "[TreeIndex]") {
	// B has no branch length.
	auto tree = miint::NewickTree::parse("(A:1,B)R;");
	REQUIRE_THROWS_WITH(TreeIndex::FromNewickTree(tree), ContainsSubstring("NULL"));
}
