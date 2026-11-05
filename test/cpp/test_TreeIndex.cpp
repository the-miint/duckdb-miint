#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "NewickTree.hpp"
#include "TreeIndex.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <string>
#include <vector>

using Catch::Matchers::ContainsSubstring;
using miint::TreeIndex;

namespace {

constexpr uint32_t NP = TreeIndex::NO_PARENT;
const double kNaN = std::numeric_limits<double>::quiet_NaN();

// A tree as plain arrays, the input TreeIndex::Build takes.
struct Arrays {
	std::vector<uint32_t> parent;
	std::vector<double> length;
};

uint32_t node_named(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

// ---------------------------------------------------------------------------
// The naive model the index is checked against: everything by walking parents.
// ---------------------------------------------------------------------------

bool naive_is_ancestor_or_self(const Arrays &t, uint32_t a, uint32_t v) {
	for (uint32_t cur = v; cur != NP; cur = t.parent[cur]) {
		if (cur == a) {
			return true;
		}
	}
	return false;
}

uint32_t naive_lca(const Arrays &t, uint32_t a, uint32_t b) {
	std::vector<bool> above_a(t.parent.size(), false);
	for (uint32_t cur = a; cur != NP; cur = t.parent[cur]) {
		above_a[cur] = true;
	}
	for (uint32_t cur = b; cur != NP; cur = t.parent[cur]) {
		if (above_a[cur]) {
			return cur;
		}
	}
	FAIL("no common ancestor");
	return NP;
}

uint32_t naive_depth(const Arrays &t, uint32_t v) {
	uint32_t d = 0;
	for (uint32_t cur = v; t.parent[cur] != NP; cur = t.parent[cur]) {
		++d;
	}
	return d;
}

// The root distance, summed root first: the order the index promises, so even
// non-dyadic lengths must match bit for bit.
double naive_root_distance(const Arrays &t, uint32_t v) {
	std::vector<uint32_t> path; // v up to (excluding) the root
	for (uint32_t cur = v; t.parent[cur] != NP; cur = t.parent[cur]) {
		path.push_back(cur);
	}
	double d = 0.0;
	for (auto it = path.rbegin(); it != path.rend(); ++it) {
		d += t.length[*it];
	}
	return d;
}

// ---------------------------------------------------------------------------
// Random trees. Structure first (each node's parent is an earlier node, and often
// the immediately preceding one, so long chains appear), then the labels are
// shuffled so the numbering is in no particular order.
// ---------------------------------------------------------------------------

Arrays random_tree(std::mt19937 &rng, uint32_t n, bool dyadic) {
	std::vector<uint32_t> structural_parent(n, NP);
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	for (uint32_t k = 1; k < n; ++k) {
		if (unit(rng) < 0.3) {
			structural_parent[k] = k - 1;
		} else {
			structural_parent[k] = std::uniform_int_distribution<uint32_t>(0, k - 1)(rng);
		}
	}
	std::vector<uint32_t> label(n);
	std::iota(label.begin(), label.end(), 0u);
	std::shuffle(label.begin(), label.end(), rng);

	const double dyadic_lengths[] = {0.0, 0.0, 0.125, 0.25, 0.5, 1.0, 2.0};
	Arrays t;
	t.parent.assign(n, NP);
	t.length.assign(n, 0.0);
	for (uint32_t k = 0; k < n; ++k) {
		if (structural_parent[k] != NP) {
			t.parent[label[k]] = label[structural_parent[k]];
		}
		t.length[label[k]] = dyadic ? dyadic_lengths[std::uniform_int_distribution<int>(0, 6)(rng)] : unit(rng);
	}
	// The root's own length is ignored, whatever it holds.
	t.length[label[0]] = kNaN;
	return t;
}

// A caterpillar `levels` deep: internal node k (index k) is the parent of internal node
// k+1 and of tip k (index levels + k); the deepest internal node also holds the last tip.
// Every edge has length 1, so root distances are exact integers.
Arrays caterpillar(uint32_t levels) {
	const uint32_t n = 2 * levels + 1;
	Arrays t;
	t.parent.assign(n, NP);
	t.length.assign(n, 1.0);
	for (uint32_t k = 1; k < levels; ++k) {
		t.parent[k] = k - 1;
	}
	for (uint32_t k = 0; k < levels; ++k) {
		t.parent[levels + k] = k;
	}
	t.parent[2 * levels] = levels - 1;
	return t;
}

void require_matches_naive(const TreeIndex &idx, const Arrays &t) {
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
