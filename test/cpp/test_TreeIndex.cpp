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
		REQUIRE(idx.is_tip(v) == naive_is_tip(t, v));
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
	auto idx = index_of(tree);

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
	auto idx = index_of(tree);
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
	// parent-by-parent climb would take O(depth) steps; the answers here are known in
	// closed form.
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
		// A tip is its own LCA.
		uint32_t expected = a == b ? tip(a) : std::min(std::min(a, b), levels - 1);
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
		                    ContainsSubstring("node_index 1002") && ContainsSubstring("NULL or NaN"));
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

TEST_CASE("TreeIndex: rejects root distances that overflow", "[TreeIndex]") {
	// Each length is finite, but their sum is not; a distance from it would be infinite or NaN.
	const double big = std::numeric_limits<double>::max() / 2 * 1.5;
	REQUIRE_THROWS_WITH(TreeIndex::Build({NP, 0, 1}, {kNaN, big, big}), ContainsSubstring("node 2") &&
	                                                                        ContainsSubstring("root distance") &&
	                                                                        ContainsSubstring("not finite"));
}

// ---------------------------------------------------------------------------
// DenseTree: a tree table's rows, in any order and with any distinct ids, numbered
// densely for TreeIndex::Build. Callers number nodes however they like (read_newick's
// postorder, an offset, a sparse key) and store rows in any order, so neither may change
// the tree that is built.
// ---------------------------------------------------------------------------

namespace {

using miint::DenseTree;

struct Rows {
	std::vector<int64_t> id;
	std::vector<int64_t> parent_id;
	std::vector<bool> has_parent;
	std::vector<double> length;
};

DenseTree dense(Rows rows) {
	return DenseTree::FromRows(std::move(rows.id), std::move(rows.parent_id), std::move(rows.has_parent),
	                           std::move(rows.length));
}

// The rows of tree t with node v given id ids[v], stored in the order `order`.
Rows rows_of(const TreeArrays &t, const std::vector<int64_t> &ids, const std::vector<uint32_t> &order) {
	Rows r;
	for (uint32_t v : order) {
		r.id.push_back(ids[v]);
		r.has_parent.push_back(t.parent[v] != NP);
		r.parent_id.push_back(t.parent[v] != NP ? ids[t.parent[v]] : 0);
		r.length.push_back(t.length[v]);
	}
	return r;
}

} // namespace

TEST_CASE("DenseTree: contiguous ids in order keep their order", "[TreeIndex]") {
	// R(1003) over C(1002), C over A(1000) and B(1001): read_newick's postorder, offset by 1000.
	auto d = dense({{1000, 1001, 1002, 1003}, {1002, 1002, 1003, 0}, {true, true, true, false}, {1.0, 2.0, 3.0, kNaN}});
	REQUIRE(d.node_ids == std::vector<int64_t> {1000, 1001, 1002, 1003});
	REQUIRE(d.parent == std::vector<uint32_t> {2, 2, 3, NP});
	REQUIRE(d.branch_length[2] == 3.0);
	REQUIRE(std::isnan(d.branch_length[3]));
	REQUIRE(d.find(1001) == 1);
	REQUIRE(d.find(999) == DenseTree::NOT_FOUND);
	REQUIRE(d.find(1004) == DenseTree::NOT_FOUND);
}

TEST_CASE("DenseTree: rows out of order with sparse ids are numbered by ascending id", "[TreeIndex]") {
	// The same tree as above, ids 40 (A), 10 (B), 30 (C), -5 (R), rows in no order. Dense
	// order is ascending id: R, B, C, A.
	auto d = dense({{30, 40, -5, 10}, {-5, 30, 0, 30}, {true, true, false, true}, {3.0, 1.0, kNaN, 2.0}});
	REQUIRE(d.node_ids == std::vector<int64_t> {-5, 10, 30, 40});
	REQUIRE(d.parent == std::vector<uint32_t> {NP, 2, 0, 2});
	REQUIRE(std::isnan(d.branch_length[0]));
	REQUIRE(d.branch_length[1] == 2.0);
	REQUIRE(d.branch_length[2] == 3.0);
	REQUIRE(d.branch_length[3] == 1.0);
	REQUIRE(d.find(30) == 2);
	REQUIRE(d.find(20) == DenseTree::NOT_FOUND);
	REQUIRE(d.find(-6) == DenseTree::NOT_FOUND);
	REQUIRE(d.find(41) == DenseTree::NOT_FOUND);
}

TEST_CASE("DenseTree: ids at both ends of the int64 range", "[TreeIndex]") {
	// The id span overflows int64, so a contiguity test done in signed arithmetic would be
	// undefined; three ids that are far from contiguous must not be taken for contiguous.
	const int64_t lo = std::numeric_limits<int64_t>::min();
	const int64_t hi = std::numeric_limits<int64_t>::max();
	auto d = dense({{hi, lo, 0}, {0, 0, 0}, {true, true, false}, {1.0, 2.0, kNaN}});
	REQUIRE(d.node_ids == std::vector<int64_t> {lo, 0, hi});
	REQUIRE(d.parent == std::vector<uint32_t> {1, NP, 1});
	REQUIRE(d.find(lo) == 0);
	REQUIRE(d.find(hi) == 2);
	REQUIRE(d.find(1) == DenseTree::NOT_FOUND);
	REQUIRE(d.find(lo + 1) == DenseTree::NOT_FOUND);
}

TEST_CASE("DenseTree: any row order and any ids give the same tree", "[TreeIndex]") {
	std::mt19937 rng(20261106);
	for (int trial = 0; trial < 40; ++trial) {
		const uint32_t n = std::uniform_int_distribution<uint32_t>(1, 150)(rng);
		auto t = random_tree(rng, n, trial % 2 == 0);
		// Distinct ids: contiguous from an offset, or sparse and signed.
		std::vector<int64_t> ids(n);
		const bool contiguous = trial % 3 == 0;
		for (uint32_t v = 0; v < n; ++v) {
			ids[v] = contiguous ? 500 + v : static_cast<int64_t>(v) * 7919 - 100000;
		}
		std::shuffle(ids.begin(), ids.end(), rng);
		std::vector<uint32_t> order(n);
		std::iota(order.begin(), order.end(), 0u);
		std::shuffle(order.begin(), order.end(), rng);
		CAPTURE(trial, n, contiguous);

		auto d = dense(rows_of(t, ids, order));
		REQUIRE(std::is_sorted(d.node_ids.begin(), d.node_ids.end()));
		for (uint32_t v = 0; v < n; ++v) {
			const uint32_t dv = d.find(ids[v]);
			REQUIRE(dv < n);
			REQUIRE(d.node_ids[dv] == ids[v]);
			if (t.parent[v] == NP) {
				REQUIRE(d.parent[dv] == NP);
			} else {
				REQUIRE(d.parent[dv] == d.find(ids[t.parent[v]]));
				REQUIRE(d.branch_length[dv] == t.length[v]);
			}
		}
		// And it is a tree the index accepts, with the same root distances.
		auto idx = TreeIndex::Build(d.parent, d.branch_length, &d.node_ids);
		for (uint32_t v = 0; v < n; ++v) {
			REQUIRE(idx.root_distance(d.find(ids[v])) == naive_root_distance(t, v));
		}
	}
}

TEST_CASE("DenseTree: rejects rows that cannot be numbered", "[TreeIndex]") {
	SECTION("no rows") {
		REQUIRE_THROWS_WITH(dense({}), ContainsSubstring("empty"));
	}
	SECTION("columns of different lengths") {
		REQUIRE_THROWS_WITH(dense({{1, 2}, {0}, {false, true}, {1.0, 1.0}}), ContainsSubstring("size"));
	}
	SECTION("a repeated id, out of order") {
		REQUIRE_THROWS_WITH(dense({{7, 3, 7}, {0, 7, 3}, {false, true, true}, {kNaN, 1.0, 1.0}}),
		                    ContainsSubstring("node_index 7 appears more than once"));
	}
	SECTION("a parent id that is no row's id") {
		REQUIRE_THROWS_WITH(dense({{0, 1, 2}, {2, 99, 0}, {true, true, false}, {1.0, 1.0, kNaN}}),
		                    ContainsSubstring("parent_index of node_index 1 (99)"));
	}
	SECTION("a parent id outside a contiguous range") {
		REQUIRE_THROWS_WITH(dense({{0, 1, 2}, {3, 2, 0}, {true, true, false}, {1.0, 1.0, kNaN}}),
		                    ContainsSubstring("parent_index of node_index 0 (3)"));
	}
}
