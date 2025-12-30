// Tests for insert_fully_resolved phylogenetic placement algorithm
// and supporting distance calculation methods
//
// Test cases inspired by improved-octo-waddle (BSD-3-Clause License)
// https://github.com/biocore/improved-octo-waddle
// bp/tests/test_insert.py

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <algorithm>
#include <cmath>
#include <set>
#include <string>
#include <vector>
#include "NewickTree.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

// ============================================================================
// Tests for distance calculation methods
// ============================================================================

TEST_CASE("NewickTree tips() returns leaf nodes", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A:1.0,B:2.0):0.5,C:3.0);");

	auto tip_indices = tree.tips();
	REQUIRE(tip_indices.size() == 3);

	std::set<std::string> names;
	for (auto idx : tip_indices) {
		names.insert(tree.name(idx));
	}
	REQUIRE(names == std::set<std::string>{"A", "B", "C"});
}

TEST_CASE("NewickTree tip_names() returns leaf names", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A:1.0,B:2.0):0.5,C:3.0);");

	auto names = tree.tip_names();
	REQUIRE(names.size() == 3);

	std::set<std::string> name_set(names.begin(), names.end());
	REQUIRE(name_set == std::set<std::string>{"A", "B", "C"});
}

TEST_CASE("NewickTree distance_to_root", "[NewickTree][distance]") {
	// Tree: ((A:1.0,B:2.0):0.5,C:3.0):0.0;
	//        root (0.0)
	//       /    \
	//   internal   C (3.0)
	//    (0.5)
	//   /    \
	//  A      B
	// (1.0)  (2.0)

	auto tree = miint::NewickTree::parse("((A:1.0,B:2.0):0.5,C:3.0):0.0;");

	auto a = tree.find_node_by_name("A").value();
	auto b = tree.find_node_by_name("B").value();
	auto c = tree.find_node_by_name("C").value();

	// A: 1.0 + 0.5 = 1.5
	REQUIRE(tree.distance_to_root(a) == Approx(1.5));

	// B: 2.0 + 0.5 = 2.5
	REQUIRE(tree.distance_to_root(b) == Approx(2.5));

	// C: 3.0
	REQUIRE(tree.distance_to_root(c) == Approx(3.0));

	// Root itself: 0
	REQUIRE(tree.distance_to_root(tree.root()) == Approx(0.0));
}

TEST_CASE("NewickTree distance_to_root handles NaN branch lengths", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A,B),C);");

	auto a = tree.find_node_by_name("A").value();
	// All branch lengths are NaN, so distance is 0
	REQUIRE(tree.distance_to_root(a) == Approx(0.0));
}

TEST_CASE("NewickTree find_lca", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A:1.0,B:2.0)AB:0.5,(C:1.0,D:2.0)CD:0.5)root;");

	auto a = tree.find_node_by_name("A").value();
	auto b = tree.find_node_by_name("B").value();
	auto c = tree.find_node_by_name("C").value();
	auto ab = tree.find_node_by_name("AB").value();

	// LCA of A and B is AB
	REQUIRE(tree.find_lca(a, b) == ab);

	// LCA of A and C is root
	REQUIRE(tree.find_lca(a, c) == tree.root());

	// LCA of A and A is A
	REQUIRE(tree.find_lca(a, a) == a);

	// LCA of A and AB is AB
	REQUIRE(tree.find_lca(a, ab) == ab);
}

TEST_CASE("NewickTree pairwise_distance", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A:1.0,B:2.0):0.5,C:3.0):0.0;");

	auto a = tree.find_node_by_name("A").value();
	auto b = tree.find_node_by_name("B").value();
	auto c = tree.find_node_by_name("C").value();

	// Distance A to B: 1.0 + 2.0 = 3.0
	REQUIRE(tree.pairwise_distance(a, b) == Approx(3.0));

	// Distance A to C: 1.0 + 0.5 + 3.0 = 4.5
	REQUIRE(tree.pairwise_distance(a, c) == Approx(4.5));

	// Distance B to C: 2.0 + 0.5 + 3.0 = 5.5
	REQUIRE(tree.pairwise_distance(b, c) == Approx(5.5));

	// Distance to self is 0
	REQUIRE(tree.pairwise_distance(a, a) == Approx(0.0));
}

TEST_CASE("NewickTree build_edge_index", "[NewickTree][distance]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	auto index = tree.build_edge_index();

	REQUIRE(index.size() == 5);
	REQUIRE(tree.name(index.at(0)) == "A");
	REQUIRE(tree.name(index.at(1)) == "B");
	REQUIRE(tree.name(index.at(3)) == "C");
}

// ============================================================================
// Tests for insert_fully_resolved - basic functionality
// ============================================================================

TEST_CASE("insert_fully_resolved single placement", "[insert][basic]") {
	// Tree: ((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};
	// Insert fragment F on edge 0 (A's edge) with distal_length=0.3, pendant_length=0.1
	//
	// Expected result:
	//   The edge to A (length 1.0) is split at 0.3 from A
	//   New internal node at distance 0.7 from parent, 0.3 to A
	//   Fragment F hangs off with pendant_length 0.1

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	REQUIRE(tree.num_tips() == 3);
	REQUIRE(tree.num_nodes() == 5);

	std::vector<miint::Placement> placements = {{.fragment_id = "F",
	                                             .edge_id = 0,
	                                             .distal_length = 0.3,
	                                             .pendant_length = 0.1,
	                                             .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	// Should have added 2 nodes: new internal + fragment F
	REQUIRE(tree.num_nodes() == 7);
	REQUIRE(tree.num_tips() == 4);

	// Check F exists
	auto f_node = tree.find_node_by_name("F");
	REQUIRE(f_node.has_value());
	REQUIRE(tree.is_tip(f_node.value()));
	REQUIRE(tree.branch_length(f_node.value()) == Approx(0.1));

	// Check A's branch length is now 0.3 (distal_length)
	auto a_node = tree.find_node_by_name("A");
	REQUIRE(a_node.has_value());
	REQUIRE(tree.branch_length(a_node.value()) == Approx(0.3));

	// Check the new internal node's branch length is 0.7 (1.0 - 0.3)
	uint32_t a_parent = tree.parent(a_node.value());
	REQUIRE(tree.branch_length(a_parent) == Approx(0.7));

	// F should be sibling of A
	uint32_t f_parent = tree.parent(f_node.value());
	REQUIRE(f_parent == a_parent);
}

TEST_CASE("insert_fully_resolved preserves tree distances", "[insert][distance]") {
	// Verify that original tip-to-tip distances are preserved after insertion
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	auto a = tree.find_node_by_name("A").value();
	auto b = tree.find_node_by_name("B").value();
	auto c = tree.find_node_by_name("C").value();

	double orig_ab = tree.pairwise_distance(a, b);
	double orig_ac = tree.pairwise_distance(a, c);
	double orig_bc = tree.pairwise_distance(b, c);

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F1", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    {.fragment_id = "F2", .edge_id = 1, .distal_length = 0.5, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	// Re-find nodes (indices may have changed)
	a = tree.find_node_by_name("A").value();
	b = tree.find_node_by_name("B").value();
	c = tree.find_node_by_name("C").value();

	// Original distances should be preserved
	REQUIRE(tree.pairwise_distance(a, b) == Approx(orig_ab));
	REQUIRE(tree.pairwise_distance(a, c) == Approx(orig_ac));
	REQUIRE(tree.pairwise_distance(b, c) == Approx(orig_bc));
}

TEST_CASE("insert_fully_resolved multiple placements same edge", "[insert][multi]") {
	// Insert multiple fragments on the same edge
	// They should form a chain sorted by distal_length (descending)
	//
	// Edge to A has length 1.0
	// Insert F1 at distal=0.7, F2 at distal=0.3
	// After sorting by distal desc: F1 (0.7), F2 (0.3)
	//
	// Result chain from original parent toward A:
	//   original_parent --0.3--> new_node_1 --0.4--> new_node_2 --0.3--> A
	//                              |                    |
	//                              F1 (0.1)             F2 (0.2)

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F1", .edge_id = 0, .distal_length = 0.7, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    {.fragment_id = "F2", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	// Should have added 4 nodes: 2 internals + 2 fragments
	REQUIRE(tree.num_nodes() == 9);
	REQUIRE(tree.num_tips() == 5);

	// Verify A's final branch length is 0.3 (smallest distal_length)
	auto a = tree.find_node_by_name("A").value();
	REQUIRE(tree.branch_length(a) == Approx(0.3));

	// Verify F2 is sibling of A (both at distal=0.3 level)
	auto f2 = tree.find_node_by_name("F2").value();
	REQUIRE(tree.parent(f2) == tree.parent(a));

	// Verify F1's branch length
	auto f1 = tree.find_node_by_name("F1").value();
	REQUIRE(tree.branch_length(f1) == Approx(0.1));

	// Verify F2's branch length
	REQUIRE(tree.branch_length(f2) == Approx(0.2));
}

TEST_CASE("insert_fully_resolved deduplicates by fragment_id", "[insert][dedup]") {
	// When same fragment has multiple placements, keep best one
	// Best = highest like_weight_ratio, then lowest pendant_length

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    // Lower quality placement (lower like_weight_ratio)
	    {.fragment_id = "F", .edge_id = 1, .distal_length = 0.5, .pendant_length = 0.1, .like_weight_ratio = 0.5},
	    // Higher quality placement (higher like_weight_ratio)
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	// Only one F should be inserted
	REQUIRE(tree.num_tips() == 4);

	auto f = tree.find_node_by_name("F").value();
	// Should have pendant_length from the better placement
	REQUIRE(tree.branch_length(f) == Approx(0.2));

	// Should be on edge 0 (A's edge), not edge 1
	auto a = tree.find_node_by_name("A").value();
	// F and A should share same parent (the new internal node on edge 0)
	REQUIRE(tree.parent(f) == tree.parent(a));
}

TEST_CASE("insert_fully_resolved dedup tiebreaker pendant_length", "[insert][dedup]") {
	// Same like_weight_ratio: prefer lower pendant_length

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.5, .like_weight_ratio = 1.0},
	    {.fragment_id = "F", .edge_id = 1, .distal_length = 0.5, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	auto f = tree.find_node_by_name("F").value();
	// Should have the lower pendant_length
	REQUIRE(tree.branch_length(f) == Approx(0.1));
}

TEST_CASE("insert_fully_resolved on root edge", "[insert][root]") {
	// Insert on the root's edge (edge_id for root itself)

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):1.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 4, .distal_length = 0.4, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 4);

	auto f = tree.find_node_by_name("F").value();
	REQUIRE(tree.branch_length(f) == Approx(0.1));
}

TEST_CASE("insert_fully_resolved empty placements", "[insert][edge]") {
	// No placements - tree should be unchanged

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");
	size_t orig_nodes = tree.num_nodes();

	std::vector<miint::Placement> placements;
	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_nodes() == orig_nodes);
}

TEST_CASE("insert_fully_resolved invalid edge_id throws", "[insert][error]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 999, .distal_length = 0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	REQUIRE_THROWS_WITH(tree.insert_fully_resolved(placements), ContainsSubstring("edge"));
}

TEST_CASE("insert_fully_resolved distal_length exceeds edge throws", "[insert][error]") {
	// Edge 0 has length 1.0, distal_length 1.5 is invalid

	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 1.5, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	REQUIRE_THROWS_WITH(tree.insert_fully_resolved(placements), ContainsSubstring("distal_length"));
}

// ============================================================================
// Scale tests
// ============================================================================

TEST_CASE("insert_fully_resolved many placements on single edge", "[insert][scale]") {
	// Insert 100 placements on the same edge
	auto tree = miint::NewickTree::parse("((A:100.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements;
	placements.reserve(100);
	for (int i = 0; i < 100; ++i) {
		placements.push_back({.fragment_id = "F" + std::to_string(i),
		                      .edge_id = 0,
		                      .distal_length = static_cast<double>(i + 1), // 1.0 to 100.0
		                      .pendant_length = 0.1,
		                      .like_weight_ratio = 1.0});
	}

	tree.insert_fully_resolved(placements);

	// 3 original tips + 100 fragments
	REQUIRE(tree.num_tips() == 103);
	// 5 original nodes + 200 new (100 internals + 100 fragments)
	REQUIRE(tree.num_nodes() == 205);
}

TEST_CASE("insert_fully_resolved many placements on different edges", "[insert][scale]") {
	// Build a larger tree and insert placements on multiple edges
	// Edge lengths: A,B,C,D = 1.0; internal edges = 0.5
	auto tree = miint::NewickTree::parse(
	    "((A:1.0{0},B:1.0{1}):0.5{2},(C:1.0{3},D:1.0{4}):0.5{5}):0.0{6};");

	std::vector<miint::Placement> placements;
	int frag_id = 0;
	// Only place on tip edges (0, 1, 3, 4) which have length 1.0
	std::vector<int64_t> tip_edges = {0, 1, 3, 4};
	for (int64_t edge : tip_edges) {
		for (int i = 0; i < 10; ++i) {
			placements.push_back({.fragment_id = "F" + std::to_string(frag_id++),
			                      .edge_id = edge,
			                      .distal_length = 0.1 * (i + 1),
			                      .pendant_length = 0.05,
			                      .like_weight_ratio = 1.0});
		}
	}

	tree.insert_fully_resolved(placements);

	// 4 original tips + 40 fragments
	REQUIRE(tree.num_tips() == 44);
}

TEST_CASE("insert_fully_resolved large scale", "[insert][scale]") {
	// Test with 10,000 placements to verify performance
	// Build a simple backbone tree
	auto tree = miint::NewickTree::parse("((A:1000.0{0},B:1000.0{1}):500.0{2},C:1000.0{3}):0.0{4};");

	std::vector<miint::Placement> placements;
	placements.reserve(10000);

	// Distribute placements across edges 0, 1, 3 (tips)
	for (int i = 0; i < 10000; ++i) {
		int64_t edge = (i % 3 == 0) ? 0 : (i % 3 == 1) ? 1 : 3;
		placements.push_back({.fragment_id = "F" + std::to_string(i),
		                      .edge_id = edge,
		                      .distal_length = static_cast<double>(i % 900 + 1),
		                      .pendant_length = 0.01,
		                      .like_weight_ratio = 1.0});
	}

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 10003);
}

// ============================================================================
// Serialization after insertion
// ============================================================================

TEST_CASE("insert_fully_resolved roundtrip serialization", "[insert][serialize]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F1", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    {.fragment_id = "F2", .edge_id = 1, .distal_length = 0.5, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	// Serialize and re-parse
	std::string newick = tree.to_newick();
	auto tree2 = miint::NewickTree::parse(newick);

	// Same structure
	REQUIRE(tree2.num_nodes() == tree.num_nodes());
	REQUIRE(tree2.num_tips() == tree.num_tips());

	// Same tip names
	auto names1 = tree.tip_names();
	auto names2 = tree2.tip_names();
	std::set<std::string> set1(names1.begin(), names1.end());
	std::set<std::string> set2(names2.begin(), names2.end());
	REQUIRE(set1 == set2);
}

// ============================================================================
// Edge case tests (from code review)
// ============================================================================

TEST_CASE("insert_fully_resolved negative pendant_length throws", "[insert][error]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.3, .pendant_length = -0.1, .like_weight_ratio = 1.0}};

	REQUIRE_THROWS_WITH(tree.insert_fully_resolved(placements), ContainsSubstring("pendant_length"));
}

TEST_CASE("insert_fully_resolved negative distal_length throws", "[insert][error]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = -0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	REQUIRE_THROWS_WITH(tree.insert_fully_resolved(placements), ContainsSubstring("distal_length"));
}

TEST_CASE("insert_fully_resolved validates non-winning placements", "[insert][error]") {
	// Ensure validation catches errors even in placements that would be filtered by dedup
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    // This one would "win" deduplication (higher like_weight_ratio)
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    // This one would "lose" but has invalid negative pendant_length - should still be caught
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.5, .pendant_length = -0.1, .like_weight_ratio = 0.5}};

	REQUIRE_THROWS_WITH(tree.insert_fully_resolved(placements), ContainsSubstring("pendant_length"));
}

TEST_CASE("insert_fully_resolved on edge with NaN length", "[insert][nan]") {
	// Tree with NaN branch length (no length specified)
	auto tree = miint::NewickTree::parse("((A{0},B{1}){2},C{3}){4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.3, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	// Should not throw - NaN edge length allows any distal_length
	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 4);
	auto f = tree.find_node_by_name("F");
	REQUIRE(f.has_value());
}

TEST_CASE("insert_fully_resolved duplicate distal_lengths creates zero-length branch", "[insert][edge]") {
	// Two placements at same distal_length should create zero-length internal branch
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F1", .edge_id = 0, .distal_length = 0.5, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    {.fragment_id = "F2", .edge_id = 0, .distal_length = 0.5, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 5);

	// One of the internal nodes should have zero-length branch
	auto f1 = tree.find_node_by_name("F1").value();
	auto f2 = tree.find_node_by_name("F2").value();

	// F1 and F2 have different parents (each has its own internal node)
	// One internal should have branch_length 0
	uint32_t p1 = tree.parent(f1);
	uint32_t p2 = tree.parent(f2);
	REQUIRE(p1 != p2);

	// Check that one has zero-length branch (the one processed second)
	double len1 = tree.branch_length(p1);
	double len2 = tree.branch_length(p2);
	REQUIRE((Approx(len1).margin(1e-9) == 0.0 || Approx(len2).margin(1e-9) == 0.0));
}

TEST_CASE("insert_fully_resolved multiple placements on root edge", "[insert][root]") {
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):1.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F1", .edge_id = 4, .distal_length = 0.8, .pendant_length = 0.1, .like_weight_ratio = 1.0},
	    {.fragment_id = "F2", .edge_id = 4, .distal_length = 0.3, .pendant_length = 0.2, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 5);

	auto f1 = tree.find_node_by_name("F1").value();
	auto f2 = tree.find_node_by_name("F2").value();
	REQUIRE(tree.branch_length(f1) == Approx(0.1));
	REQUIRE(tree.branch_length(f2) == Approx(0.2));
}

TEST_CASE("insert_fully_resolved at distal_length = 0", "[insert][edge]") {
	// Placement exactly at the child end of the edge
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 0.0, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 4);

	// A should now have branch_length 0 (distal_length)
	auto a = tree.find_node_by_name("A").value();
	REQUIRE(tree.branch_length(a) == Approx(0.0));

	// F should be sibling of A
	auto f = tree.find_node_by_name("F").value();
	REQUIRE(tree.parent(f) == tree.parent(a));
}

TEST_CASE("insert_fully_resolved at distal_length = edge_length", "[insert][edge]") {
	// Placement exactly at the parent end of the edge
	auto tree = miint::NewickTree::parse("((A:1.0{0},B:2.0{1}):0.5{2},C:3.0{3}):0.0{4};");

	std::vector<miint::Placement> placements = {
	    {.fragment_id = "F", .edge_id = 0, .distal_length = 1.0, .pendant_length = 0.1, .like_weight_ratio = 1.0}};

	tree.insert_fully_resolved(placements);

	REQUIRE(tree.num_tips() == 4);

	// A should retain its full branch_length of 1.0
	auto a = tree.find_node_by_name("A").value();
	REQUIRE(tree.branch_length(a) == Approx(1.0));

	// The new internal node should have branch_length 0 (at parent position)
	auto f = tree.find_node_by_name("F").value();
	uint32_t internal = tree.parent(f);
	REQUIRE(tree.branch_length(internal) == Approx(0.0));
}
