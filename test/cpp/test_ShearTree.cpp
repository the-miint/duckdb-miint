#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <string>
#include <unordered_set>
#include "NewickTree.hpp"

using Catch::Matchers::ContainsSubstring;

namespace {

// Resolve a node index by name (fails the test if the name is absent).
uint32_t idx(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

// A binary tree with distinct branch lengths on every edge:
//
//                 root
//                /    \
//              H:5     I:8
//             /  \     /  \
//          G:3   C:4 D:6  E:7
//         /  \
//       A:1  B:2
const char *kTreeWithLengths = "(((A:1,B:2)G:3,C:4)H:5,(D:6,E:7)I:8)root;";

// Same topology, no branch lengths (a cladogram / topology-only tree).
const char *kCladogram = "(((A,B)G,C)H,(D,E)I)root;";

std::unordered_set<std::string> names_of(std::initializer_list<const char *> xs) {
	std::unordered_set<std::string> s;
	for (const char *x : xs) {
		s.insert(x);
	}
	return s;
}

} // namespace

// ============================================================================
// collapse=false : preserve every internal node on the kept paths
// ============================================================================

TEST_CASE("shear preserve keeps all ancestors of kept tips", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// Keep A, B, C -> the entire left subtree stays; I/D/E are dropped.
	auto out = tree.shear(names_of({"A", "B", "C"}), /*collapse=*/false, /*ignore_missing=*/false);

	// A,B,G,C,H,root == 6 nodes; 3 tips.
	REQUIRE(out.num_nodes() == 6);
	REQUIRE(out.num_tips() == 3);

	// Root ("root") is retained even though it now has a single child (H).
	// WHY: collapse=false must not remove unifurcations — callers ask for this
	// mode precisely when they need every internal node preserved.
	auto root = out.root();
	REQUIRE(out.name(root) == "root");
	REQUIRE(out.children(root).size() == 1);
	REQUIRE(out.find_node_by_name("I") == std::nullopt);
	REQUIRE(out.find_node_by_name("D") == std::nullopt);

	// Branch lengths on the surviving edges are unchanged (no summation).
	REQUIRE(out.branch_length(idx(out, "A")) == Catch::Approx(1.0));
	REQUIRE(out.branch_length(idx(out, "G")) == Catch::Approx(3.0));
	REQUIRE(out.branch_length(idx(out, "H")) == Catch::Approx(5.0));
	REQUIRE(out.branch_length(idx(out, "C")) == Catch::Approx(4.0));
}

TEST_CASE("shear preserve keeps the unifurcation chain above the LCA", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// Keep A, B: their LCA is G, but root->H->G is a single-child chain that
	// collapse=false must retain.
	auto out = tree.shear(names_of({"A", "B"}), /*collapse=*/false, /*ignore_missing=*/false);

	// A,B,G,H,root == 5 nodes.
	REQUIRE(out.num_nodes() == 5);
	REQUIRE(out.num_tips() == 2);
	REQUIRE(out.name(out.root()) == "root");
	REQUIRE(out.children(out.root()).size() == 1);    // root -> H only
	REQUIRE(out.children(idx(out, "H")).size() == 1); // H -> G only
}

// ============================================================================
// collapse=true : remove single-child ancestors, sum branch lengths
// ============================================================================

TEST_CASE("shear collapse removes single-descendant ancestors", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	auto out = tree.shear(names_of({"A", "B", "C"}), /*collapse=*/true, /*ignore_missing=*/false);

	// root had a single kept child (H) so it collapses away; H becomes the new
	// root. A,B,C,G,H == 5 nodes.
	REQUIRE(out.num_nodes() == 5);
	REQUIRE(out.num_tips() == 3);

	auto root = out.root();
	REQUIRE(out.name(root) == "H");
	// New root has no incoming edge -> branch length dropped to NaN.
	REQUIRE(std::isnan(out.branch_length(root)));

	// G (a real bifurcation) survives with its own branch length.
	REQUIRE(out.branch_length(idx(out, "G")) == Catch::Approx(3.0));
	REQUIRE(out.parent(idx(out, "G")) == root);
	REQUIRE(out.branch_length(idx(out, "C")) == Catch::Approx(4.0));
	REQUIRE(out.parent(idx(out, "C")) == root);
}

TEST_CASE("shear collapse sums branch lengths across the merged chain", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// Keep A and D: their LCA is the root. A's path A->G->H->root collapses to a
	// single edge of length 1+3+5=9; D's path D->I->root collapses to 6+8=14.
	// WHY: the whole point of collapse is that tip-to-root distances are
	// preserved after unifurcations are removed.
	auto out = tree.shear(names_of({"A", "D"}), /*collapse=*/true, /*ignore_missing=*/false);

	REQUIRE(out.num_nodes() == 3); // root, A, D
	REQUIRE(out.num_tips() == 2);
	REQUIRE(out.name(out.root()) == "root");
	REQUIRE(out.branch_length(idx(out, "A")) == Catch::Approx(9.0));
	REQUIRE(out.branch_length(idx(out, "D")) == Catch::Approx(14.0));
	REQUIRE(out.parent(idx(out, "A")) == out.root());
	REQUIRE(out.parent(idx(out, "D")) == out.root());
}

TEST_CASE("shear collapse of a topology-only tree keeps branch lengths NaN", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kCladogram);

	// Every edge is unspecified (NaN). Collapsing must not fabricate a 0.0
	// length — an all-NaN chain stays NaN so a cladogram round-trips as a
	// cladogram.
	auto out = tree.shear(names_of({"A", "D"}), /*collapse=*/true, /*ignore_missing=*/false);

	REQUIRE(out.num_nodes() == 3);
	REQUIRE(std::isnan(out.branch_length(idx(out, "A"))));
	REQUIRE(std::isnan(out.branch_length(idx(out, "D"))));
}

TEST_CASE("shear collapse drops nodes above the LCA", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// Keep A, B: LCA is G, so H and root (above G) disappear and G is the root.
	auto out = tree.shear(names_of({"A", "B"}), /*collapse=*/true, /*ignore_missing=*/false);

	REQUIRE(out.num_nodes() == 3); // G, A, B
	REQUIRE(out.name(out.root()) == "G");
	REQUIRE(std::isnan(out.branch_length(out.root())));
	REQUIRE(out.find_node_by_name("H") == std::nullopt);
	REQUIRE(out.find_node_by_name("root") == std::nullopt);
}

// ============================================================================
// Multifurcations, degenerate cases
// ============================================================================

TEST_CASE("shear collapse leaves genuine multifurcations intact", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:3)root;");

	SECTION("keeping all three tips is an identity") {
		auto out = tree.shear(names_of({"A", "B", "C"}), true, false);
		REQUIRE(out.num_nodes() == 4);
		REQUIRE(out.children(out.root()).size() == 3);
	}

	SECTION("keeping two of three still has a real bifurcation at the root") {
		auto out = tree.shear(names_of({"A", "B"}), true, false);
		REQUIRE(out.num_nodes() == 3);
		REQUIRE(out.name(out.root()) == "root");
		REQUIRE(out.children(out.root()).size() == 2);
	}
}

TEST_CASE("shear to a single tip yields a one-node tree", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	auto out = tree.shear(names_of({"C"}), /*collapse=*/true, /*ignore_missing=*/false);
	REQUIRE(out.num_nodes() == 1);
	REQUIRE(out.num_tips() == 1);
	REQUIRE(out.name(out.root()) == "C");
	REQUIRE(std::isnan(out.branch_length(out.root())));
}

TEST_CASE("shear keeping every tip preserves the full tree", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);
	auto out = tree.shear(names_of({"A", "B", "C", "D", "E"}), true, false);
	REQUIRE(out.num_nodes() == 9); // nothing collapses
	REQUIRE(out.num_tips() == 5);
}

// ============================================================================
// Unrooted / trifurcating root (the shape phylogeny_fasttree emits)
// ============================================================================

TEST_CASE("shear collapse handles a trifurcating (unrooted) root", "[NewickTree][shear]") {
	// root has THREE children F, G, I (an unrooted tree).
	auto tree = miint::NewickTree::parse("((A:1,B:2)F:3,(C:4,D:5)G:6,(E:7,X:8)I:9)root;");

	SECTION("keeping one tip from each of the 3 subtrees keeps the trifurcation") {
		// F,G,I each become unifurcations and collapse; root retains 3 kept
		// children, so the trifurcating root survives.
		auto out = tree.shear(names_of({"A", "C", "E"}), /*collapse=*/true, /*ignore_missing=*/false);
		REQUIRE(out.num_nodes() == 4); // root + A + C + E
		REQUIRE(out.name(out.root()) == "root");
		REQUIRE(out.children(out.root()).size() == 3);
		// Branch lengths summed across each collapsed subtree.
		REQUIRE(out.branch_length(idx(out, "A")) == Catch::Approx(4.0));  // 1+3
		REQUIRE(out.branch_length(idx(out, "C")) == Catch::Approx(10.0)); // 4+6
		REQUIRE(out.branch_length(idx(out, "E")) == Catch::Approx(16.0)); // 7+9
	}

	SECTION("keeping tips from only 2 subtrees leaves a bifurcating root") {
		auto out = tree.shear(names_of({"A", "C"}), /*collapse=*/true, /*ignore_missing=*/false);
		REQUIRE(out.num_nodes() == 3); // root + A + C
		REQUIRE(out.name(out.root()) == "root");
		REQUIRE(out.children(out.root()).size() == 2);
		REQUIRE(out.find_node_by_name("I") == std::nullopt);
	}

	SECTION("keeping tips from a single subtree collapses the root away") {
		// Only F's subtree is kept -> root has one kept child and collapses; F
		// (a real bifurcation) becomes the new root.
		auto out = tree.shear(names_of({"A", "B"}), /*collapse=*/true, /*ignore_missing=*/false);
		REQUIRE(out.num_nodes() == 3); // F + A + B
		REQUIRE(out.name(out.root()) == "F");
		REQUIRE(std::isnan(out.branch_length(out.root())));
		REQUIRE(out.children(out.root()).size() == 2);
		REQUIRE(out.find_node_by_name("root") == std::nullopt);
	}

	SECTION("keeping every tip is an identity") {
		auto out = tree.shear(names_of({"A", "B", "C", "D", "E", "X"}), /*collapse=*/true, /*ignore_missing=*/false);
		REQUIRE(out.num_nodes() == 10); // 6 tips + F,G,I,root
		REQUIRE(out.children(out.root()).size() == 3);
	}
}

// ============================================================================
// edge_id handling across a collapse
// ============================================================================

TEST_CASE("shear collapse keeps the surviving node's edge_id", "[NewickTree][shear][jplace]") {
	// jplace-style edge ids in {n} notation.
	auto tree = miint::NewickTree::parse("(((A:1{0},B:2{1})G:3{2},C:4{3})H:5{4},(D:6{5},E:7{6})I:8{7})root;");

	auto out = tree.shear(names_of({"A", "B"}), /*collapse=*/true, /*ignore_missing=*/false);
	// A and B survive with their own edge ids; the collapsed H/root ids vanish.
	REQUIRE(out.edge_id(idx(out, "A")) == std::optional<int64_t>(0));
	REQUIRE(out.edge_id(idx(out, "B")) == std::optional<int64_t>(1));
}

// ============================================================================
// Missing tips policy
// ============================================================================

TEST_CASE("shear throws when a requested tip is absent (strict)", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	REQUIRE_THROWS_WITH(tree.shear(names_of({"A", "ZZZ"}), true, /*ignore_missing=*/false), ContainsSubstring("ZZZ"));
}

TEST_CASE("shear ignore_missing skips absent tips", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// A matches, ZZZ does not -> with ignore_missing we still get A.
	auto out = tree.shear(names_of({"A", "ZZZ"}), true, /*ignore_missing=*/true);
	REQUIRE(out.num_nodes() == 1);
	REQUIRE(out.name(out.root()) == "A");
}

TEST_CASE("shear throws when no requested tip matches", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// Even with ignore_missing, an empty result is an error (can't shear to
	// nothing).
	REQUIRE_THROWS_WITH(tree.shear(names_of({"ZZZ", "QQQ"}), true, /*ignore_missing=*/true),
	                    ContainsSubstring("no tips"));
}

TEST_CASE("shear does not match internal-node labels", "[NewickTree][shear]") {
	auto tree = miint::NewickTree::parse(kTreeWithLengths);

	// "G" is an internal label, not a tip -> treated as missing.
	REQUIRE_THROWS_WITH(tree.shear(names_of({"G"}), true, /*ignore_missing=*/false), ContainsSubstring("G"));
}
