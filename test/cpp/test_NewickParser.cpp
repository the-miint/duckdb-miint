#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <string>
#include "NewickTree.hpp"

using Catch::Matchers::ContainsSubstring;

// ============================================================================
// Parsing: Basic structure
// ============================================================================

TEST_CASE("NewickTree parse single leaf", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("A;");

	REQUIRE(tree.num_nodes() == 1);
	REQUIRE(tree.num_tips() == 1);
	REQUIRE(tree.name(tree.root()) == "A");
	REQUIRE(tree.is_tip(tree.root()));
}

TEST_CASE("NewickTree parse simple two-leaf tree", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A,B);");

	REQUIRE(tree.num_nodes() == 3);
	REQUIRE(tree.num_tips() == 2);

	auto root = tree.root();
	REQUIRE(tree.name(root) == "");
	REQUIRE_FALSE(tree.is_tip(root));

	auto children = tree.children(root);
	REQUIRE(children.size() == 2);

	// Children should be A and B (order preserved)
	REQUIRE(tree.name(children[0]) == "A");
	REQUIRE(tree.name(children[1]) == "B");
	REQUIRE(tree.is_tip(children[0]));
	REQUIRE(tree.is_tip(children[1]));

	// Parent relationships
	REQUIRE(tree.parent(children[0]) == root);
	REQUIRE(tree.parent(children[1]) == root);
}

TEST_CASE("NewickTree parse nested tree", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("((A,B),(C,D));");

	REQUIRE(tree.num_nodes() == 7);
	REQUIRE(tree.num_tips() == 4);

	auto root = tree.root();
	auto root_children = tree.children(root);
	REQUIRE(root_children.size() == 2);

	// First subtree (A,B)
	auto left = root_children[0];
	REQUIRE_FALSE(tree.is_tip(left));
	auto left_children = tree.children(left);
	REQUIRE(left_children.size() == 2);
	REQUIRE(tree.name(left_children[0]) == "A");
	REQUIRE(tree.name(left_children[1]) == "B");

	// Second subtree (C,D)
	auto right = root_children[1];
	REQUIRE_FALSE(tree.is_tip(right));
	auto right_children = tree.children(right);
	REQUIRE(right_children.size() == 2);
	REQUIRE(tree.name(right_children[0]) == "C");
	REQUIRE(tree.name(right_children[1]) == "D");
}

TEST_CASE("NewickTree parse multifurcating tree", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A,B,C,D);");

	REQUIRE(tree.num_nodes() == 5);
	REQUIRE(tree.num_tips() == 4);

	auto children = tree.children(tree.root());
	REQUIRE(children.size() == 4);
	REQUIRE(tree.name(children[0]) == "A");
	REQUIRE(tree.name(children[1]) == "B");
	REQUIRE(tree.name(children[2]) == "C");
	REQUIRE(tree.name(children[3]) == "D");
}

// ============================================================================
// Parsing: Branch lengths
// ============================================================================

TEST_CASE("NewickTree parse branch lengths", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A:0.1,B:0.2):0.3;");

	auto root = tree.root();
	REQUIRE(tree.branch_length(root) == Catch::Approx(0.3));

	auto children = tree.children(root);
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.1));
	REQUIRE(tree.branch_length(children[1]) == Catch::Approx(0.2));
}

TEST_CASE("NewickTree parse missing branch lengths are NaN", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A,B);");

	auto root = tree.root();
	REQUIRE(std::isnan(tree.branch_length(root)));

	auto children = tree.children(root);
	REQUIRE(std::isnan(tree.branch_length(children[0])));
	REQUIRE(std::isnan(tree.branch_length(children[1])));
}

TEST_CASE("NewickTree parse zero branch length", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A:0.0,B:0):0;");

	auto root = tree.root();
	REQUIRE(tree.branch_length(root) == Catch::Approx(0.0));

	auto children = tree.children(root);
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.0));
	REQUIRE(tree.branch_length(children[1]) == Catch::Approx(0.0));
}

TEST_CASE("NewickTree parse scientific notation branch length", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A:1.5e-3,B:2E+2);");

	auto children = tree.children(tree.root());
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.0015));
	REQUIRE(tree.branch_length(children[1]) == Catch::Approx(200.0));
}

// ============================================================================
// Parsing: Internal node labels
// ============================================================================

TEST_CASE("NewickTree parse internal node labels", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(A,B)Root;");

	REQUIRE(tree.name(tree.root()) == "Root");
}

TEST_CASE("NewickTree parse internal labels with branch lengths", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("((A,B)AB:0.5,(C,D)CD:0.6)Root:0.1;");

	REQUIRE(tree.name(tree.root()) == "Root");
	REQUIRE(tree.branch_length(tree.root()) == Catch::Approx(0.1));

	auto children = tree.children(tree.root());
	REQUIRE(tree.name(children[0]) == "AB");
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.5));
	REQUIRE(tree.name(children[1]) == "CD");
	REQUIRE(tree.branch_length(children[1]) == Catch::Approx(0.6));
}

// ============================================================================
// Parsing: Edge identifiers (jplace format)
// ============================================================================

TEST_CASE("NewickTree parse edge identifiers", "[NewickTree][parse][jplace]") {
	auto tree = miint::NewickTree::parse("(A:0.1{0},B:0.2{1}):0.3{2};");

	auto root = tree.root();
	REQUIRE(tree.edge_id(root).has_value());
	REQUIRE(tree.edge_id(root).value() == 2);

	auto children = tree.children(root);
	REQUIRE(tree.edge_id(children[0]).has_value());
	REQUIRE(tree.edge_id(children[0]).value() == 0);
	REQUIRE(tree.edge_id(children[1]).has_value());
	REQUIRE(tree.edge_id(children[1]).value() == 1);
}

TEST_CASE("NewickTree parse edge identifiers without branch lengths", "[NewickTree][parse][jplace]") {
	auto tree = miint::NewickTree::parse("(A{0},B{1}){2};");

	auto root = tree.root();
	REQUIRE(tree.edge_id(root).value() == 2);
	REQUIRE(std::isnan(tree.branch_length(root)));

	auto children = tree.children(root);
	REQUIRE(tree.edge_id(children[0]).value() == 0);
	REQUIRE(tree.edge_id(children[1]).value() == 1);
}

TEST_CASE("NewickTree parse mixed edge identifiers", "[NewickTree][parse][jplace]") {
	// Some nodes have edge IDs, some don't
	auto tree = miint::NewickTree::parse("(A:0.1{0},B:0.2):0.3{2};");

	auto root = tree.root();
	REQUIRE(tree.edge_id(root).value() == 2);

	auto children = tree.children(root);
	REQUIRE(tree.edge_id(children[0]).value() == 0);
	REQUIRE_FALSE(tree.edge_id(children[1]).has_value());
}

TEST_CASE("NewickTree parse complex jplace tree", "[NewickTree][parse][jplace]") {
	// Real-world-like jplace tree
	auto tree = miint::NewickTree::parse("((A:0.1{0},B:0.2{1})AB:0.3{2},(C:0.4{3},D:0.5{4})CD:0.6{5})root:0.0{6};");

	REQUIRE(tree.num_nodes() == 7);
	REQUIRE(tree.num_tips() == 4);

	auto root = tree.root();
	REQUIRE(tree.name(root) == "root");
	REQUIRE(tree.edge_id(root).value() == 6);

	// Verify all edge IDs are present and unique
	std::set<int64_t> edge_ids;
	for (size_t i = 0; i < tree.num_nodes(); ++i) {
		auto eid = tree.edge_id(i);
		REQUIRE(eid.has_value());
		REQUIRE(edge_ids.find(eid.value()) == edge_ids.end());
		edge_ids.insert(eid.value());
	}
	REQUIRE(edge_ids.size() == 7);
}

// ============================================================================
// Parsing: Quoted labels
// ============================================================================

TEST_CASE("NewickTree parse quoted labels", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("('Species A','Species B');");

	auto children = tree.children(tree.root());
	REQUIRE(tree.name(children[0]) == "Species A");
	REQUIRE(tree.name(children[1]) == "Species B");
}

TEST_CASE("NewickTree parse quoted labels with special characters", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("('A:0.1',\"B,C\");");

	auto children = tree.children(tree.root());
	REQUIRE(tree.name(children[0]) == "A:0.1");
	REQUIRE(tree.name(children[1]) == "B,C");
}

TEST_CASE("NewickTree parse quoted labels with semicolons", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("('foo; bar','baz;qux');");

	auto children = tree.children(tree.root());
	REQUIRE(children.size() == 2);
	REQUIRE(tree.name(children[0]) == "foo; bar");
	REQUIRE(tree.name(children[1]) == "baz;qux");
}

TEST_CASE("NewickTree parse quoted labels with escaped quotes", "[NewickTree][parse]") {
	// Note: ('label'); creates a tree with internal root and one child
	auto tree = miint::NewickTree::parse("('It''s a test');");

	auto children = tree.children(tree.root());
	REQUIRE(children.size() == 1);
	REQUIRE(tree.name(children[0]) == "It's a test");
}

// ============================================================================
// Parsing: Whitespace handling
// ============================================================================

TEST_CASE("NewickTree parse with whitespace", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("( A , B ) ;");

	REQUIRE(tree.num_nodes() == 3);
	auto children = tree.children(tree.root());
	REQUIRE(tree.name(children[0]) == "A");
	REQUIRE(tree.name(children[1]) == "B");
}

TEST_CASE("NewickTree parse with newlines", "[NewickTree][parse]") {
	auto tree = miint::NewickTree::parse("(\n  A:0.1,\n  B:0.2\n);");

	REQUIRE(tree.num_nodes() == 3);
}

// ============================================================================
// Parsing: Error handling (strict mode)
// ============================================================================

TEST_CASE("NewickTree parse error: empty string", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse(""), ContainsSubstring("empty"));
}

TEST_CASE("NewickTree parse error: missing semicolon", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A,B)"), ContainsSubstring("semicolon"));
}

TEST_CASE("NewickTree parse error: unmatched opening parenthesis", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A,B;"), ContainsSubstring("parenthes"));
}

TEST_CASE("NewickTree parse error: unexpected content after node", "[NewickTree][parse][error]") {
	// "A,B);" has extra content after the first node - commas only valid inside parens
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("A,B);"), ContainsSubstring("semicolon"));
}

TEST_CASE("NewickTree parse error: unclosed brace", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A:0.1{0,B:0.2{1});"), ContainsSubstring("brace"));
}

TEST_CASE("NewickTree parse error: invalid branch length", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A:abc,B);"), ContainsSubstring("branch length"));
}

TEST_CASE("NewickTree parse error: invalid edge identifier", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A{abc},B);"), ContainsSubstring("edge"));
}

TEST_CASE("NewickTree parse error: unclosed quote", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("('A,B);"), ContainsSubstring("quote"));
}

TEST_CASE("NewickTree parse error: unclosed comment", "[NewickTree][parse][error]") {
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A[unclosed comment,B);"), ContainsSubstring("Unclosed comment"));
}

TEST_CASE("NewickTree parse error: malformed branch length", "[NewickTree][parse][error]") {
	// Multiple decimal points
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A:1.2.3,B);"), ContainsSubstring("branch length"));
	// Multiple minus signs
	REQUIRE_THROWS_WITH(miint::NewickTree::parse("(A:--5,B);"), ContainsSubstring("branch length"));
}

TEST_CASE("NewickTree modification error: invalid index", "[NewickTree][modify][error]") {
	auto tree = miint::NewickTree::parse("(A,B);");

	SECTION("set_parent with invalid child") {
		REQUIRE_THROWS_AS(tree.set_parent(999, 0), std::out_of_range);
	}

	SECTION("set_parent with invalid parent") {
		REQUIRE_THROWS_AS(tree.set_parent(0, 999), std::out_of_range);
	}

	SECTION("set_parent with self-reference") {
		REQUIRE_THROWS_AS(tree.set_parent(0, 0), std::invalid_argument);
	}

	SECTION("remove_child with invalid parent") {
		REQUIRE_THROWS_AS(tree.remove_child(999, 0), std::out_of_range);
	}
}

// ============================================================================
// Serialization: to_newick()
// ============================================================================

TEST_CASE("NewickTree to_newick single leaf", "[NewickTree][serialize]") {
	auto tree = miint::NewickTree::parse("A;");
	REQUIRE(tree.to_newick() == "A;");
}

TEST_CASE("NewickTree to_newick simple tree", "[NewickTree][serialize]") {
	auto tree = miint::NewickTree::parse("(A,B);");
	REQUIRE(tree.to_newick() == "(A,B);");
}

TEST_CASE("NewickTree to_newick with branch lengths", "[NewickTree][serialize]") {
	auto tree = miint::NewickTree::parse("(A:0.1,B:0.2):0.3;");
	auto newick = tree.to_newick();
	// Re-parse and verify (avoids floating point formatting issues)
	auto tree2 = miint::NewickTree::parse(newick);
	REQUIRE(tree2.branch_length(tree2.root()) == Catch::Approx(0.3));
}

TEST_CASE("NewickTree to_newick with edge identifiers", "[NewickTree][serialize]") {
	auto tree = miint::NewickTree::parse("(A:0.1{0},B:0.2{1}):0.3{2};");
	auto newick = tree.to_newick();
	auto tree2 = miint::NewickTree::parse(newick);

	REQUIRE(tree2.edge_id(tree2.root()).value() == 2);
	auto children = tree2.children(tree2.root());
	REQUIRE(tree2.edge_id(children[0]).value() == 0);
	REQUIRE(tree2.edge_id(children[1]).value() == 1);
}

TEST_CASE("NewickTree to_newick roundtrip complex tree", "[NewickTree][serialize]") {
	std::string original = "((A:0.1{0},B:0.2{1})AB:0.3{2},(C:0.4{3},D:0.5{4})CD:0.6{5})root:0.0{6};";
	auto tree = miint::NewickTree::parse(original);
	auto newick = tree.to_newick();
	auto tree2 = miint::NewickTree::parse(newick);

	REQUIRE(tree2.num_nodes() == tree.num_nodes());
	REQUIRE(tree2.num_tips() == tree.num_tips());

	// Verify structure is identical
	for (size_t i = 0; i < tree.num_nodes(); ++i) {
		REQUIRE(tree2.name(i) == tree.name(i));
		if (std::isnan(tree.branch_length(i))) {
			REQUIRE(std::isnan(tree2.branch_length(i)));
		} else {
			REQUIRE(tree2.branch_length(i) == Catch::Approx(tree.branch_length(i)));
		}
		REQUIRE(tree2.edge_id(i) == tree.edge_id(i));
	}
}

// ============================================================================
// Traversal
// ============================================================================

TEST_CASE("NewickTree postorder traversal", "[NewickTree][traversal]") {
	auto tree = miint::NewickTree::parse("((A,B),(C,D));");
	auto postorder = tree.postorder();

	REQUIRE(postorder.size() == 7);

	// Root should be last in postorder
	REQUIRE(postorder.back() == tree.root());

	// All children should appear before their parents
	for (size_t i = 0; i < postorder.size(); ++i) {
		auto node = postorder[i];
		for (auto child : tree.children(node)) {
			// Find child's position
			auto child_pos = std::find(postorder.begin(), postorder.end(), child) - postorder.begin();
			REQUIRE(child_pos < static_cast<ptrdiff_t>(i));
		}
	}
}

TEST_CASE("NewickTree preorder traversal", "[NewickTree][traversal]") {
	auto tree = miint::NewickTree::parse("((A,B),(C,D));");
	auto preorder = tree.preorder();

	REQUIRE(preorder.size() == 7);

	// Root should be first in preorder
	REQUIRE(preorder.front() == tree.root());

	// All parents should appear before their children
	for (size_t i = 0; i < preorder.size(); ++i) {
		auto node = preorder[i];
		for (auto child : tree.children(node)) {
			auto child_pos = std::find(preorder.begin(), preorder.end(), child) - preorder.begin();
			REQUIRE(child_pos > static_cast<ptrdiff_t>(i));
		}
	}
}

// ============================================================================
// Tree modification (for insert_fully_resolved)
// ============================================================================

TEST_CASE("NewickTree add_node", "[NewickTree][modify]") {
	auto tree = miint::NewickTree::parse("(A,B);");
	size_t initial_count = tree.num_nodes();

	auto new_node = tree.add_node("C", 0.5, std::nullopt);

	REQUIRE(tree.num_nodes() == initial_count + 1);
	REQUIRE(tree.name(new_node) == "C");
	REQUIRE(tree.branch_length(new_node) == Catch::Approx(0.5));
	REQUIRE_FALSE(tree.edge_id(new_node).has_value());
}

TEST_CASE("NewickTree set_parent and remove_child", "[NewickTree][modify]") {
	auto tree = miint::NewickTree::parse("(A,B);");
	auto root = tree.root();
	auto children = tree.children(root);
	auto a = children[0];

	// Add new internal node
	auto new_internal = tree.add_node("", 0.1, std::nullopt);

	// Remove A from root
	tree.remove_child(root, a);
	REQUIRE(tree.children(root).size() == 1);

	// Set new_internal as child of root
	tree.set_parent(new_internal, root);
	REQUIRE(tree.children(root).size() == 2);
	REQUIRE(tree.parent(new_internal) == root);

	// Set A as child of new_internal
	tree.set_parent(a, new_internal);
	REQUIRE(tree.parent(a) == new_internal);
	REQUIRE(tree.children(new_internal).size() == 1);
}

TEST_CASE("NewickTree set_branch_length", "[NewickTree][modify]") {
	auto tree = miint::NewickTree::parse("(A:0.1,B:0.2);");
	auto children = tree.children(tree.root());

	tree.set_branch_length(children[0], 0.5);
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.5));
}

TEST_CASE("NewickTree set_edge_id", "[NewickTree][modify]") {
	auto tree = miint::NewickTree::parse("(A,B);");
	auto children = tree.children(tree.root());

	REQUIRE_FALSE(tree.edge_id(children[0]).has_value());

	tree.set_edge_id(children[0], 42);
	REQUIRE(tree.edge_id(children[0]).value() == 42);

	tree.set_edge_id(children[0], std::nullopt);
	REQUIRE_FALSE(tree.edge_id(children[0]).has_value());
}

// ============================================================================
// Scalability tests
// ============================================================================

TEST_CASE("NewickTree parse large balanced tree", "[NewickTree][scale]") {
	// Build a tree string with 1024 tips (10 levels deep)
	std::function<std::string(int, int)> build_subtree = [&](int depth, int id) -> std::string {
		if (depth == 0) {
			return "t" + std::to_string(id);
		}
		int left_id = id * 2;
		int right_id = id * 2 + 1;
		return "(" + build_subtree(depth - 1, left_id) + "," + build_subtree(depth - 1, right_id) + ")";
	};

	std::string newick = build_subtree(10, 1) + ";";

	auto tree = miint::NewickTree::parse(newick);

	REQUIRE(tree.num_tips() == 1024);
	REQUIRE(tree.num_nodes() == 2047); // 2^11 - 1
}

TEST_CASE("NewickTree parse caterpillar tree", "[NewickTree][scale]") {
	// Caterpillar (maximally unbalanced) tree with 1000 tips
	// ((((A,B),C),D),E)...
	std::string newick = "t0";
	for (int i = 1; i < 1000; ++i) {
		newick = "(" + newick + ",t" + std::to_string(i) + ")";
	}
	newick += ";";

	auto tree = miint::NewickTree::parse(newick);

	REQUIRE(tree.num_tips() == 1000);
	REQUIRE(tree.num_nodes() == 1999); // 2n - 1 for caterpillar
}

TEST_CASE("NewickTree serialize deep caterpillar tree", "[NewickTree][serialize][scale]") {
	// Test that iterative serialization handles deep trees without stack overflow
	// Build a caterpillar tree with 5000 tips (would overflow stack if recursive)
	std::string newick = "t0:0.1";
	for (int i = 1; i < 5000; ++i) {
		newick = "(" + newick + ",t" + std::to_string(i) + ":0.1)";
	}
	newick += ";";

	auto tree = miint::NewickTree::parse(newick);
	REQUIRE(tree.num_tips() == 5000);

	// Serialize and verify it can be re-parsed
	std::string serialized = tree.to_newick();
	auto tree2 = miint::NewickTree::parse(serialized);

	REQUIRE(tree2.num_nodes() == tree.num_nodes());
	REQUIRE(tree2.num_tips() == tree.num_tips());
}

TEST_CASE("NewickTree find_node_by_name", "[NewickTree][query]") {
	auto tree = miint::NewickTree::parse("((A,B)AB,(C,D)CD)root;");

	auto a = tree.find_node_by_name("A");
	REQUIRE(a.has_value());
	REQUIRE(tree.name(a.value()) == "A");

	auto ab = tree.find_node_by_name("AB");
	REQUIRE(ab.has_value());
	REQUIRE_FALSE(tree.is_tip(ab.value()));

	auto missing = tree.find_node_by_name("X");
	REQUIRE_FALSE(missing.has_value());
}

TEST_CASE("NewickTree find_node_by_edge_id", "[NewickTree][query][jplace]") {
	auto tree = miint::NewickTree::parse("((A{0},B{1}){2},(C{3},D{4}){5}){6};");

	for (int64_t i = 0; i <= 6; ++i) {
		auto node = tree.find_node_by_edge_id(i);
		REQUIRE(node.has_value());
		REQUIRE(tree.edge_id(node.value()).value() == i);
	}

	auto missing = tree.find_node_by_edge_id(99);
	REQUIRE_FALSE(missing.has_value());
}

// ============================================================================
// Tests inspired by scikit-bio (BSD-3-Clause License)
// https://github.com/scikit-bio/scikit-bio
// scikit-bio/skbio/io/format/tests/test_newick.py
// ============================================================================

TEST_CASE("NewickTree parse blank tree structures", "[NewickTree][parse][skbio]") {
	// Trees with unnamed nodes
	SECTION("nested empty nodes") {
		auto tree = miint::NewickTree::parse("(,,(,));");
		REQUIRE(tree.num_nodes() == 6);
		REQUIRE(tree.num_tips() == 4);
		// All nodes should have empty names
		for (size_t i = 0; i < tree.num_nodes(); ++i) {
			REQUIRE(tree.name(i).empty());
		}
	}

	SECTION("variant ordering") {
		auto tree = miint::NewickTree::parse("(,(,),);");
		REQUIRE(tree.num_nodes() == 6);
		REQUIRE(tree.num_tips() == 4);
	}

	SECTION("different nesting") {
		auto tree = miint::NewickTree::parse("((,),,);");
		REQUIRE(tree.num_nodes() == 6);
		REQUIRE(tree.num_tips() == 4);
	}

	SECTION("balanced binary blank") {
		auto tree = miint::NewickTree::parse("((,),(,));");
		REQUIRE(tree.num_nodes() == 7);
		REQUIRE(tree.num_tips() == 4);
	}
}

TEST_CASE("NewickTree parse single node trees", "[NewickTree][parse][skbio]") {
	SECTION("empty root only") {
		auto tree = miint::NewickTree::parse(";");
		REQUIRE(tree.num_nodes() == 1);
		REQUIRE(tree.num_tips() == 1);
		REQUIRE(tree.name(tree.root()).empty());
		REQUIRE(std::isnan(tree.branch_length(tree.root())));
	}

	SECTION("named only") {
		auto tree = miint::NewickTree::parse("athing;");
		REQUIRE(tree.num_nodes() == 1);
		REQUIRE(tree.name(tree.root()) == "athing");
		REQUIRE(std::isnan(tree.branch_length(tree.root())));
	}

	SECTION("distance only") {
		auto tree = miint::NewickTree::parse(":200.0;");
		REQUIRE(tree.num_nodes() == 1);
		REQUIRE(tree.name(tree.root()).empty());
		REQUIRE(tree.branch_length(tree.root()) == Catch::Approx(200.0));
	}

	SECTION("quoted name with special chars and distance") {
		auto tree = miint::NewickTree::parse("'[a]':200.0;");
		REQUIRE(tree.num_nodes() == 1);
		REQUIRE(tree.name(tree.root()) == "[a]");
		REQUIRE(tree.branch_length(tree.root()) == Catch::Approx(200.0));
	}
}

TEST_CASE("NewickTree parse deeply nested empty", "[NewickTree][parse][skbio]") {
	auto tree = miint::NewickTree::parse("((((()))));");
	// Structure: ((((()))));
	// Each level of parens creates a new internal node wrapping the inner
	// () = 1, (()) = 2, ((())) = 3, (((()))) = 4, ((((())))  = 5, (((((())))) = 6
	REQUIRE(tree.num_nodes() == 6);
	REQUIRE(tree.num_tips() == 1);
}

TEST_CASE("NewickTree parse comments", "[NewickTree][parse][skbio]") {
	SECTION("comment before tree") {
		auto tree = miint::NewickTree::parse("[comment](A,B);");
		REQUIRE(tree.num_tips() == 2);
		auto children = tree.children(tree.root());
		REQUIRE(tree.name(children[0]) == "A");
		REQUIRE(tree.name(children[1]) == "B");
	}

	SECTION("comment within tree") {
		// Structure: ((,COMMENT),) -> root has 2 children: inner (with 2 empty) and empty
		auto tree = miint::NewickTree::parse("((,[this is a comment]),);");
		REQUIRE(tree.num_nodes() == 5);
	}

	SECTION("complex nested comment") {
		// Same structure with nested brackets in comment
		auto tree = miint::NewickTree::parse("((,[i_can_do_this[0] or escape unmatched]),);");
		REQUIRE(tree.num_nodes() == 5);
	}

	SECTION("comment between tokens") {
		auto tree = miint::NewickTree::parse("(A[comment]:0.1,B);");
		auto children = tree.children(tree.root());
		REQUIRE(tree.name(children[0]) == "A");
		REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.1));
	}

	SECTION("comment after closing paren") {
		auto tree = miint::NewickTree::parse("(A,B)[root comment];");
		REQUIRE(tree.num_tips() == 2);
	}
}

TEST_CASE("NewickTree parse linked list chain", "[NewickTree][parse][skbio]") {
	SECTION("with distances") {
		auto tree = miint::NewickTree::parse("((((:0.4):0.3):0.2):0.1):0.0;");
		REQUIRE(tree.num_nodes() == 5);
		REQUIRE(tree.num_tips() == 1);

		// Verify chain structure
		auto root = tree.root();
		REQUIRE(tree.branch_length(root) == Catch::Approx(0.0));
		REQUIRE(tree.children(root).size() == 1);

		auto n1 = tree.children(root)[0];
		REQUIRE(tree.branch_length(n1) == Catch::Approx(0.1));

		auto n2 = tree.children(n1)[0];
		REQUIRE(tree.branch_length(n2) == Catch::Approx(0.2));
	}

	SECTION("with names") {
		auto tree = miint::NewickTree::parse("((((aaa)bbb)ccc)ddd)eee;");
		REQUIRE(tree.num_nodes() == 5);
		REQUIRE(tree.name(tree.root()) == "eee");

		auto node = tree.children(tree.root())[0];
		REQUIRE(tree.name(node) == "ddd");
	}
}

TEST_CASE("NewickTree parse distance only nodes", "[NewickTree][parse][skbio]") {
	auto tree = miint::NewickTree::parse("(:0.1,:0.2,(:0.3,:0.4):0.5);");

	REQUIRE(tree.num_nodes() == 6);
	REQUIRE(tree.num_tips() == 4);

	auto root = tree.root();
	REQUIRE(std::isnan(tree.branch_length(root)));

	auto children = tree.children(root);
	REQUIRE(children.size() == 3);
	REQUIRE(tree.branch_length(children[0]) == Catch::Approx(0.1));
	REQUIRE(tree.branch_length(children[1]) == Catch::Approx(0.2));
	REQUIRE(tree.branch_length(children[2]) == Catch::Approx(0.5));

	// Check nested children
	auto nested = tree.children(children[2]);
	REQUIRE(tree.branch_length(nested[0]) == Catch::Approx(0.3));
	REQUIRE(tree.branch_length(nested[1]) == Catch::Approx(0.4));
}

TEST_CASE("NewickTree parse root with distance", "[NewickTree][parse][skbio]") {
	auto tree = miint::NewickTree::parse("(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;");

	REQUIRE(tree.branch_length(tree.root()) == Catch::Approx(0.0));
}

TEST_CASE("NewickTree parse mixed quoted and unquoted names", "[NewickTree][parse][skbio]") {
	auto tree = miint::NewickTree::parse("('a_',b,(c,d));");

	auto root_children = tree.children(tree.root());
	REQUIRE(root_children.size() == 3);
	REQUIRE(tree.name(root_children[0]) == "a_");
	REQUIRE(tree.name(root_children[1]) == "b");

	auto nested = tree.children(root_children[2]);
	REQUIRE(tree.name(nested[0]) == "c");
	REQUIRE(tree.name(nested[1]) == "d");
}

TEST_CASE("NewickTree parse names with colons", "[NewickTree][parse][skbio]") {
	// Colons in names must be quoted
	auto tree = miint::NewickTree::parse("((c:3.0,d:4.0)a:1.0,(e:5.0,'f:f''f':6.0)b:2.0)g:0.0;");

	REQUIRE(tree.name(tree.root()) == "g");
	REQUIRE(tree.branch_length(tree.root()) == Catch::Approx(0.0));

	auto root_children = tree.children(tree.root());
	REQUIRE(tree.name(root_children[0]) == "a");
	REQUIRE(tree.name(root_children[1]) == "b");

	// Find the node with colon in name
	auto node = tree.find_node_by_name("f:f'f");
	REQUIRE(node.has_value());
	REQUIRE(tree.branch_length(node.value()) == Catch::Approx(6.0));
}

TEST_CASE("NewickTree parse special characters in quoted names", "[NewickTree][parse][skbio]") {
	SECTION("brackets in name") {
		auto tree = miint::NewickTree::parse("(a,b,(c,'[whaaat!'']')e)f;");
		auto node = tree.find_node_by_name("[whaaat!']");
		REQUIRE(node.has_value());
	}

	SECTION("parentheses in name") {
		auto tree = miint::NewickTree::parse("('(nested)',normal);");
		auto children = tree.children(tree.root());
		REQUIRE(tree.name(children[0]) == "(nested)");
	}
}

TEST_CASE("NewickTree parse underscores", "[NewickTree][parse][skbio]") {
	// Our parser preserves underscores (no conversion to spaces)
	auto tree = miint::NewickTree::parse("(species_one,species_two);");
	auto children = tree.children(tree.root());
	REQUIRE(tree.name(children[0]) == "species_one");
	REQUIRE(tree.name(children[1]) == "species_two");
}

TEST_CASE("NewickTree parse with whitespace variations", "[NewickTree][parse][skbio]") {
	SECTION("spaces everywhere") {
		auto tree = miint::NewickTree::parse("   (   (   ,   )   ,   ,   )   ;   ");
		REQUIRE(tree.num_nodes() == 6);
	}

	SECTION("tabs and newlines") {
		auto tree = miint::NewickTree::parse("(\n\tA,\n\tB\n);");
		auto children = tree.children(tree.root());
		REQUIRE(tree.name(children[0]) == "A");
		REQUIRE(tree.name(children[1]) == "B");
	}
}

TEST_CASE("NewickTree roundtrip with complex structure", "[NewickTree][serialize][skbio]") {
	std::vector<std::string> test_cases = {
	    "(,,(,));",
	    "((,),(,));",
	    "(A:0.1,B:0.2,C:0.3);",
	    "((A,B)AB,(C,D)CD)root;",
	    "(:0.1,:0.2,(:0.3,:0.4):0.5):0.0;",
	};

	for (const auto &original : test_cases) {
		auto tree = miint::NewickTree::parse(original);
		auto serialized = tree.to_newick();
		auto reparsed = miint::NewickTree::parse(serialized);

		REQUIRE(reparsed.num_nodes() == tree.num_nodes());
		REQUIRE(reparsed.num_tips() == tree.num_tips());
	}
}

// ============================================================================
// Build from node data
// ============================================================================

TEST_CASE("NewickTree build simple tree", "[NewickTree][build]") {
	// Build: ((A,B),C)
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "", .branch_length = std::numeric_limits<double>::quiet_NaN(), .edge_id = std::nullopt},  // root
	    {.node_id = 1, .parent_id = 0, .name = "", .branch_length = 0.5, .edge_id = std::nullopt},  // internal
	    {.node_id = 2, .parent_id = 1, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},
	    {.node_id = 3, .parent_id = 1, .name = "B", .branch_length = 0.2, .edge_id = std::nullopt},
	    {.node_id = 4, .parent_id = 0, .name = "C", .branch_length = 0.3, .edge_id = std::nullopt},
	};

	auto tree = miint::NewickTree::build(nodes);

	REQUIRE(tree.num_nodes() == 5);
	REQUIRE(tree.num_tips() == 3);

	auto tip_names = tree.tip_names();
	std::set<std::string> names(tip_names.begin(), tip_names.end());
	REQUIRE(names == std::set<std::string>{"A", "B", "C"});
}

TEST_CASE("NewickTree build with edge IDs", "[NewickTree][build]") {
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "root", .branch_length = 0.0, .edge_id = 4},
	    {.node_id = 1, .parent_id = 0, .name = "A", .branch_length = 0.1, .edge_id = 0},
	    {.node_id = 2, .parent_id = 0, .name = "B", .branch_length = 0.2, .edge_id = 1},
	};

	auto tree = miint::NewickTree::build(nodes);

	REQUIRE(tree.num_nodes() == 3);

	auto a_node = tree.find_node_by_name("A");
	REQUIRE(a_node.has_value());
	REQUIRE(tree.edge_id(a_node.value()) == 0);

	auto root_node = tree.find_node_by_name("root");
	REQUIRE(root_node.has_value());
	REQUIRE(tree.edge_id(root_node.value()) == 4);
}

TEST_CASE("NewickTree build preserves structure for serialization", "[NewickTree][build]") {
	// Build a tree and serialize it
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 100, .parent_id = std::nullopt, .name = "", .branch_length = std::numeric_limits<double>::quiet_NaN(), .edge_id = std::nullopt},
	    {.node_id = 101, .parent_id = 100, .name = "A", .branch_length = 1.5, .edge_id = std::nullopt},
	    {.node_id = 102, .parent_id = 100, .name = "B", .branch_length = 2.5, .edge_id = std::nullopt},
	};

	auto tree = miint::NewickTree::build(nodes);
	auto newick = tree.to_newick();

	// Parse the serialized tree
	auto reparsed = miint::NewickTree::parse(newick);

	REQUIRE(reparsed.num_nodes() == 3);
	REQUIRE(reparsed.num_tips() == 2);

	auto tip_names = reparsed.tip_names();
	std::set<std::string> names(tip_names.begin(), tip_names.end());
	REQUIRE(names == std::set<std::string>{"A", "B"});
}

TEST_CASE("NewickTree build empty list throws", "[NewickTree][build][error]") {
	std::vector<miint::NodeInput> nodes;
	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("empty"));
}

TEST_CASE("NewickTree build no root throws", "[NewickTree][build][error]") {
	// All nodes have parents - no root
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = 1, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},
	    {.node_id = 1, .parent_id = 0, .name = "B", .branch_length = 0.2, .edge_id = std::nullopt},
	};

	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("root"));
}

TEST_CASE("NewickTree build multiple roots throws", "[NewickTree][build][error]") {
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},
	    {.node_id = 1, .parent_id = std::nullopt, .name = "B", .branch_length = 0.2, .edge_id = std::nullopt},
	};

	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("Multiple roots"));
}

TEST_CASE("NewickTree build invalid parent throws", "[NewickTree][build][error]") {
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "root", .branch_length = 0.0, .edge_id = std::nullopt},
	    {.node_id = 1, .parent_id = 999, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},  // Invalid parent
	};

	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("non-existent parent"));
}

TEST_CASE("NewickTree build duplicate node_id throws", "[NewickTree][build][error]") {
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "root", .branch_length = 0.0, .edge_id = std::nullopt},
	    {.node_id = 0, .parent_id = 0, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},  // Duplicate
	};

	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("Duplicate"));
}

TEST_CASE("NewickTree build disconnected node throws", "[NewickTree][build][error]") {
	// Node 2 is not connected to root
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 0, .parent_id = std::nullopt, .name = "root", .branch_length = 0.0, .edge_id = std::nullopt},
	    {.node_id = 1, .parent_id = 0, .name = "A", .branch_length = 0.1, .edge_id = std::nullopt},
	    {.node_id = 2, .parent_id = 3, .name = "B", .branch_length = 0.2, .edge_id = std::nullopt},  // Parent exists but not connected
	    {.node_id = 3, .parent_id = 2, .name = "C", .branch_length = 0.3, .edge_id = std::nullopt},  // Creates isolated cycle
	};

	REQUIRE_THROWS_WITH(miint::NewickTree::build(nodes), ContainsSubstring("not reachable"));
}

TEST_CASE("NewickTree build single node tree", "[NewickTree][build]") {
	std::vector<miint::NodeInput> nodes = {
	    {.node_id = 42, .parent_id = std::nullopt, .name = "lonely", .branch_length = 1.0, .edge_id = 0},
	};

	auto tree = miint::NewickTree::build(nodes);

	REQUIRE(tree.num_nodes() == 1);
	REQUIRE(tree.num_tips() == 1);
	REQUIRE(tree.name(tree.root()) == "lonely");
	REQUIRE(tree.edge_id(tree.root()) == 0);
}
