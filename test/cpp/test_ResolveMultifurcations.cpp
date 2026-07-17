#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <vector>
#include "NewickTree.hpp"

namespace {

// Resolve a node index by name (fails the test if the name is absent).
uint32_t idx(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

// Postcondition of resolve_multifurcations: no node has more than two children.
// (Unifurcations are intentionally left alone, so the check is <= 2, not == 2.)
void require_no_multifurcations(const miint::NewickTree &t) {
	for (uint32_t i = 0; i < t.num_nodes(); ++i) {
		if (!t.is_tip(i)) {
			INFO("internal node " << i << " name='" << t.name(i) << "' has " << t.children(i).size() << " children");
			REQUIRE(t.children(i).size() <= 2);
		}
	}
}

// The one child of `node` whose index is not `other` (for a 2-child node).
uint32_t sibling_of(const miint::NewickTree &t, uint32_t node, uint32_t other) {
	const auto &ch = t.children(node);
	REQUIRE(ch.size() == 2);
	return ch[0] == other ? ch[1] : ch[0];
}

// A fully bifurcating tree with distinct branch lengths on every edge.
const char *kBinaryTree = "(((A:1,B:2)G:3,C:4)H:5,(D:6,E:7)I:8)root;";

} // namespace

// ============================================================================
// Pass-through: nothing to resolve
// ============================================================================

TEST_CASE("resolve leaves an already-bifurcating tree byte-identical", "[NewickTree][resolve]") {
	auto tree = miint::NewickTree::parse(kBinaryTree);
	auto before = tree.to_newick();

	auto out = tree.resolve_multifurcations();

	// WHY: a tree with no polytomy must round-trip unchanged — original indices,
	// child order, names, and branch lengths are all preserved, so the Newick
	// serialization is identical.
	REQUIRE(out.num_nodes() == tree.num_nodes());
	REQUIRE(out.num_tips() == tree.num_tips());
	REQUIRE(out.to_newick() == before);
	require_no_multifurcations(out);
}

TEST_CASE("resolve leaves a single tip unchanged", "[NewickTree][resolve]") {
	auto tree = miint::NewickTree::parse("A;");
	auto out = tree.resolve_multifurcations();
	REQUIRE(out.num_nodes() == 1);
	REQUIRE(out.num_tips() == 1);
	REQUIRE(out.is_tip(out.root()));
}

TEST_CASE("resolve leaves a unifurcation alone (polytomies only)", "[NewickTree][resolve]") {
	// root has a single child X, which is a real bifurcation. resolve() must not
	// collapse the unifurcation (that is shear()'s job), only split polytomies.
	auto tree = miint::NewickTree::parse("((A:1,B:2)X:3)root;");
	auto out = tree.resolve_multifurcations();
	REQUIRE(out.num_nodes() == tree.num_nodes());
	REQUIRE(out.children(out.root()).size() == 1);
}

// ============================================================================
// Resolving polytomies into a deterministic left-comb
// ============================================================================

TEST_CASE("resolve trifurcation -> bifurcation with a zero-length connector", "[NewickTree][resolve]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:3)root;");
	auto out = tree.resolve_multifurcations();

	REQUIRE(out.num_tips() == 3);
	REQUIRE(out.num_nodes() == 5); // root, A, B, C, + 1 connector (m-2 = 1)
	require_no_multifurcations(out);

	auto root = out.root();
	REQUIRE(out.name(root) == "root");
	REQUIRE(out.children(root).size() == 2);

	uint32_t a = idx(out, "A");
	uint32_t b = idx(out, "B");
	uint32_t c = idx(out, "C");

	// Left-comb, child order preserved: root keeps A (first child) plus the new
	// connector, in that order.
	REQUIRE(out.children(root)[0] == a);
	REQUIRE(out.parent(a) == root);
	REQUIRE(out.branch_length(a) == Catch::Approx(1.0));

	uint32_t connector = sibling_of(out, root, a);
	REQUIRE_FALSE(out.is_tip(connector));
	REQUIRE(out.name(connector).empty());
	REQUIRE_FALSE(out.edge_id(connector).has_value());
	REQUIRE(out.branch_length(connector) == Catch::Approx(0.0).margin(1e-12));

	// The connector holds the remaining two children with lengths preserved.
	REQUIRE(out.parent(b) == connector);
	REQUIRE(out.parent(c) == connector);
	REQUIRE(out.branch_length(b) == Catch::Approx(2.0));
	REQUIRE(out.branch_length(c) == Catch::Approx(3.0));
	REQUIRE(out.children(connector)[0] == b); // order preserved
}

TEST_CASE("resolve four-way multifurcation into a left-comb", "[NewickTree][resolve]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:3,D:4)root;");
	auto out = tree.resolve_multifurcations();

	REQUIRE(out.num_tips() == 4);
	REQUIRE(out.num_nodes() == 7); // 4 tips + root + (m-2 = 2) connectors
	require_no_multifurcations(out);

	uint32_t a = idx(out, "A"), b = idx(out, "B"), c = idx(out, "C"), d = idx(out, "D");
	auto root = out.root();

	// root -> [A, N1]
	REQUIRE(out.children(root)[0] == a);
	uint32_t n1 = sibling_of(out, root, a);
	REQUIRE_FALSE(out.is_tip(n1));
	REQUIRE(out.branch_length(n1) == Catch::Approx(0.0).margin(1e-12));

	// N1 -> [B, N2]
	REQUIRE(out.parent(b) == n1);
	REQUIRE(out.children(n1)[0] == b);
	uint32_t n2 = sibling_of(out, n1, b);
	REQUIRE_FALSE(out.is_tip(n2));
	REQUIRE(out.branch_length(n2) == Catch::Approx(0.0).margin(1e-12));

	// N2 -> [C, D], with C ahead of D (order preserved at the final connector too).
	REQUIRE(out.parent(c) == n2);
	REQUIRE(out.parent(d) == n2);
	REQUIRE(out.children(n2)[0] == c);
	REQUIRE(out.branch_length(c) == Catch::Approx(3.0));
	REQUIRE(out.branch_length(d) == Catch::Approx(4.0));
}

TEST_CASE("resolve a cladogram polytomy inserts explicit zero-length connectors", "[NewickTree][resolve]") {
	// A topology-only tree (no branch lengths anywhere). Resolving the root
	// polytomy inserts connectors with an EXPLICIT 0.0 length (not NaN): an
	// arbitrary resolution asserts "no evolutionary distance here", and PIC (the
	// downstream consumer) requires finite lengths. WHY pin this: it deliberately
	// differs from shear()'s all-NaN-preserving collapse, so it must be a
	// conscious choice, not an accident.
	auto tree = miint::NewickTree::parse("(A,B,C)root;");
	auto out = tree.resolve_multifurcations();

	REQUIRE(out.num_tips() == 3);
	REQUIRE(out.num_nodes() == 5);
	require_no_multifurcations(out);

	// Original tip edges stay unspecified (NaN); the new connector is 0.0.
	REQUIRE(std::isnan(out.branch_length(idx(out, "A"))));
	uint32_t connector = sibling_of(out, out.root(), idx(out, "A"));
	REQUIRE_FALSE(out.is_tip(connector));
	REQUIRE_FALSE(std::isnan(out.branch_length(connector)));
	REQUIRE(out.branch_length(connector) == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("resolve an internal polytomy, leaving its siblings untouched", "[NewickTree][resolve]") {
	// X has 3 children; root is already bifurcating (X, D).
	auto tree = miint::NewickTree::parse("((A:1,B:2,C:3)X:4,D:5)root;");
	auto out = tree.resolve_multifurcations();

	REQUIRE(out.num_tips() == 4);
	REQUIRE(out.num_nodes() == 7); // A,B,C,D + root + X + 1 connector
	require_no_multifurcations(out);

	// X keeps its name and branch length and now has exactly two children.
	uint32_t x = idx(out, "X");
	REQUIRE(out.branch_length(x) == Catch::Approx(4.0));
	REQUIRE(out.children(x).size() == 2);

	// The sibling side (D) is completely unaffected.
	uint32_t d = idx(out, "D");
	REQUIRE(out.parent(d) == out.root());
	REQUIRE(out.branch_length(d) == Catch::Approx(5.0));
}

// ============================================================================
// Determinism and scale
// ============================================================================

TEST_CASE("resolve is deterministic across runs", "[NewickTree][resolve]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:3,D:4,E:5)root;");
	auto first = tree.resolve_multifurcations().to_newick();
	auto second = tree.resolve_multifurcations().to_newick();
	REQUIRE(first == second);
}

TEST_CASE("resolve a large star tree fully bifurcates it", "[NewickTree][resolve]") {
	const int T = 1000;
	std::vector<miint::NodeInput> in;
	in.push_back({0, std::nullopt, "root", std::numeric_limits<double>::quiet_NaN(), std::nullopt});
	for (int i = 1; i <= T; ++i) {
		in.push_back({static_cast<int64_t>(i), std::optional<int64_t>(0), "t" + std::to_string(i), 1.0, std::nullopt});
	}
	auto tree = miint::NewickTree::build(in);

	auto out = tree.resolve_multifurcations();

	REQUIRE(out.num_tips() == static_cast<size_t>(T));
	// tips + root + (T-2) connectors.
	REQUIRE(out.num_nodes() == static_cast<size_t>(T + 1 + (T - 2)));
	require_no_multifurcations(out);
}
