#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include "NewickTree.hpp"

using Catch::Matchers::ContainsSubstring;

namespace {

uint32_t idx(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

// The contrast computed at a given internal node (fails if none).
miint::IndependentContrast at(const std::vector<miint::IndependentContrast> &cs, uint32_t node) {
	for (const auto &c : cs) {
		if (c.node == node) {
			return c;
		}
	}
	FAIL("no contrast at node " << node);
	return {}; // unreachable
}

// Hand-derived golden tree (Felsenstein 1985 recurrence, computed by hand):
//
//              root
//             /    \
//         E:2       C:4
//        /   \
//      A:1   B:3
//
// The internal edge E carries branch length 2; the asymmetric tip lengths
// (A:1, B:3) make the ancestral-value WEIGHTING observable: X_A is weighted by
// v_B and X_B by v_A, so a swapped weighting would be caught.
const char *kGolden = "((A:1,B:3)E:2,C:4)root;";

} // namespace

// ============================================================================
// Golden values (Felsenstein 1985), derived by hand from the recurrence:
//   contrast_k       = (X_i - X_j) / sqrt(v_i + v_j)
//   ancestral X_k    = (X_i*v_j + X_j*v_i) / (v_i + v_j)
//   extended v_k'    = v_k + (v_i*v_j)/(v_i + v_j)
//   contrast_var     = v_i + v_j
// ============================================================================

TEST_CASE("independent_contrasts matches the hand-derived golden", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse(kGolden);

	// Trait X: A=2, B=6, C=1
	// E: contrast = (2-6)/sqrt(1+3) = -4/2       = -2.0
	//    X_E      = (2*3 + 6*1)/4  = 12/4        =  3.0
	//    var      = 1+3                          =  4.0
	//    v_E'     = 2 + (1*3)/4    = 2 + 0.75    =  2.75
	// root: contrast = (3-1)/sqrt(2.75+4) = 2/sqrt(6.75) = 0.769800359...
	//       X_root   = (3*4 + 1*2.75)/6.75 = 14.75/6.75  = 2.185185185...
	//       var      = 6.75
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	auto cs = tree.independent_contrasts(x);

	// n=3 tips -> exactly n-1 = 2 contrasts (one per internal node).
	REQUIRE(cs.size() == 2);

	auto e = at(cs, idx(tree, "E"));
	REQUIRE(e.contrast == Catch::Approx(-2.0).margin(1e-9));
	REQUIRE(e.ancestral_estimate == Catch::Approx(3.0).margin(1e-9));
	REQUIRE(e.contrast_variance == Catch::Approx(4.0).margin(1e-9));

	auto r = at(cs, idx(tree, "root"));
	REQUIRE(r.contrast == Catch::Approx(0.769800359).margin(1e-9));
	REQUIRE(r.ancestral_estimate == Catch::Approx(2.185185185).margin(1e-9));
	REQUIRE(r.contrast_variance == Catch::Approx(6.75).margin(1e-9));
}

TEST_CASE("independent_contrasts: contrast_variance is trait-independent", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse(kGolden);

	// Trait Y: A=5, B=1, C=0
	// E:   contrast = (5-1)/2 = 2.0 ; X_E = (5*3+1*1)/4 = 4.0 ; var = 4.0
	// root:contrast = (4-0)/sqrt(6.75) = 1.539600718... ; X = 16/6.75 = 2.370370370 ; var = 6.75
	std::unordered_map<std::string, double> y({{"A", 5.0}, {"B", 1.0}, {"C", 0.0}});
	auto cs = tree.independent_contrasts(y);

	auto e = at(cs, idx(tree, "E"));
	auto r = at(cs, idx(tree, "root"));

	REQUIRE(e.contrast == Catch::Approx(2.0).margin(1e-9));
	REQUIRE(e.ancestral_estimate == Catch::Approx(4.0).margin(1e-9));
	REQUIRE(r.contrast == Catch::Approx(1.539600718).margin(1e-9));
	REQUIRE(r.ancestral_estimate == Catch::Approx(2.370370370).margin(1e-9));

	// WHY: variance depends only on the tree, not the trait. E stays 4.0 and root
	// stays 6.75 exactly as for trait X — this is what lets it be precomputed once.
	REQUIRE(e.contrast_variance == Catch::Approx(4.0).margin(1e-9));
	REQUIRE(r.contrast_variance == Catch::Approx(6.75).margin(1e-9));
}

TEST_CASE("independent_contrasts sign follows child order (deterministic)", "[NewickTree][pic]") {
	// A is the first child of E, so the contrast is X_A - X_B (negative here), not
	// the reverse. The sign is arbitrary in theory but must be stable in practice.
	auto tree = miint::NewickTree::parse(kGolden);
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	auto cs = tree.independent_contrasts(x);
	REQUIRE(at(cs, idx(tree, "E")).contrast < 0.0);
}

TEST_CASE("independent_contrasts count is n-1 for a larger bifurcating tree", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse("(((A:1,B:1)F:1,C:1)G:1,(D:1,E:1)H:1)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}, {"D", 4.0}, {"E", 5.0}});
	auto cs = tree.independent_contrasts(t);
	REQUIRE(cs.size() == 4); // 5 tips -> 4 contrasts
}

TEST_CASE("independent_contrasts batch matches per-trait single calls", "[NewickTree][pic]") {
	// The batch overload hoists the trait-independent work but must produce exactly
	// the same numbers as calling the single-trait method once per trait.
	auto tree = miint::NewickTree::parse(kGolden);
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	std::unordered_map<std::string, double> y({{"A", 5.0}, {"B", 1.0}, {"C", 0.0}});

	auto batch = tree.independent_contrasts(std::vector<std::unordered_map<std::string, double>> {x, y});
	REQUIRE(batch.size() == 2);

	auto same = [](const std::vector<miint::IndependentContrast> &a, const std::vector<miint::IndependentContrast> &b) {
		REQUIRE(a.size() == b.size());
		for (size_t i = 0; i < a.size(); i++) {
			REQUIRE(a[i].node == b[i].node);
			REQUIRE(a[i].contrast == Catch::Approx(b[i].contrast).margin(1e-12));
			REQUIRE(a[i].ancestral_estimate == Catch::Approx(b[i].ancestral_estimate).margin(1e-12));
			REQUIRE(a[i].contrast_variance == Catch::Approx(b[i].contrast_variance).margin(1e-12));
		}
	};
	same(batch[0], tree.independent_contrasts(x));
	same(batch[1], tree.independent_contrasts(y));
}

TEST_CASE("independent_contrasts batch validates each trait's completeness", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse(kGolden);
	std::unordered_map<std::string, double> ok({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	std::unordered_map<std::string, double> bad({{"A", 2.0}, {"B", 6.0}}); // C missing
	REQUIRE_THROWS_WITH(tree.independent_contrasts(std::vector<std::unordered_map<std::string, double>> {ok, bad}),
	                    ContainsSubstring("C"));
}

// ============================================================================
// Zero-length internal edges are allowed (the tree_resolve_multifurcations
// output must feed PIC). Only a zero *contrast variance* is an error.
// ============================================================================

TEST_CASE("independent_contrasts accepts a zero-length internal edge", "[NewickTree][pic]") {
	// E's own edge is 0 (as a resolver connector would be). Its children still
	// have positive lengths, so v_E' = 0 + (1*3)/4 = 0.75 and the root variance
	// 0.75+4 = 4.75 > 0 — valid.
	auto tree = miint::NewickTree::parse("((A:1,B:3)E:0,C:4)root;");
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});

	std::vector<miint::IndependentContrast> cs;
	REQUIRE_NOTHROW(cs = tree.independent_contrasts(x));
	REQUIRE(cs.size() == 2);
	REQUIRE(at(cs, idx(tree, "E")).contrast == Catch::Approx(-2.0).margin(1e-9));
	// The zero internal edge was actually used: root variance = 0.75 + 4 = 4.75.
	REQUIRE(at(cs, idx(tree, "root")).contrast_variance == Catch::Approx(4.75).margin(1e-9));
}

// ============================================================================
// Error paths (fail loud)
// ============================================================================

TEST_CASE("independent_contrasts rejects a polytomy, naming the resolver", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:3)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("tree_resolve_multifurcations"));
}

TEST_CASE("independent_contrasts rejects a unifurcation, naming shear_tree", "[NewickTree][pic]") {
	// root has a single child X: not bifurcating, but the remedy is collapsing the
	// unifurcation (shear_tree), not resolving a polytomy.
	auto tree = miint::NewickTree::parse("((A:1,B:2)X:3)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("shear_tree"));
}

TEST_CASE("independent_contrasts rejects a missing trait value", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse(kGolden);
	std::unordered_map<std::string, double> t({{"A", 2.0}, {"B", 6.0}}); // C missing
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("C"));
}

TEST_CASE("independent_contrasts rejects a trait value for a non-tip name", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse(kGolden);
	std::unordered_map<std::string, double> t({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}, {"Z", 9.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("not a tip"));
}

TEST_CASE("independent_contrasts rejects duplicate tip names", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse("((A:1,A:2)E:3,C:4)root;");
	std::unordered_map<std::string, double> t({{"A", 2.0}, {"C", 1.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("duplicate"));
}

TEST_CASE("independent_contrasts rejects NaN (unspecified) branch lengths", "[NewickTree][pic]") {
	// A cladogram has no branch lengths — PIC is undefined.
	auto tree = miint::NewickTree::parse("((A,B)E,C)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("branch length"));
}

TEST_CASE("independent_contrasts rejects infinite branch lengths", "[NewickTree][pic]") {
	// +Inf is finite-looking to an isnan()/negative check but poisons sqrt()/the
	// weighting into silent 0 / NaN output; it must be rejected up front. Built via
	// NodeInput so the infinity reaches validation regardless of parser behavior.
	std::vector<miint::NodeInput> in = {
	    {0, std::optional<int64_t>(2), "A", std::numeric_limits<double>::infinity(), std::nullopt},
	    {1, std::optional<int64_t>(2), "B", 1.0, std::nullopt},
	    {2, std::optional<int64_t>(4), "E", 1.0, std::nullopt},
	    {3, std::optional<int64_t>(4), "C", 1.0, std::nullopt},
	    {4, std::nullopt, "root", std::numeric_limits<double>::quiet_NaN(), std::nullopt},
	};
	auto tree = miint::NewickTree::build(in);
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("finite"));
}

TEST_CASE("independent_contrasts rejects negative branch lengths", "[NewickTree][pic]") {
	auto tree = miint::NewickTree::parse("((A:1,B:-1)E:2,C:4)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("negative"));
}

TEST_CASE("independent_contrasts rejects a zero contrast variance", "[NewickTree][pic]") {
	// Both of E's children have zero length -> v_i + v_j = 0 -> division by zero.
	auto tree = miint::NewickTree::parse("((A:0,B:0)E:1,C:1)root;");
	std::unordered_map<std::string, double> t({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}});
	REQUIRE_THROWS_WITH(tree.independent_contrasts(t), ContainsSubstring("zero"));
}
