#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>
#include <string>
#include <unordered_map>
#include <vector>
#include "NewickTree.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

namespace {

uint32_t idx(const miint::NewickTree &tree, const std::string &name) {
	auto found = tree.find_node_by_name(name);
	REQUIRE(found.has_value());
	return found.value();
}

// The ancestral state at a given internal node (fails if none).
miint::AncestralStateBM at(const std::vector<miint::AncestralStateBM> &v, uint32_t node) {
	for (const auto &s : v) {
		if (s.node == node) {
			return s;
		}
	}
	FAIL("no ancestral state at node " << node);
	return {}; // unreachable
}

// The PIC ancestral_estimate (down-pass) at a node, for cross-checks.
double pic_estimate(const std::vector<miint::IndependentContrast> &cs, uint32_t node) {
	for (const auto &c : cs) {
		if (c.node == node) {
			return c.ancestral_estimate;
		}
	}
	FAIL("no PIC contrast at node " << node);
	return 0.0;
}

constexpr double kZ975 = 1.959963984540054; // qnorm(0.975)

} // namespace

// Golden tree (reused from the PIC test so the root estimate cross-checks):
//
//              root
//             /    \
//         E:2       C:4
//        /   \
//      A:1   B:3
//
// Trait X: A=2, B=6, C=1. Hand-derived (Felsenstein 1985 down-pass + up-pass):
//   Down-pass:  m_E = 3,      V_E = 3/4  ; m_root = 59/27, V_root = 44/27.
//   Up message to E: N(mean=1, var=6)  => E marginal:
//       estimate_E = 25/9  = 2.777... , var_rel_E = 2/3   = 0.666...
//   Root (no parent) marginal == its down-pass:
//       estimate_root = 59/27 = 2.185..., var_rel_root = 44/27 = 1.629...
// E's marginal (25/9) DIFFERS from its down-pass value (3.0): the up-pass matters.
TEST_CASE("ancestral_states_bm full marginal on the PIC golden tree", "[NewickTree][asr]") {
	auto tree = miint::NewickTree::parse("((A:1,B:3)E:2,C:4)root;");
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	auto states = tree.ancestral_states_bm(x);

	// n=3 tips -> exactly 2 internal nodes (E, root).
	REQUIRE(states.size() == 2);

	auto e = at(states, idx(tree, "E"));
	auto r = at(states, idx(tree, "root"));

	// Estimates are convention-free (no dependence on the sigma^2 divisor).
	REQUIRE(e.estimate == Approx(25.0 / 9.0).margin(1e-9));
	REQUIRE(r.estimate == Approx(59.0 / 27.0).margin(1e-9));

	// Cross-checks against PIC on the same tree:
	auto pic = tree.independent_contrasts(x);
	//  - the root marginal equals PIC's root ancestral_estimate (free oracle);
	REQUIRE(r.estimate == Approx(pic_estimate(pic, idx(tree, "root"))).margin(1e-9));
	//  - but E's marginal is NOT PIC's E down-pass value (3.0): proves the up-pass runs.
	REQUIRE(pic_estimate(pic, idx(tree, "E")) == Approx(3.0).margin(1e-9));
	REQUIRE(e.estimate != Approx(3.0).margin(1e-6));

	// Structural-variance ratio is convention-free (sigma^2 cancels):
	//   var_root / var_E = (44/27) / (2/3) = 22/9.
	REQUIRE(r.variance / e.variance == Approx(22.0 / 9.0).margin(1e-9));

	// CIs are symmetric about the estimate and use the 95% normal quantile.
	for (const auto &s : {e, r}) {
		REQUIRE(s.ci_low < s.estimate);
		REQUIRE(s.estimate < s.ci_high);
		REQUIRE((s.ci_high - s.estimate) == Approx(s.estimate - s.ci_low).margin(1e-12));
		REQUIRE((s.ci_high - s.estimate) == Approx(kZ975 * std::sqrt(s.variance)).margin(1e-9));
	}
}

// Multifurcation: a 3-way star. The root is the only internal node; its estimate
// is the PRECISION-weighted mean of the tips (weight 1/branch_length). This is
// the single most important guard against the m=2 "cross-multiply by the sibling
// variance" shortcut, which would give a different (wrong) answer here.
TEST_CASE("ancestral_states_bm precision-weighted mean on a multifurcating star", "[NewickTree][asr]") {
	auto tree = miint::NewickTree::parse("(A:1,B:2,C:4)root;");
	std::unordered_map<std::string, double> x({{"A", 4.0}, {"B", 6.0}, {"C", 8.0}});
	auto states = tree.ancestral_states_bm(x);

	REQUIRE(states.size() == 1);
	// (4/1 + 6/2 + 8/4) / (1/1 + 1/2 + 1/4) = 9 / 1.75 = 36/7.
	REQUIRE(at(states, idx(tree, "root")).estimate == Approx(36.0 / 7.0).margin(1e-9));
}

// Equal-length star -> the plain arithmetic mean (special case of precision weighting).
TEST_CASE("ancestral_states_bm equal-length star gives the arithmetic mean", "[NewickTree][asr]") {
	auto tree = miint::NewickTree::parse("(A:1,B:1,C:1,D:1)root;");
	std::unordered_map<std::string, double> x({{"A", 1.0}, {"B", 2.0}, {"C", 3.0}, {"D", 10.0}});
	auto states = tree.ancestral_states_bm(x);

	REQUIRE(states.size() == 1);
	REQUIRE(at(states, idx(tree, "root")).estimate == Approx(4.0).margin(1e-9)); // (1+2+3+10)/4
}

// A unifurcation (internal node with one child) carries no information and would
// leave the leave-one-out cavity empty; it is rejected with a shear_tree hint.
TEST_CASE("ancestral_states_bm rejects a unifurcation", "[NewickTree][asr]") {
	// E has a single child A.
	auto tree = miint::NewickTree::parse("((A:1)E:2,B:3)root;");
	std::unordered_map<std::string, double> x({{"A", 1.0}, {"B", 2.0}});
	REQUIRE_THROWS_WITH(tree.ancestral_states_bm(x), ContainsSubstring("shear_tree"));
}

// A NaN/infinite trait value would silently poison every estimate and the shared
// REML rate (the whole output), so it must fail loud.
TEST_CASE("ancestral_states_bm rejects a non-finite trait value", "[NewickTree][asr]") {
	auto tree = miint::NewickTree::parse("((A:1,B:3)E:2,C:4)root;");
	std::unordered_map<std::string, double> nan_x(
	    {{"A", std::numeric_limits<double>::quiet_NaN()}, {"B", 6.0}, {"C", 1.0}});
	REQUIRE_THROWS_WITH(tree.ancestral_states_bm(nan_x), ContainsSubstring("non-finite"));
	std::unordered_map<std::string, double> inf_x(
	    {{"A", 2.0}, {"B", std::numeric_limits<double>::infinity()}, {"C", 1.0}});
	REQUIRE_THROWS_WITH(tree.ancestral_states_bm(inf_x), ContainsSubstring("non-finite"));
}

// The batch overload (one scaffold reused across traits) must produce identical
// results to the single-trait overload called per trait.
TEST_CASE("ancestral_states_bm batch matches per-trait single calls", "[NewickTree][asr]") {
	auto tree = miint::NewickTree::parse("((A:1,B:3)E:2,C:4)root;");
	std::unordered_map<std::string, double> x({{"A", 2.0}, {"B", 6.0}, {"C", 1.0}});
	std::unordered_map<std::string, double> y({{"A", 5.0}, {"B", 1.0}, {"C", 0.0}});

	auto batch = tree.ancestral_states_bm(std::vector<std::unordered_map<std::string, double>> {x, y});
	REQUIRE(batch.size() == 2);

	auto sx = tree.ancestral_states_bm(x);
	auto sy = tree.ancestral_states_bm(y);
	REQUIRE(batch[0].size() == sx.size());
	REQUIRE(batch[1].size() == sy.size());
	for (size_t i = 0; i < sx.size(); i++) {
		REQUIRE(batch[0][i].node == sx[i].node);
		REQUIRE(batch[0][i].estimate == Approx(sx[i].estimate).margin(1e-12));
		REQUIRE(batch[0][i].variance == Approx(sx[i].variance).margin(1e-12));
	}
	for (size_t i = 0; i < sy.size(); i++) {
		REQUIRE(batch[1][i].node == sy[i].node);
		REQUIRE(batch[1][i].estimate == Approx(sy[i].estimate).margin(1e-12));
	}
}
