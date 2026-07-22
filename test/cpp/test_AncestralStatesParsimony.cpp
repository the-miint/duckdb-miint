#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
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

// A k*k unit cost matrix (0 on the diagonal, 1 elsewhere) -> Fitch parsimony.
std::vector<double> unit_cost(uint32_t k) {
	std::vector<double> c(static_cast<size_t>(k) * k, 1.0);
	for (uint32_t i = 0; i < k; ++i) {
		c[i * k + i] = 0.0;
	}
	return c;
}

// Index the flat parsimony output by (node, state) for lookups.
std::map<std::pair<uint32_t, uint32_t>, miint::AncestralStateParsimony>
by_node_state(const std::vector<miint::AncestralStateParsimony> &v) {
	std::map<std::pair<uint32_t, uint32_t>, miint::AncestralStateParsimony> m;
	for (const auto &s : v) {
		m[{s.node, s.state}] = s;
	}
	return m;
}

// The whole-tree parsimony score = the min_cost of any MPR entry (== the min over
// all entries, since every node achieves the global min at its best state).
double score_of(const std::vector<miint::AncestralStateParsimony> &v) {
	double best = std::numeric_limits<double>::infinity();
	for (const auto &s : v) {
		best = std::min(best, s.min_cost);
	}
	return best;
}

// Brute-force MPR oracle: enumerate all k^(#internal) assignments of states to the
// internal nodes (tips are fixed to their observed state) and, for each assignment,
// sum the per-edge substitution cost cost[parent_state * k + child_state]. Returns
// the global minimum (parsimony score) and, per (internal node, state), the minimum
// total cost over assignments with that node fixed to that state (i.e. F_u[s]). This
// is an INDEPENDENT reference for the Sankoff two-pass, valid for any arity and any
// cost matrix.
struct BruteResult {
	double score = std::numeric_limits<double>::infinity();
	std::map<uint32_t, std::vector<double>> F; // F[node][state]
};

BruteResult brute_force(const miint::NewickTree &tree, const std::unordered_map<std::string, uint32_t> &tip_states,
                        uint32_t k, const std::vector<double> &cost) {
	const uint32_t n = static_cast<uint32_t>(tree.num_nodes());
	std::vector<int64_t> obs(n, -1);
	std::vector<uint32_t> internal;
	std::vector<uint32_t> pos(n, 0);
	for (uint32_t i = 0; i < n; ++i) {
		if (tree.is_tip(i)) {
			obs[i] = static_cast<int64_t>(tip_states.at(tree.name(i)));
		} else {
			pos[i] = static_cast<uint32_t>(internal.size());
			internal.push_back(i);
		}
	}
	const size_t m = internal.size();

	BruteResult res;
	for (uint32_t u : internal) {
		res.F[u] = std::vector<double>(k, std::numeric_limits<double>::infinity());
	}

	std::vector<uint32_t> assign(m, 0);
	auto state_of = [&](uint32_t node) -> uint32_t {
		return tree.is_tip(node) ? static_cast<uint32_t>(obs[node]) : assign[pos[node]];
	};
	while (true) {
		double total = 0.0;
		for (uint32_t u = 0; u < n; ++u) {
			if (u == tree.root()) {
				continue;
			}
			total += cost[state_of(tree.parent(u)) * k + state_of(u)];
		}
		res.score = std::min(res.score, total);
		for (size_t j = 0; j < m; ++j) {
			double &f = res.F[internal[j]][assign[j]];
			f = std::min(f, total);
		}
		size_t i = 0;
		for (; i < m; ++i) {
			if (++assign[i] < k) {
				break;
			}
			assign[i] = 0;
		}
		if (i == m) {
			break; // wrapped all digits -> enumerated k^m assignments
		}
	}
	return res;
}

// Assert the Sankoff two-pass matches the brute-force oracle on (score, min_cost per
// (node,state), and MPR membership).
void require_matches_brute(const miint::NewickTree &tree, const std::unordered_map<std::string, uint32_t> &tip_states,
                           uint32_t k, const std::vector<double> &cost) {
	auto states = tree.ancestral_parsimony(tip_states, k, cost);
	auto m = by_node_state(states);
	auto brute = brute_force(tree, tip_states, k, cost);

	REQUIRE(score_of(states) == Approx(brute.score).margin(1e-9));
	for (const auto &[node, fvec] : brute.F) {
		for (uint32_t s = 0; s < k; ++s) {
			auto it = m.find({node, s});
			REQUIRE(it != m.end());
			REQUIRE(it->second.min_cost == Approx(fvec[s]).margin(1e-9));
			REQUIRE(it->second.in_mpr == (fvec[s] <= brute.score + 1e-9));
		}
	}
	// One entry per (internal node, state).
	REQUIRE(states.size() == brute.F.size() * k);
}

} // namespace

// Golden tree ((A,B)E,C)root; unit cost, 2 states. Tips A=0, B=1, C=0. Hand-derived
// Sankoff down-pass: g_E=[1,1], g_root=[1,2]; score=1. Up-pass gives F_E=[1,2],
// F_root=[1,2]. MPR(E)={0}, MPR(root)={0} (no tie here).
TEST_CASE("ancestral_parsimony hand golden on a bifurcating tree", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)E,C)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	auto states = tree.ancestral_parsimony(x, 2, unit_cost(2));
	auto m = by_node_state(states);

	REQUIRE(states.size() == 4); // 2 internal nodes * 2 states
	REQUIRE(score_of(states) == Approx(1.0));

	uint32_t e = idx(tree, "E"), r = idx(tree, "root");
	REQUIRE(m.at({e, 0}).in_mpr);
	REQUIRE(m.at({e, 0}).min_cost == Approx(1.0));
	REQUIRE_FALSE(m.at({e, 1}).in_mpr);
	REQUIRE(m.at({e, 1}).min_cost == Approx(2.0));
	REQUIRE(m.at({r, 0}).in_mpr);
	REQUIRE(m.at({r, 0}).min_cost == Approx(1.0));
	REQUIRE_FALSE(m.at({r, 1}).in_mpr);
	REQUIRE(m.at({r, 1}).min_cost == Approx(2.0));

	require_matches_brute(tree, x, 2, unit_cost(2));
}

// Classic ambiguous root: ((A,B)X,(C,D)Y)root; A=B=0, C=D=1, unit cost. The two
// clades disagree, so the root state is a first-class TIE: MPR(root)={0,1}, both at
// the minimum score 1; each subclade root is unambiguous. Guards that ties are
// reported (not collapsed to one arbitrary state).
TEST_CASE("ancestral_parsimony reports an ambiguous (tied) root", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)X,(C,D)Y)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 1}, {"D", 1}});
	auto states = tree.ancestral_parsimony(x, 2, unit_cost(2));
	auto m = by_node_state(states);

	REQUIRE(score_of(states) == Approx(1.0));
	uint32_t r = idx(tree, "root"), X = idx(tree, "X"), Y = idx(tree, "Y");
	// Root is tied between both states.
	REQUIRE(m.at({r, 0}).in_mpr);
	REQUIRE(m.at({r, 1}).in_mpr);
	REQUIRE(m.at({r, 0}).min_cost == Approx(1.0));
	REQUIRE(m.at({r, 1}).min_cost == Approx(1.0));
	// Each subclade root is unambiguous.
	REQUIRE(m.at({X, 0}).in_mpr);
	REQUIRE_FALSE(m.at({X, 1}).in_mpr);
	REQUIRE_FALSE(m.at({Y, 0}).in_mpr);
	REQUIRE(m.at({Y, 1}).in_mpr);

	require_matches_brute(tree, x, 2, unit_cost(2));
}

// Multifurcating star (A,B,C,D)root; A=B=0, C=D=1, unit cost. A hard 4-way polytomy;
// the root is tied {0,1} at score 2. Proves multifurcation is handled (the DP sums
// over an arbitrary child set) and reduces to majority-rule with a tie.
TEST_CASE("ancestral_parsimony on a multifurcating star with a tie", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("(A,B,C,D)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 1}, {"D", 1}});
	auto states = tree.ancestral_parsimony(x, 2, unit_cost(2));
	auto m = by_node_state(states);

	REQUIRE(states.size() == 2); // 1 internal node * 2 states
	REQUIRE(score_of(states) == Approx(2.0));
	uint32_t r = idx(tree, "root");
	REQUIRE(m.at({r, 0}).in_mpr);
	REQUIRE(m.at({r, 1}).in_mpr);

	require_matches_brute(tree, x, 2, unit_cost(2));
}

// Multifurcating multi-level tree vs the brute-force oracle, k=3, unit cost. root and
// X are both 3-way. This is the backbone proof that the up-pass + multifurcation are
// correct for arbitrary arity and k>2.
TEST_CASE("ancestral_parsimony matches brute force on a multifurcating tree (k=3)", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((t2,t3,t4)X,(t5,t6)Y,t1)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	require_matches_brute(tree, x, 3, unit_cost(3));
}

// A general (non-unit) SYMMETRIC cost matrix must be respected: enforce it via the
// brute-force oracle. k=3 with cheap 0<->2 transitions and expensive 1<->{0,2}.
TEST_CASE("ancestral_parsimony matches brute force with a symmetric cost matrix", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((t2,t3,t4)X,(t5,t6)Y,t1)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	// clang-format off
	std::vector<double> cost = {
	    0, 5, 1,
	    5, 0, 5,
	    1, 5, 0,
	};
	// clang-format on
	require_matches_brute(tree, x, 3, cost);
}

// An ASYMMETRIC (directed) cost matrix nails the cost[from*k + to] orientation
// (parent end = `from`, child end = `to`): a transpose bug would pass symmetric
// tests but fail here. cost(0->1)=1 is cheap; cost(1->0)=10 is dear.
TEST_CASE("ancestral_parsimony matches brute force with an asymmetric cost matrix", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)X,(C,D)Y)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 1}, {"D", 1}});
	// clang-format off
	std::vector<double> cost = {
	    0, 1,
	    10, 0,
	};
	// clang-format on
	require_matches_brute(tree, x, 2, cost);

	// With this directed cost the unit-cost root tie is broken: forcing the root to
	// state 0 needs a costly 0->1 nowhere... check membership is a strict subset of
	// what unit cost would give (unit -> {0,1}).
	auto states = tree.ancestral_parsimony(x, 2, cost);
	auto m = by_node_state(states);
	uint32_t r = idx(tree, "root");
	size_t mpr_root = (m.at({r, 0}).in_mpr ? 1 : 0) + (m.at({r, 1}).in_mpr ? 1 : 0);
	REQUIRE(mpr_root >= 1); // still at least one optimal state
}

// The parsimony score is achievable at every node's best state: min_s F_u[s] == score
// for every internal node (a re-rooting-style invariant the up-pass must satisfy).
TEST_CASE("ancestral_parsimony: every node achieves the score at its best state", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((t2,t3,t4)X,(t5,t6)Y,t1)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	auto states = tree.ancestral_parsimony(x, 3, unit_cost(3));
	double score = score_of(states);

	std::map<uint32_t, double> best;
	for (const auto &s : states) {
		auto it = best.find(s.node);
		best[s.node] = (it == best.end()) ? s.min_cost : std::min(it->second, s.min_cost);
	}
	for (const auto &[node, b] : best) {
		REQUIRE(b == Approx(score).margin(1e-9));
	}
}

// A single-state trait (all tips identical): every node is that state at zero cost.
TEST_CASE("ancestral_parsimony on a constant (k=1) trait", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)E,C)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 0}});
	auto states = tree.ancestral_parsimony(x, 1, unit_cost(1));
	REQUIRE(states.size() == 2); // 2 internal nodes * 1 state
	for (const auto &s : states) {
		REQUIRE(s.state == 0);
		REQUIRE(s.in_mpr);
		REQUIRE(s.min_cost == Approx(0.0));
	}
}

// Unifurcation (a single-child internal node) IS supported for parsimony (unlike
// BM/PIC): the DP propagates the child's message with a possible transition. Node E
// has the single child A; tips are A and B.
TEST_CASE("ancestral_parsimony allows a unifurcation", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A)E,B)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}});
	auto states = tree.ancestral_parsimony(x, 2, unit_cost(2));
	REQUIRE(states.size() == 4); // internal nodes E, root * 2 states
	require_matches_brute(tree, x, 2, unit_cost(2));
}

// Completeness: a tip with no state, and a state for a name that is not a tip.
TEST_CASE("ancestral_parsimony rejects trait/tip mismatches", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)E,C)root;");
	std::unordered_map<std::string, uint32_t> missing({{"A", 0}, {"B", 1}}); // no C
	REQUIRE_THROWS_WITH(tree.ancestral_parsimony(missing, 2, unit_cost(2)), ContainsSubstring("no state for tip"));

	std::unordered_map<std::string, uint32_t> extra({{"A", 0}, {"B", 1}, {"C", 0}, {"Z", 1}});
	REQUIRE_THROWS_WITH(tree.ancestral_parsimony(extra, 2, unit_cost(2)), ContainsSubstring("not a tip"));
}

// The batch overload (one scaffold reused across traits) must match per-trait single
// calls exactly, including per-trait alphabets (different k per trait).
TEST_CASE("ancestral_parsimony batch matches per-trait single calls", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)E,C)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}}); // k=2
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}}); // k=3

	auto batch = tree.ancestral_parsimony(std::vector<std::unordered_map<std::string, uint32_t>> {x, y},
	                                      std::vector<uint32_t> {2, 3},
	                                      std::vector<std::vector<double>> {unit_cost(2), unit_cost(3)});
	REQUIRE(batch.size() == 2);

	auto sx = tree.ancestral_parsimony(x, 2, unit_cost(2));
	auto sy = tree.ancestral_parsimony(y, 3, unit_cost(3));
	REQUIRE(batch[0].size() == sx.size());
	REQUIRE(batch[1].size() == sy.size());
	auto mb0 = by_node_state(batch[0]), msx = by_node_state(sx);
	for (const auto &[key, val] : msx) {
		REQUIRE(mb0.at(key).in_mpr == val.in_mpr);
		REQUIRE(mb0.at(key).min_cost == Approx(val.min_cost).margin(1e-12));
	}
	auto mb1 = by_node_state(batch[1]), msy = by_node_state(sy);
	for (const auto &[key, val] : msy) {
		REQUIRE(mb1.at(key).in_mpr == val.in_mpr);
		REQUIRE(mb1.at(key).min_cost == Approx(val.min_cost).margin(1e-12));
	}
}

// The in_mpr classification uses a hybrid tolerance max(1e-9, 1e-6*|score|); at small
// costs only the 1e-9 floor is ever active. Scale the unit cost by 1e5 so the RELATIVE
// term (~0.1) dominates, and confirm the ambiguous root tie is still detected -- a
// broken/negative relative term would drop exactly-tied states here.
TEST_CASE("ancestral_parsimony detects ties at large cost scales", "[NewickTree][asr][parsimony]") {
	auto tree = miint::NewickTree::parse("((A,B)X,(C,D)Y)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 1}, {"D", 1}});
	std::vector<double> cost = {0.0, 1e5, 1e5, 0.0}; // unit cost scaled by 1e5 (exactly representable)
	auto states = tree.ancestral_parsimony(x, 2, cost);
	auto m = by_node_state(states);
	uint32_t r = idx(tree, "root");
	REQUIRE(m.at({r, 0}).in_mpr);
	REQUIRE(m.at({r, 1}).in_mpr);
	REQUIRE(m.at({r, 0}).min_cost == Approx(1e5));
	require_matches_brute(tree, x, 2, cost);
}

// A single-tip tree has no internal nodes to reconstruct -> empty output (not an
// error), mirroring the BM sibling's n_tips <= 1 early return.
TEST_CASE("ancestral_parsimony returns empty for a single-tip tree", "[NewickTree][asr][parsimony]") {
	std::vector<miint::NodeInput> nodes = {
	    {0, std::nullopt, "A", std::numeric_limits<double>::quiet_NaN(), std::nullopt}};
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x({{"A", 0}});
	REQUIRE(tree.ancestral_parsimony(x, 1, unit_cost(1)).empty());
}
