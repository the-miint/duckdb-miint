#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <random>
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

std::map<std::pair<uint32_t, uint32_t>, double> post_by_node_state(const miint::AncestralMLResult &r) {
	std::map<std::pair<uint32_t, uint32_t>, double> m;
	for (const auto &s : r.states) {
		m[{s.node, s.state}] = s.probability;
	}
	return m;
}

// Closed-form ER transition probability P(t)[a][b] = 1/k + (delta_ab - 1/k) e^{-k r t}.
double er_P(uint32_t k, double r, double t, uint32_t a, uint32_t b) {
	double E = std::exp(-static_cast<double>(k) * r * t);
	return 1.0 / k + ((a == b ? 1.0 : 0.0) - 1.0 / k) * E;
}

// Dense k*k (row-major) matrix product Z = X*Y (test-only helper).
std::vector<double> mat_mul(const std::vector<double> &X, const std::vector<double> &Y, uint32_t k) {
	std::vector<double> Z(static_cast<size_t>(k) * k, 0.0);
	for (uint32_t i = 0; i < k; ++i) {
		for (uint32_t l = 0; l < k; ++l) {
			double x = X[i * k + l];
			for (uint32_t j = 0; j < k; ++j) {
				Z[i * k + j] += x * Y[l * k + j];
			}
		}
	}
	return Z;
}

// Test-only dense matrix exponential exp(M) via scaling-and-squaring + a Taylor series.
// This is INDEPENDENT of the production SelfAdjointEigenSolver path (it validates it):
// the brute-force SYM oracle below builds P(t) = exp(Q*t) with this, not with Eigen.
std::vector<double> mat_expm(std::vector<double> M, uint32_t k) {
	double norm = 0.0;
	for (double v : M) {
		norm = std::max(norm, std::abs(v));
	}
	int s = 0;
	while (norm > 0.5) {
		norm *= 0.5;
		++s;
	}
	double scale = std::ldexp(1.0, -s); // 1 / 2^s
	for (double &v : M) {
		v *= scale;
	}
	std::vector<double> E(static_cast<size_t>(k) * k, 0.0), term(static_cast<size_t>(k) * k, 0.0);
	for (uint32_t i = 0; i < k; ++i) {
		E[i * k + i] = 1.0;
		term[i * k + i] = 1.0;
	}
	for (int order = 1; order <= 24; ++order) {
		term = mat_mul(term, M, k);
		for (double &v : term) {
			v /= order;
		}
		for (size_t i = 0; i < E.size(); ++i) {
			E[i] += term[i];
		}
	}
	for (int i = 0; i < s; ++i) {
		E = mat_mul(E, E, k);
	}
	return E;
}

// Build a symmetric Mk (SYM) rate matrix Q from k(k-1)/2 off-diagonal rates in pair
// order (i,j), i<j; diagonal = negated row sum. Test-only reference construction.
std::vector<double> sym_Q(const std::vector<double> &rates, uint32_t k) {
	std::vector<double> Q(static_cast<size_t>(k) * k, 0.0);
	size_t p = 0;
	for (uint32_t i = 0; i < k; ++i) {
		for (uint32_t j = i + 1; j < k; ++j) {
			Q[i * k + j] = Q[j * k + i] = rates[p++];
		}
	}
	for (uint32_t i = 0; i < k; ++i) {
		double sum = 0.0;
		for (uint32_t j = 0; j < k; ++j) {
			if (j != i) {
				sum += Q[i * k + j];
			}
		}
		Q[i * k + i] = -sum;
	}
	return Q;
}

// Exact brute-force ML oracle: enumerate all k^(#internal) internal-state assignments
// (tips fixed to observed), sum the joint probability pi[root] * prod_edges
// P(parent_state, child_state; t) with a UNIFORM root prior pi=1/k. Returns log(Lik) and
// the exact per-(internal node, state) marginal posteriors. Independent of the pruning
// engine; `P(a,b,t)` is the caller's transition probability.
struct MLBrute {
	double logL = 0.0;
	std::map<uint32_t, std::vector<double>> post; // post[node][state]
};

template <typename Pfun>
MLBrute brute_force_ml_generic(const miint::NewickTree &tree,
                               const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k, Pfun &&P) {
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

	MLBrute res;
	for (uint32_t u : internal) {
		res.post[u] = std::vector<double>(k, 0.0);
	}

	std::vector<uint32_t> assign(m, 0);
	auto state_of = [&](uint32_t node) -> uint32_t {
		return tree.is_tip(node) ? static_cast<uint32_t>(obs[node]) : assign[pos[node]];
	};
	double Lik = 0.0;
	while (true) {
		double p = 1.0 / k; // pi[root], uniform
		for (uint32_t u = 0; u < n; ++u) {
			if (u == tree.root()) {
				continue;
			}
			p *= P(state_of(tree.parent(u)), state_of(u), tree.branch_length(u));
		}
		Lik += p;
		for (size_t j = 0; j < m; ++j) {
			res.post[internal[j]][assign[j]] += p;
		}
		size_t i = 0;
		for (; i < m; ++i) {
			if (++assign[i] < k) {
				break;
			}
			assign[i] = 0;
		}
		if (i == m) {
			break;
		}
	}
	for (uint32_t u : internal) {
		for (uint32_t s = 0; s < k; ++s) {
			res.post[u][s] /= Lik;
		}
	}
	res.logL = std::log(Lik);
	return res;
}

MLBrute brute_force_ml(const miint::NewickTree &tree, const std::unordered_map<std::string, uint32_t> &tip_states,
                       uint32_t k, double r) {
	return brute_force_ml_generic(tree, tip_states, k,
	                              [&](uint32_t a, uint32_t b, double t) { return er_P(k, r, t, a, b); });
}

// Build a general (ARD) rate matrix Q from k*(k-1) off-diagonal rates in row order (i,j),
// i!=j; diagonal = negated row sum. Test-only reference construction (may be asymmetric).
std::vector<double> ard_Q(const std::vector<double> &rates, uint32_t k) {
	std::vector<double> Q(static_cast<size_t>(k) * k, 0.0);
	size_t p = 0;
	for (uint32_t i = 0; i < k; ++i) {
		for (uint32_t j = 0; j < k; ++j) {
			if (j != i) {
				Q[i * k + j] = rates[p++];
			}
		}
	}
	for (uint32_t i = 0; i < k; ++i) {
		double sum = 0.0;
		for (uint32_t j = 0; j < k; ++j) {
			if (j != i) {
				sum += Q[i * k + j];
			}
		}
		Q[i * k + i] = -sum;
	}
	return Q;
}

// Brute force with P(t) = exp(Q*t) via the INDEPENDENT test expm, cached per branch length.
// Validates the production pruning (eigen-based SYM or Eigen-expm ARD) against a separate
// matrix-exponential reference. Works for symmetric (SYM) and general (ARD) Q alike.
MLBrute brute_force_ml_from_Q(const miint::NewickTree &tree,
                              const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
                              const std::vector<double> &Q) {
	std::map<double, std::vector<double>> cache;
	auto Pt = [&](double t) -> const std::vector<double> & {
		auto it = cache.find(t);
		if (it == cache.end()) {
			std::vector<double> Qt(Q);
			for (double &v : Qt) {
				v *= t;
			}
			it = cache.emplace(t, mat_expm(std::move(Qt), k)).first;
		}
		return it->second;
	};
	return brute_force_ml_generic(tree, tip_states, k,
	                              [&](uint32_t a, uint32_t b, double t) { return Pt(t)[a * k + b]; });
}

MLBrute brute_force_ml_sym(const miint::NewickTree &tree, const std::unordered_map<std::string, uint32_t> &tip_states,
                           uint32_t k, const std::vector<double> &rates) {
	return brute_force_ml_from_Q(tree, tip_states, k, sym_Q(rates, k));
}

MLBrute brute_force_ml_ard(const miint::NewickTree &tree, const std::unordered_map<std::string, uint32_t> &tip_states,
                           uint32_t k, const std::vector<double> &rates) {
	return brute_force_ml_from_Q(tree, tip_states, k, ard_Q(rates, k));
}

// Assert the ER pruning engine matches the brute-force oracle at a FIXED rate.
void require_matches_brute_ml(const miint::NewickTree &tree,
                              const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k, double r) {
	auto out = tree.ancestral_ml(tip_states, k, r);
	auto m = post_by_node_state(out);
	auto brute = brute_force_ml(tree, tip_states, k, r);

	REQUIRE(out.rate == Approx(r));
	REQUIRE(out.log_likelihood == Approx(brute.logL).margin(1e-9));
	for (const auto &[node, pvec] : brute.post) {
		double sum = 0.0;
		for (uint32_t s = 0; s < k; ++s) {
			auto it = m.find({node, s});
			REQUIRE(it != m.end());
			REQUIRE(it->second == Approx(pvec[s]).margin(1e-9));
			sum += it->second;
		}
		REQUIRE(sum == Approx(1.0).margin(1e-9)); // per-node posteriors normalize
	}
	REQUIRE(out.states.size() == brute.post.size() * k);
}

// Assert the SYM pruning engine matches the brute-force expm oracle at the FITTED rates:
// the oracle uses the fitted `rates` but computes P(t) independently (test expm), so this
// validates the eigen-based P(t) + matrix pruning + matrix marginal up-pass. Independent
// of optimizer quality (both sides evaluate at the same fitted rates).
void require_matches_brute_sym(const miint::NewickTree &tree,
                               const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k) {
	auto out = tree.ancestral_ml_sym(tip_states, k);
	REQUIRE(out.rates.size() == static_cast<size_t>(k) * (k - 1) / 2);
	REQUIRE(std::isnan(out.rate)); // SYM has no single scalar rate
	auto m = post_by_node_state(out);
	auto brute = brute_force_ml_sym(tree, tip_states, k, out.rates);

	REQUIRE(out.log_likelihood == Approx(brute.logL).margin(1e-9));
	for (const auto &[node, pvec] : brute.post) {
		double sum = 0.0;
		for (uint32_t s = 0; s < k; ++s) {
			auto it = m.find({node, s});
			REQUIRE(it != m.end());
			REQUIRE(it->second == Approx(pvec[s]).margin(1e-9));
			sum += it->second;
		}
		REQUIRE(sum == Approx(1.0).margin(1e-9));
	}
	REQUIRE(out.states.size() == brute.post.size() * k);
}

// Assert the ARD pruning engine matches the brute-force expm oracle at the FITTED rates.
// This validates the general (non-symmetric) P(t) = exp(Q*t) plus the asymmetric matrix
// messages -- crucially the up-pass uses P(t)^T (backward), so a P-vs-P^T mix-up would
// surface here whenever the fitted rates are asymmetric. Independent of optimizer quality.
void require_matches_brute_ard(const miint::NewickTree &tree,
                               const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k) {
	auto out = tree.ancestral_ml_ard(tip_states, k);
	REQUIRE(out.rates.size() == static_cast<size_t>(k) * (k - 1));
	REQUIRE(std::isnan(out.rate));
	auto m = post_by_node_state(out);
	auto brute = brute_force_ml_ard(tree, tip_states, k, out.rates);

	REQUIRE(out.log_likelihood == Approx(brute.logL).margin(1e-9));
	for (const auto &[node, pvec] : brute.post) {
		double sum = 0.0;
		for (uint32_t st = 0; st < k; ++st) {
			auto it = m.find({node, st});
			REQUIRE(it != m.end());
			REQUIRE(it->second == Approx(pvec[st]).margin(1e-9));
			sum += it->second;
		}
		REQUIRE(sum == Approx(1.0).margin(1e-9));
	}
	REQUIRE(out.states.size() == brute.post.size() * k);
}

std::vector<miint::NodeInput> make_balanced(int n_nodes, double bl) {
	std::vector<miint::NodeInput> nodes;
	for (int i = 0; i < n_nodes; ++i) {
		std::optional<int64_t> parent = (i == 0) ? std::nullopt : std::optional<int64_t>((i - 1) / 2);
		double b = (i == 0) ? std::numeric_limits<double>::quiet_NaN() : bl;
		nodes.push_back({i, parent, "n" + std::to_string(i), b, std::nullopt});
	}
	return nodes;
}

// A caterpillar (maximally unbalanced) tree with `n_tips` tips: internals 0..n_tips-2
// form a path (0 = root); each non-last internal has one tip child, the last has two.
std::vector<miint::NodeInput> make_caterpillar(int n_tips, double bl) {
	std::vector<miint::NodeInput> nodes;
	const int n_internal = n_tips - 1;
	for (int j = 0; j < n_internal; ++j) {
		std::optional<int64_t> parent = (j == 0) ? std::nullopt : std::optional<int64_t>(j - 1);
		double b = (j == 0) ? std::numeric_limits<double>::quiet_NaN() : bl;
		nodes.push_back({j, parent, "i" + std::to_string(j), b, std::nullopt});
	}
	int tipid = n_internal;
	for (int j = 0; j < n_internal - 1; ++j) {
		nodes.push_back({tipid, static_cast<int64_t>(j), "t" + std::to_string(tipid), bl, std::nullopt});
		tipid++;
	}
	nodes.push_back({tipid, static_cast<int64_t>(n_internal - 1), "t" + std::to_string(tipid), bl, std::nullopt});
	tipid++;
	nodes.push_back({tipid, static_cast<int64_t>(n_internal - 1), "t" + std::to_string(tipid), bl, std::nullopt});
	return nodes;
}

} // namespace

// Pruning + posteriors match the exact brute-force oracle at a fixed rate, on a
// bifurcating tree with unequal branch lengths.
TEST_CASE("ancestral_ml matches brute force on a bifurcating tree (fixed rate)", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	require_matches_brute_ml(tree, x, 2, 0.4);
	require_matches_brute_ml(tree, x, 2, 1.7);
}

// Multifurcating tree, k=3: proves the pruning + marginal up-pass are correct for
// arbitrary arity and k>2 against the exact oracle.
TEST_CASE("ancestral_ml matches brute force on a multifurcating tree (k=3)", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((t2:0.4,t3:0.9,t4:0.6)X:0.5,(t5:0.7,t6:1.1)Y:0.8,t1:0.3)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	require_matches_brute_ml(tree, x, 3, 0.6);
}

// Unifurcation (single-child internal node) is supported for ML.
TEST_CASE("ancestral_ml allows a unifurcation", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5)E:0.7,B:1.0)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}});
	require_matches_brute_ml(tree, x, 2, 0.8);
}

// At the FITTED rate (optimizer engaged), the reported posteriors + logL must still
// match the brute force evaluated at that same fitted rate: validates the whole
// pipeline (optimize -> prune -> posteriors) is self-consistent.
TEST_CASE("ancestral_ml fitted-rate output is consistent with brute force", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	auto out = tree.ancestral_ml(x, 2);
	REQUIRE(out.rate > 0.0);
	auto brute = brute_force_ml(tree, x, 2, out.rate);
	REQUIRE(out.log_likelihood == Approx(brute.logL).margin(1e-9));
	auto m = post_by_node_state(out);
	for (const auto &[node, pvec] : brute.post) {
		for (uint32_t s = 0; s < 2; ++s) {
			REQUIRE(m.at({node, s}) == Approx(pvec[s]).margin(1e-9));
		}
	}
}

// Analytical limit: monomorphic tips (all in one state) -> the fitted rate collapses
// toward 0 and every internal node's posterior concentrates on that state.
TEST_CASE("ancestral_ml monomorphic tips -> posterior ~1 everywhere", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 0}});
	auto out = tree.ancestral_ml(x, 2);
	auto m = post_by_node_state(out);
	REQUIRE(m.at({idx(tree, "E"), 0}) == Approx(1.0).margin(1e-4));
	REQUIRE(m.at({idx(tree, "root"), 0}) == Approx(1.0).margin(1e-4));
}

// Analytical limit: as rate -> infinity, P(t) -> uniform and every posterior -> the
// stationary prior (1/k). Exercised via a fixed huge rate.
TEST_CASE("ancestral_ml rate -> infinity -> posteriors approach the uniform prior", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	auto out = tree.ancestral_ml(x, 2, 1e6);
	for (const auto &s : out.states) {
		REQUIRE(s.probability == Approx(0.5).margin(1e-3));
	}
}

// Underflow / scaling: a deep caterpillar (2000 tips) must yield a finite logL and
// per-node posteriors that sum to 1. Without per-node rescaling the likelihood would
// underflow to 0 (logL = -inf) and this would fail.
TEST_CASE("ancestral_ml deep caterpillar stays finite (scaling works)", "[NewickTree][asr][ml]") {
	auto nodes = make_caterpillar(2000, 0.05);
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x;
	uint32_t t = 0;
	for (const auto &tip : tree.tip_names()) {
		x[tip] = (t++ % 2);
	}
	auto out = tree.ancestral_ml(x, 2, 0.3); // fixed rate (no need to optimize here)
	REQUIRE(std::isfinite(out.log_likelihood));
	std::map<uint32_t, double> sum;
	for (const auto &s : out.states) {
		sum[s.node] += s.probability;
	}
	for (const auto &[node, v] : sum) {
		REQUIRE(v == Approx(1.0).margin(1e-9));
	}
}

// The optimizer recovers a known rate from data SIMULATED under the ER model on a
// balanced 1024-tip tree (fixed seed -> deterministic). ML has finite-sample variance,
// so the tolerance is generous (within a factor of two).
TEST_CASE("ancestral_ml optimizer recovers a simulated rate", "[NewickTree][asr][ml]") {
	auto nodes = make_balanced(2047, 1.0);
	auto tree = miint::NewickTree::build(nodes);
	const uint32_t k = 2;
	const double r_true = 0.3;

	std::mt19937 rng(12345);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	const uint32_t n = static_cast<uint32_t>(tree.num_nodes());
	std::vector<uint32_t> state(n, 0);
	for (uint32_t u : tree.preorder()) {
		if (u == tree.root()) {
			state[u] = 0;
			continue;
		}
		uint32_t ps = state[tree.parent(u)];
		double E = std::exp(-static_cast<double>(k) * r_true * tree.branch_length(u));
		double p_stay = 1.0 / k + (static_cast<double>(k) - 1.0) / k * E;
		state[u] = (unif(rng) < p_stay) ? ps : (1 - ps);
	}
	std::unordered_map<std::string, uint32_t> x;
	for (uint32_t i = 0; i < n; ++i) {
		if (tree.is_tip(i)) {
			x[tree.name(i)] = state[i];
		}
	}

	auto out = tree.ancestral_ml(x, k);
	REQUIRE(out.rate > 0.5 * r_true);
	REQUIRE(out.rate < 2.0 * r_true);
}

// The fitted rate is a (local) maximum of the log-likelihood: logL(r_hat) is >= logL
// at half and double the fitted rate.
TEST_CASE("ancestral_ml fitted rate maximizes the log-likelihood", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	auto out = tree.ancestral_ml(x, 2);
	double lo = tree.ancestral_ml(x, 2, 0.5 * out.rate).log_likelihood;
	double hi = tree.ancestral_ml(x, 2, 2.0 * out.rate).log_likelihood;
	REQUIRE(out.log_likelihood >= lo - 1e-9);
	REQUIRE(out.log_likelihood >= hi - 1e-9);
}

// A constant (k=1) trait: one state, posterior 1 everywhere, logL 0.
TEST_CASE("ancestral_ml on a constant (k=1) trait", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 0}});
	auto out = tree.ancestral_ml(x, 1);
	REQUIRE(out.states.size() == 2); // 2 internal nodes * 1 state
	for (const auto &s : out.states) {
		REQUIRE(s.state == 0);
		REQUIRE(s.probability == Approx(1.0));
	}
	REQUIRE(out.log_likelihood == Approx(0.0).margin(1e-12));
}

// Single-tip tree -> no internal nodes -> empty output (not an error).
TEST_CASE("ancestral_ml returns empty for a single-tip tree", "[NewickTree][asr][ml]") {
	std::vector<miint::NodeInput> nodes = {
	    {0, std::nullopt, "A", std::numeric_limits<double>::quiet_NaN(), std::nullopt}};
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x({{"A", 0}});
	REQUIRE(tree.ancestral_ml(x, 1).states.empty());
}

// Rejections: incomplete/extra traits, and the branch-length policy.
TEST_CASE("ancestral_ml rejects bad inputs", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> missing({{"A", 0}, {"B", 1}}); // no C
	REQUIRE_THROWS_WITH(tree.ancestral_ml(missing, 2), ContainsSubstring("no state for tip"));

	std::unordered_map<std::string, uint32_t> extra({{"A", 0}, {"B", 1}, {"C", 0}, {"Z", 1}});
	REQUIRE_THROWS_WITH(tree.ancestral_ml(extra, 2), ContainsSubstring("not a tip"));

	// Zero-length tip edge is rejected (a hard zero would break the leave-one-out).
	auto ztip = miint::NewickTree::parse("((A:0.0,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 0}});
	REQUIRE_THROWS_WITH(ztip.ancestral_ml(y, 2), ContainsSubstring("positive"));

	// Non-finite (missing) branch length on a non-root edge is rejected.
	auto nanbl = miint::NewickTree::parse("((A,B:1.3)E:0.7,C:2.1)root;");
	REQUIRE_THROWS_WITH(nanbl.ancestral_ml(y, 2), ContainsSubstring("non-finite"));
}

// Batch matches per-trait single calls (shared scaffold; per-trait alphabet).
TEST_CASE("ancestral_ml batch matches per-trait single calls", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}}); // k=2
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}}); // k=3

	auto batch =
	    tree.ancestral_ml(std::vector<std::unordered_map<std::string, uint32_t>> {x, y}, std::vector<uint32_t> {2, 3});
	REQUIRE(batch.size() == 2);

	auto sx = tree.ancestral_ml(x, 2);
	auto sy = tree.ancestral_ml(y, 3);
	REQUIRE(batch[0].rate == Approx(sx.rate));
	REQUIRE(batch[0].log_likelihood == Approx(sx.log_likelihood));
	REQUIRE(batch[1].rate == Approx(sy.rate));
	auto mb0 = post_by_node_state(batch[0]), msx = post_by_node_state(sx);
	for (const auto &[key, val] : msx) {
		REQUIRE(mb0.at(key) == Approx(val).margin(1e-12));
	}
}

// An effectively-zero (but validly > 0) branch length whose exp(-k r t) rounds to
// EXACTLY 1.0 acts like a zero-length edge, hard-pinning a state and driving a message
// or the whole-tree likelihood to zero. Both failure shapes must fail LOUD (throw), not
// silently emit a NaN probability or a -inf log-likelihood.
TEST_CASE("ancestral_ml rejects effectively-zero branch lengths", "[NewickTree][asr][ml]") {
	// Shape (a): conflicting sibling tips with ~machine-zero branches at a FIXED small
	// rate -> the whole-tree likelihood is zero (caught after the final pruning pass).
	// (With a fitted rate the optimizer escapes to a huge rate that keeps the likelihood
	// tiny-but-finite -- a degenerate-but-valid ML fit -- so a fixed rate forces the
	// zero-likelihood regime the guard is there to catch.)
	auto ta = miint::NewickTree::parse("(A:1e-18,B:1e-18,C:1.0)root;");
	std::unordered_map<std::string, uint32_t> xa({{"A", 0}, {"B", 1}, {"C", 0}});
	REQUIRE_THROWS_WITH(ta.ancestral_ml(xa, 2, 0.5), ContainsSubstring("effectively-zero-length"));

	// Shape (b): a tiny tip branch + tiny internal branch zero out one state's message
	// while the log-likelihood stays finite -> caught at the up-pass division site.
	auto tb = miint::NewickTree::parse("((A:1e-18,D:1.0)E:1e-18,C:1.0)root;");
	std::unordered_map<std::string, uint32_t> xb({{"A", 0}, {"D", 0}, {"C", 0}});
	REQUIRE_THROWS_WITH(tb.ancestral_ml(xb, 2), ContainsSubstring("effectively-zero-length"));
}

// A zero-length INTERNAL edge is allowed (P(0) = identity); results still match the
// brute-force oracle.
TEST_CASE("ancestral_ml allows a zero-length internal edge", "[NewickTree][asr][ml]") {
	auto tree = miint::NewickTree::parse("((A:1.0,B:1.0)E:0.0,C:1.0)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});
	require_matches_brute_ml(tree, x, 2, 0.5);
}

// A mid-size caterpillar (18 tips, depth 17) checked EXACTLY against the brute-force
// oracle at a fixed rate: exercises the per-node rescaling path with an independent
// numeric reference (small-tree oracle cases barely trigger a rescale, and the deep
// caterpillar test only checks finiteness / sum-to-1, not the absolute logL value).
TEST_CASE("ancestral_ml matches brute force on a mid-size caterpillar (rescaling)", "[NewickTree][asr][ml]") {
	auto nodes = make_caterpillar(18, 0.3);
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x;
	uint32_t t = 0;
	for (const auto &tip : tree.tip_names()) {
		x[tip] = (t++ % 2);
	}
	require_matches_brute_ml(tree, x, 2, 0.5);
}

// ===========================================================================
// SYM (symmetric-rates) Mk-ML — Phase 3b
// ===========================================================================

// The general (eigen-based) SYM pruning + marginal up-pass matches the INDEPENDENT
// matrix-exponential brute-force oracle at the fitted rates, on a k=3 bifurcating tree
// with unequal branch lengths. This is the core correctness proof for the SYM P(t) =
// U diag(e^{lambda t}) U^T machinery and the k*k matrix messages.
TEST_CASE("ancestral_ml_sym matches brute force on a bifurcating k=3 tree", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	require_matches_brute_sym(tree, x, 3);
}

// Multifurcating k=3: proves the matrix pruning + marginal up-pass are correct for
// arbitrary arity against the exact expm oracle.
TEST_CASE("ancestral_ml_sym matches brute force on a multifurcating k=3 tree", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("((t2:0.4,t3:0.9,t4:0.6)X:0.5,(t5:0.7,t6:1.1)Y:0.8,t1:0.3)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	require_matches_brute_sym(tree, x, 3);
}

// Unifurcation (single-child internal node) is supported for SYM, like ER.
TEST_CASE("ancestral_ml_sym allows a unifurcation", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.0)E:0.4)F:0.7,C:1.2)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}});
	require_matches_brute_sym(tree, x, 3);
}

// A zero-length INTERNAL edge is allowed (P(0) = identity); results match the oracle.
TEST_CASE("ancestral_ml_sym allows a zero-length internal edge", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:1.0,B:1.0)E:0.0,C:1.0)F:0.5,D:1.3)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	require_matches_brute_sym(tree, x, 3);
}

// Analytical invariant: at k=2 there is exactly one free rate, so SYM reduces EXACTLY to
// ER. The fitted rate, log-likelihood, and posteriors must all match the ER overload.
TEST_CASE("ancestral_ml_sym reduces to ER at k=2", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}});

	auto sym = tree.ancestral_ml_sym(x, 2);
	auto er = tree.ancestral_ml(x, 2);

	REQUIRE(sym.rates.size() == 1);
	// Two different optimizers (SYM's Nelder-Mead vs ER's golden-section) maximize the
	// same 1-parameter likelihood; agreement to ~1e-4 confirms the reduction.
	REQUIRE(sym.rates[0] == Approx(er.rate).epsilon(1e-4));
	REQUIRE(sym.log_likelihood == Approx(er.log_likelihood).margin(1e-7));

	auto ms = post_by_node_state(sym), me = post_by_node_state(er);
	for (const auto &[key, val] : me) {
		REQUIRE(ms.at(key) == Approx(val).margin(1e-4));
	}
}

// Nested-model invariant: SYM has ER as a submodel (all rates equal), so its maximized
// log-likelihood must be >= the ER maximum on the same data. For k=3 the extra free rates
// generally let SYM fit strictly better.
TEST_CASE("ancestral_ml_sym log-likelihood dominates ER", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	auto sym = tree.ancestral_ml_sym(x, 3);
	auto er = tree.ancestral_ml(x, 3);
	REQUIRE(sym.log_likelihood >= er.log_likelihood - 1e-7);
}

// The Nelder-Mead fit recovers known (unequal) SYM rates from data SIMULATED under a
// known symmetric Q on a balanced 1024-tip tree (fixed seed -> deterministic). ML has
// finite-sample variance, so tolerances are generous.
TEST_CASE("ancestral_ml_sym optimizer recovers simulated rates", "[NewickTree][asr][ml][sym]") {
	auto nodes = make_balanced(2047, 1.0);
	auto tree = miint::NewickTree::build(nodes);
	const uint32_t k = 3;
	// True off-diagonal rates for pairs (0,1),(0,2),(1,2): deliberately unequal.
	std::vector<double> r_true = {0.20, 0.60, 0.35};
	std::vector<double> Q = sym_Q(r_true, k);

	std::mt19937 rng(24680);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	const uint32_t n = static_cast<uint32_t>(tree.num_nodes());
	std::vector<uint32_t> state(n, 0);
	std::map<double, std::vector<double>> Pcache;
	auto Pt = [&](double t) -> const std::vector<double> & {
		auto it = Pcache.find(t);
		if (it == Pcache.end()) {
			std::vector<double> Qt(Q);
			for (double &v : Qt) {
				v *= t;
			}
			it = Pcache.emplace(t, mat_expm(std::move(Qt), k)).first;
		}
		return it->second;
	};
	for (uint32_t u : tree.preorder()) {
		if (u == tree.root()) {
			state[u] = 0;
			continue;
		}
		uint32_t ps = state[tree.parent(u)];
		const std::vector<double> &P = Pt(tree.branch_length(u));
		double roll = unif(rng), cum = 0.0;
		uint32_t chosen = k - 1;
		for (uint32_t s = 0; s < k; ++s) {
			cum += P[ps * k + s];
			if (roll < cum) {
				chosen = s;
				break;
			}
		}
		state[u] = chosen;
	}
	std::unordered_map<std::string, uint32_t> x;
	for (uint32_t i = 0; i < n; ++i) {
		if (tree.is_tip(i)) {
			x[tree.name(i)] = state[i];
		}
	}

	auto out = tree.ancestral_ml_sym(x, k);
	REQUIRE(out.rates.size() == 3);
	for (size_t p = 0; p < 3; ++p) {
		REQUIRE(out.rates[p] > 0.4 * r_true[p]);
		REQUIRE(out.rates[p] < 2.5 * r_true[p]);
	}
	// (No brute-force cross-check here: 3^1023 internal assignments is infeasible; machinery
	// correctness is covered by the small-tree require_matches_brute_sym cases above.)
}

// k=1 (constant) trait: posterior 1 everywhere, logL 0, no rates.
TEST_CASE("ancestral_ml_sym on a constant (k=1) trait", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 0}});
	auto out = tree.ancestral_ml_sym(x, 1);
	REQUIRE(out.states.size() == 2); // 2 internal nodes * 1 state
	REQUIRE(out.rates.empty());
	for (const auto &s : out.states) {
		REQUIRE(s.state == 0);
		REQUIRE(s.probability == Approx(1.0));
	}
	REQUIRE(out.log_likelihood == Approx(0.0).margin(1e-12));
}

// Single-tip tree -> no internal nodes -> empty output (not an error), and rate is NaN
// (the SYM invariant; a symmetric rate matrix has no single scalar rate).
TEST_CASE("ancestral_ml_sym returns empty for a single-tip tree", "[NewickTree][asr][ml][sym]") {
	std::vector<miint::NodeInput> nodes = {
	    {0, std::nullopt, "A", std::numeric_limits<double>::quiet_NaN(), std::nullopt}};
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x({{"A", 0}});
	auto out = tree.ancestral_ml_sym(x, 1);
	REQUIRE(out.states.empty());
	REQUIRE(std::isnan(out.rate));
}

// Rejections mirror ER: incomplete/extra traits and the branch-length policy.
TEST_CASE("ancestral_ml_sym rejects bad inputs", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> missing({{"A", 0}, {"B", 1}, {"C", 2}}); // no D
	REQUIRE_THROWS_WITH(tree.ancestral_ml_sym(missing, 3), ContainsSubstring("no state for tip"));

	std::unordered_map<std::string, uint32_t> extra({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}, {"Z", 1}});
	REQUIRE_THROWS_WITH(tree.ancestral_ml_sym(extra, 3), ContainsSubstring("not a tip"));

	auto ztip = miint::NewickTree::parse("(((A:0.0,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	REQUIRE_THROWS_WITH(ztip.ancestral_ml_sym(y, 3), ContainsSubstring("positive"));

	auto nanbl = miint::NewickTree::parse("(((A,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	REQUIRE_THROWS_WITH(nanbl.ancestral_ml_sym(y, 3), ContainsSubstring("non-finite"));
}

// Batch matches per-trait single calls (shared scaffold; per-trait alphabet). `rate` is
// NaN for SYM so compare `rates`, logL, and posteriors instead.
TEST_CASE("ancestral_ml_sym batch matches per-trait single calls", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}, {"D", 1}}); // k=2
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}}); // k=3

	auto batch = tree.ancestral_ml_sym(std::vector<std::unordered_map<std::string, uint32_t>> {x, y},
	                                   std::vector<uint32_t> {2, 3});
	REQUIRE(batch.size() == 2);

	auto sx = tree.ancestral_ml_sym(x, 2);
	auto sy = tree.ancestral_ml_sym(y, 3);
	REQUIRE(batch[0].log_likelihood == Approx(sx.log_likelihood).margin(1e-9));
	REQUIRE(batch[1].log_likelihood == Approx(sy.log_likelihood).margin(1e-9));
	REQUIRE(batch[0].rates.size() == 1);
	REQUIRE(batch[0].rates[0] == Approx(sx.rates[0]).epsilon(1e-9));
	REQUIRE(batch[1].rates.size() == 3);
	auto mb1 = post_by_node_state(batch[1]), msy = post_by_node_state(sy);
	for (const auto &[key, val] : msy) {
		REQUIRE(mb1.at(key) == Approx(val).margin(1e-9));
	}
}

// Regression (review #1): the eigenbasis P(t)=U diag(e^{lambda t}) U^T matvec is not a sum
// of non-negative terms like ER's closed form, so at zero-length INTERNAL edges (P(0)=U U^T
// is only ~I to Eigen's orthonormality residual) a should-be-zero likelihood/message entry
// can round to a tiny NEGATIVE and leak into an emitted probability. sym_apply clamps its
// output to >= 0. Assert no emitted SYM probability is negative (and posteriors still sum to
// 1) on a tree that chains several zero-length internal edges with a k=4 trait -- the exact
// P(0) code path the clamp protects. Also cross-check the exact (expm) oracle, which is
// non-negative by construction.
TEST_CASE("ancestral_ml_sym never emits a negative probability (zero-length internal edges)",
          "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("((((A:1.0,B:1.0)I3:0.0,C:1.0)I2:0.0,D:1.0)I1:0.0,E:1.0)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}, {"E", 0}});
	auto out = tree.ancestral_ml_sym(x, 4);
	REQUIRE(!out.states.empty());
	for (const auto &s : out.states) {
		REQUIRE(s.probability >= 0.0);
		REQUIRE(s.probability <= 1.0 + 1e-9);
	}
	std::map<uint32_t, double> sum;
	for (const auto &s : out.states) {
		sum[s.node] += s.probability;
	}
	for (const auto &[node, v] : sum) {
		REQUIRE(v == Approx(1.0).margin(1e-9));
	}
	require_matches_brute_sym(tree, x, 4);
}

// Review #3: SYM fits k*(k-1)/2 rates by simplex search, so the state count is capped
// (k > 8 -> 28+ free rates degrade Nelder-Mead); such traits error, pointing to ER.
TEST_CASE("ancestral_ml_sym rejects too many states", "[NewickTree][asr][ml][sym]") {
	auto tree = miint::NewickTree::parse("(A:1,B:1,C:1,D:1,E:1,F:1,G:1,H:1,J:1)root;");
	std::unordered_map<std::string, uint32_t> x(
	    {{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}, {"E", 4}, {"F", 5}, {"G", 6}, {"H", 7}, {"J", 8}});
	REQUIRE_THROWS_WITH(tree.ancestral_ml_sym(x, 9), ContainsSubstring("model 'ER'"));
}

// ===========================================================================
// ARD (all-rates-different) Mk-ML -- Phase 3c
// ===========================================================================

// Core correctness: the general (non-symmetric) ARD pruning + marginal up-pass match the
// INDEPENDENT matrix-exponential brute-force oracle at the fitted rates. Because P(t) is
// not symmetric, the up-pass uses P(t)^T; a P-vs-P^T mix-up would surface here.
TEST_CASE("ancestral_ml_ard matches brute force on a bifurcating k=3 tree", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	require_matches_brute_ard(tree, x, 3);
}

// Multifurcating k=3 against the exact oracle (arbitrary arity, general Q).
TEST_CASE("ancestral_ml_ard matches brute force on a multifurcating k=3 tree", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("((t2:0.4,t3:0.9,t4:0.6)X:0.5,(t5:0.7,t6:1.1)Y:0.8,t1:0.3)root;");
	std::unordered_map<std::string, uint32_t> x({{"t1", 0}, {"t2", 1}, {"t3", 2}, {"t4", 0}, {"t5", 1}, {"t6", 2}});
	require_matches_brute_ard(tree, x, 3);
}

// Unifurcation and a zero-length internal edge are supported (and no negative probability).
TEST_CASE("ancestral_ml_ard allows a unifurcation", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.0)E:0.4)F:0.7,C:1.2)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}});
	require_matches_brute_ard(tree, x, 3);
}

TEST_CASE("ancestral_ml_ard allows a zero-length internal edge (no negative probability)",
          "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:1.0,B:1.0)E:0.0,C:1.0)F:0.5,D:1.3)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	auto out = tree.ancestral_ml_ard(x, 3);
	for (const auto &st : out.states) {
		REQUIRE(st.probability >= 0.0);
	}
	require_matches_brute_ard(tree, x, 3);
}

// Nested-model invariant: ER subset SYM subset ARD, so ARD's maximized log-likelihood must
// dominate ER's (and, up to optimizer slack, SYM's) on the same data.
TEST_CASE("ancestral_ml_ard log-likelihood dominates ER and SYM", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	auto ard = tree.ancestral_ml_ard(x, 3);
	auto sym = tree.ancestral_ml_sym(x, 3);
	auto er = tree.ancestral_ml(x, 3);
	REQUIRE(ard.log_likelihood >= er.log_likelihood - 1e-6);
	REQUIRE(ard.log_likelihood >= sym.log_likelihood - 1e-6);
}

// ARD-specific: the fit recovers ASYMMETRY. Data simulated under a strongly asymmetric Q
// on a balanced 1024-tip tree (state 0->others fast, reverse slow) must yield fitted rates
// that are correspondingly asymmetric (mean 0->* rate > mean *->0 rate), which SYM (forced
// q_ij=q_ji) cannot represent. Fixed seed -> deterministic.
TEST_CASE("ancestral_ml_ard recovers rate asymmetry", "[NewickTree][asr][ml][ard]") {
	auto nodes = make_balanced(2047, 1.0);
	auto tree = miint::NewickTree::build(nodes);
	const uint32_t k = 3;
	// True Q: leaving state 0 is fast (0.9), returning to 0 is slow (0.1); 1<->2 moderate.
	std::vector<double> r_true = {/*0->1*/ 0.9, /*0->2*/ 0.9, /*1->0*/ 0.1,
	                              /*1->2*/ 0.3, /*2->0*/ 0.1, /*2->1*/ 0.3};
	std::vector<double> Q = ard_Q(r_true, k);

	std::mt19937 rng(97531);
	std::uniform_real_distribution<double> unif(0.0, 1.0);
	const uint32_t n = static_cast<uint32_t>(tree.num_nodes());
	std::vector<uint32_t> state(n, 0);
	std::map<double, std::vector<double>> Pcache;
	auto Pt = [&](double t) -> const std::vector<double> & {
		auto it = Pcache.find(t);
		if (it == Pcache.end()) {
			std::vector<double> Qt(Q);
			for (double &v : Qt) {
				v *= t;
			}
			it = Pcache.emplace(t, mat_expm(std::move(Qt), k)).first;
		}
		return it->second;
	};
	for (uint32_t u : tree.preorder()) {
		if (u == tree.root()) {
			state[u] = 0;
			continue;
		}
		uint32_t ps = state[tree.parent(u)];
		const std::vector<double> &P = Pt(tree.branch_length(u));
		double roll = unif(rng), cum = 0.0;
		uint32_t chosen = k - 1;
		for (uint32_t st = 0; st < k; ++st) {
			cum += P[ps * k + st];
			if (roll < cum) {
				chosen = st;
				break;
			}
		}
		state[u] = chosen;
	}
	std::unordered_map<std::string, uint32_t> x;
	for (uint32_t i = 0; i < n; ++i) {
		if (tree.is_tip(i)) {
			x[tree.name(i)] = state[i];
		}
	}

	auto out = tree.ancestral_ml_ard(x, k);
	REQUIRE(out.rates.size() == 6);
	// Fitted rate order matches ard_Q: {0->1,0->2,1->0,1->2,2->0,2->1}.
	double leave0 = 0.5 * (out.rates[0] + out.rates[1]); // 0->1, 0->2
	double enter0 = 0.5 * (out.rates[2] + out.rates[4]); // 1->0, 2->0
	REQUIRE(leave0 > enter0);                            // the true asymmetry (fast out of 0, slow back) is recovered
}

// k=1 (constant) trait: posterior 1 everywhere, logL 0, no rates.
TEST_CASE("ancestral_ml_ard on a constant (k=1) trait", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("((A:0.5,B:1.3)E:0.7,C:2.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 0}, {"C", 0}});
	auto out = tree.ancestral_ml_ard(x, 1);
	REQUIRE(out.states.size() == 2);
	REQUIRE(out.rates.empty());
	for (const auto &st : out.states) {
		REQUIRE(st.probability == Approx(1.0));
	}
	REQUIRE(out.log_likelihood == Approx(0.0).margin(1e-12));
}

// Single-tip -> empty output, rate NaN.
TEST_CASE("ancestral_ml_ard returns empty for a single-tip tree", "[NewickTree][asr][ml][ard]") {
	std::vector<miint::NodeInput> nodes = {
	    {0, std::nullopt, "A", std::numeric_limits<double>::quiet_NaN(), std::nullopt}};
	auto tree = miint::NewickTree::build(nodes);
	std::unordered_map<std::string, uint32_t> x({{"A", 0}});
	auto out = tree.ancestral_ml_ard(x, 1);
	REQUIRE(out.states.empty());
	REQUIRE(std::isnan(out.rate));
}

// Rejections mirror the other models.
TEST_CASE("ancestral_ml_ard rejects bad inputs", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> missing({{"A", 0}, {"B", 1}, {"C", 2}});
	REQUIRE_THROWS_WITH(tree.ancestral_ml_ard(missing, 3), ContainsSubstring("no state for tip"));

	auto ztip = miint::NewickTree::parse("(((A:0.0,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}});
	REQUIRE_THROWS_WITH(ztip.ancestral_ml_ard(y, 3), ContainsSubstring("positive"));
}

// ARD fits k*(k-1) rates, so its state ceiling is tighter than SYM (k > 6 errors).
TEST_CASE("ancestral_ml_ard rejects too many states", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(A:1,B:1,C:1,D:1,E:1,F:1,G:1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 3}, {"E", 4}, {"F", 5}, {"G", 6}});
	REQUIRE_THROWS_WITH(tree.ancestral_ml_ard(x, 7), ContainsSubstring("model 'ER'"));
}

// Batch matches per-trait single calls (rate is NaN for ARD; compare rates/logL/posteriors).
TEST_CASE("ancestral_ml_ard batch matches per-trait single calls", "[NewickTree][asr][ml][ard]") {
	auto tree = miint::NewickTree::parse("(((A:0.5,B:1.3)E:0.7,C:2.1)F:0.6,D:1.1)root;");
	std::unordered_map<std::string, uint32_t> x({{"A", 0}, {"B", 1}, {"C", 0}, {"D", 1}}); // k=2
	std::unordered_map<std::string, uint32_t> y({{"A", 0}, {"B", 1}, {"C", 2}, {"D", 0}}); // k=3

	auto batch = tree.ancestral_ml_ard(std::vector<std::unordered_map<std::string, uint32_t>> {x, y},
	                                   std::vector<uint32_t> {2, 3});
	REQUIRE(batch.size() == 2);
	auto sx = tree.ancestral_ml_ard(x, 2);
	auto sy = tree.ancestral_ml_ard(y, 3);
	REQUIRE(batch[0].log_likelihood == Approx(sx.log_likelihood).margin(1e-9));
	REQUIRE(batch[1].log_likelihood == Approx(sy.log_likelihood).margin(1e-9));
	REQUIRE(batch[1].rates.size() == 6);
	auto mb1 = post_by_node_state(batch[1]), msy = post_by_node_state(sy);
	for (const auto &[key, val] : msy) {
		REQUIRE(mb1.at(key) == Approx(val).margin(1e-9));
	}
}
