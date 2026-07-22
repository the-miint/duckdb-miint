#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace miint {

// Forward declarations
class NewickTree;
class NewickParser;

// Placement data for insert_fully_resolved
// Represents where a query sequence should be placed on a reference tree
struct Placement {
	std::string fragment_id;  // Name of the query sequence
	int64_t edge_id;          // Edge to place on (from jplace)
	double distal_length;     // Distance from child end of edge
	double pendant_length;    // Branch length of the new fragment
	double like_weight_ratio; // Quality score (higher is better, for deduplication)
};

// Node data stored in contiguous array for cache efficiency
struct NewickNode {
	std::string name;               // Node label (empty for unlabeled internal nodes)
	double branch_length;           // Branch length (NaN if not specified)
	std::optional<int64_t> edge_id; // Edge identifier from {n} syntax (jplace format)
	uint32_t parent;                // Index of parent (UINT32_MAX for root)
	std::vector<uint32_t> children; // Indices of children (empty for tips)

	NewickNode() : branch_length(std::numeric_limits<double>::quiet_NaN()), parent(UINT32_MAX) {
	}

	NewickNode(std::string name_, double branch_length_, std::optional<int64_t> edge_id_)
	    : name(std::move(name_)), branch_length(branch_length_), edge_id(edge_id_), parent(UINT32_MAX) {
	}
};

// Input data for building a tree programmatically
struct NodeInput {
	int64_t node_id;                  // Unique identifier for this node
	std::optional<int64_t> parent_id; // Parent's node_id (nullopt for root)
	std::string name;                 // Node label (may be empty)
	double branch_length;             // Branch length (NaN if not specified)
	std::optional<int64_t> edge_id;   // Edge identifier (nullopt if not specified)
};

// One standardized phylogenetic independent contrast (Felsenstein 1985) at an
// internal node, for a single trait.
struct IndependentContrast {
	uint32_t node;             // Internal node index (this tree's index)
	double contrast;           // Standardized contrast (X_i - X_j) / sqrt(v_i + v_j)
	double ancestral_estimate; // Reconstructed trait value X_k at the node
	double contrast_variance;  // v_i + v_j (the contrast's expected variance)
};

// One reconstructed ancestral state under a Brownian-motion model, at an internal
// node, for a single trait. Unlike IndependentContrast's `ancestral_estimate`
// (the subtree-only down-pass value, exact only at the root), `estimate` here is
// the full marginal reconstruction from a post-order down-pass plus a pre-order
// up-pass (Felsenstein 1985; Schluter et al. 1997).
struct AncestralStateBM {
	uint32_t node;   // Internal node index (this tree's dense index)
	double estimate; // Marginal ML ancestral estimate (mean) at the node
	double variance; // Variance of the estimate (sigma^2_hat * structural variance)
	double ci_low;   // 95% CI lower bound: estimate - 1.96*sqrt(variance)
	double ci_high;  // 95% CI upper bound: estimate + 1.96*sqrt(variance)
};

// One (internal node, state) entry of a discrete-trait parsimony reconstruction
// (Sankoff 1975; Fitch 1971 is the unit-cost special case). `min_cost` is the
// minimum total tree cost achievable with this node fixed to this state; a state
// is in the most-parsimonious reconstruction (MPR) set iff min_cost equals the
// whole-tree parsimony score (`in_mpr`). One entry is emitted per (internal node,
// state) so the full per-node cost profile and any ties are first-class.
struct AncestralStateParsimony {
	uint32_t node;   // Internal node index (this tree's dense index)
	uint32_t state;  // State index (0..k-1) into the trait's sorted state alphabet
	bool in_mpr;     // Whether this state is in the node's MPR set (min_cost == score)
	double min_cost; // Min total tree cost with this node fixed to this state
};

// One (internal node, state) marginal posterior of a discrete-trait maximum-
// likelihood reconstruction under the Mk model (Felsenstein pruning + a marginal
// two-pass; Yang, Kumar & Nei 1995). `probability` is P(node = state | data) and the
// per-node probabilities over states sum to 1.
struct AncestralStateML {
	uint32_t node;      // Internal node index (this tree's dense index)
	uint32_t state;     // State index (0..k-1) into the trait's sorted state alphabet
	double probability; // Marginal posterior P(node = state | data)
};

// Result of one Mk-ML reconstruction for a single trait: the per-(internal node,
// state) posteriors plus the fitted (or fixed) rate and the model log-likelihood at
// that rate.
struct AncestralMLResult {
	std::vector<AncestralStateML> states;
	double rate;               // ER: fitted (ML) or user-fixed scalar rate. SYM/ARD: NaN (a
	                           // rate matrix has no single scalar; see `rates`).
	double log_likelihood;     // Model log-likelihood at the fitted/fixed parameters
	std::vector<double> rates; // Fitted off-diagonal rates; empty for ER. SYM: k(k-1)/2 in
	                           // pair order (i,j), i<j. ARD: k(k-1) in row order (i,j), i!=j.
};

class NewickTree {
public:
	// Constants
	static constexpr uint32_t NO_PARENT = UINT32_MAX;
	static constexpr uint32_t MAX_NODES = UINT32_MAX - 1;

	// Parse a Newick string into a tree
	// Throws on parse errors (strict mode)
	static NewickTree parse(std::string_view newick);

	// Build a tree from node input data
	// Validates:
	// - Exactly one root (node with nullopt parent_id)
	// - All parent references are valid
	// - No cycles
	// Throws on validation errors
	static NewickTree build(const std::vector<NodeInput> &nodes);

	// Serialize tree back to Newick format
	// Includes edge identifiers if present
	std::string to_newick() const;

	// ========================================================================
	// Navigation
	// ========================================================================

	// Get the root node index
	uint32_t root() const {
		return root_;
	}

	// Get parent of a node (returns NO_PARENT for root)
	uint32_t parent(uint32_t node) const {
		return nodes_[node].parent;
	}

	// Get children of a node (empty vector for tips)
	const std::vector<uint32_t> &children(uint32_t node) const {
		return nodes_[node].children;
	}

	// Check if node is a tip (leaf)
	bool is_tip(uint32_t node) const {
		return nodes_[node].children.empty();
	}

	// ========================================================================
	// Node properties
	// ========================================================================

	// Get node name (may be empty for internal nodes)
	const std::string &name(uint32_t node) const {
		return nodes_[node].name;
	}

	// Get branch length (NaN if not specified)
	double branch_length(uint32_t node) const {
		return nodes_[node].branch_length;
	}

	// Get edge identifier (nullopt if not specified)
	std::optional<int64_t> edge_id(uint32_t node) const {
		return nodes_[node].edge_id;
	}

	// ========================================================================
	// Tree statistics
	// ========================================================================

	// Total number of nodes
	size_t num_nodes() const {
		return nodes_.size();
	}

	// Number of tip (leaf) nodes
	size_t num_tips() const;

	// ========================================================================
	// Traversal
	// ========================================================================

	// Post-order traversal (children before parents)
	std::vector<uint32_t> postorder() const;

	// Pre-order traversal (parents before children)
	std::vector<uint32_t> preorder() const;

	// ========================================================================
	// Queries
	// ========================================================================

	// Find node by name (returns nullopt if not found)
	// Note: O(n) scan, consider building index for frequent lookups
	std::optional<uint32_t> find_node_by_name(std::string_view name) const;

	// Find node by edge ID (returns nullopt if not found)
	// Note: O(n) scan, consider building index for frequent lookups
	std::optional<uint32_t> find_node_by_edge_id(int64_t edge_id) const;

	// Build a map from edge_id to node index for O(1) lookups
	// Call once before repeated find_node_by_edge_id calls
	std::unordered_map<int64_t, uint32_t> build_edge_index() const;

	// ========================================================================
	// Distance calculations
	// ========================================================================

	// Get all tip node indices
	std::vector<uint32_t> tips() const;

	// Get all tip names
	std::vector<std::string> tip_names() const;

	// Calculate distance from a node to the root (sum of branch lengths)
	// NaN branch lengths are treated as 0
	double distance_to_root(uint32_t node) const;

	// Find the lowest common ancestor of two nodes
	uint32_t find_lca(uint32_t a, uint32_t b) const;

	// Calculate pairwise distance between two nodes via their LCA
	// NaN branch lengths are treated as 0
	double pairwise_distance(uint32_t a, uint32_t b) const;

	// ========================================================================
	// Phylogenetic placement
	// ========================================================================

	// Insert query sequences into the tree at their placement positions
	// Creates a fully resolved tree where each placement gets its own edge
	//
	// Algorithm:
	// 1. Deduplicates placements by fragment_id (keeps highest like_weight_ratio,
	//    then lowest pendant_length)
	// 2. Groups placements by edge_id
	// 3. Sorts each group by distal_length descending
	// 4. For each edge, creates a chain of new internal nodes with fragments
	//
	// Throws if edge_id not found or distal_length exceeds edge length
	void insert_fully_resolved(const std::vector<Placement> &placements);

	// ========================================================================
	// Subsetting (shear / prune to a set of tips)
	// ========================================================================

	// Return a new tree restricted to the tips named in `keep_names`.
	//
	// - collapse=false: keep every ancestor of a kept tip. Original parent
	//   links, branch lengths, names, and edge ids are preserved; only nodes
	//   that lie on no kept root-path are dropped. Internal nodes left with a
	//   single child are retained (unifurcations preserved).
	// - collapse=true: additionally remove single-child internal nodes, summing
	//   their branch lengths onto the surviving descendant edge. The lowest
	//   common ancestor of the kept tips becomes the new root (nodes above it
	//   are dropped). A collapsed edge keeps the surviving (lower) node's
	//   edge_id; merged intermediate edge_ids are dropped. Branch-length
	//   summation treats NaN (unspecified) as 0, but a chain that is entirely
	//   NaN stays NaN (topology-only trees are preserved).
	//
	// `keep_names` is matched against tip names only (internal-node labels do
	// not count). A name matching no tip is "missing": if ignore_missing is
	// false, throws std::runtime_error listing the missing names; otherwise the
	// missing names are skipped. Throws std::runtime_error if no tip matches.
	//
	// The returned tree's node indices are reassigned 0-based (as with build()).
	NewickTree shear(const std::unordered_set<std::string> &keep_names, bool collapse, bool ignore_missing) const;

	// ========================================================================
	// Resolution (make multifurcations bifurcating)
	// ========================================================================

	// Return a new tree in which every node with more than two children is
	// resolved into a series of bifurcations. For a node with children
	// c0,c1,...,c(m-1) (m >= 3) in their existing order, the resolution is a
	// deterministic left-comb: the node keeps c0 and a new internal node N1;
	// N1 keeps c1 and N2; ...; N(m-2) keeps c(m-2) and c(m-1). The m-2 new
	// internal nodes are unnamed, have branch length 0, and no edge id, so tip
	// counts and every original edge length are preserved and root-to-tip
	// distances are unchanged.
	//
	// Nodes with two or fewer children are left untouched, including
	// single-child unifurcations (this resolves polytomies only; use shear() to
	// collapse unifurcations). The returned tree's node indices are reassigned
	// 0-based (original nodes keep their index; new nodes follow).
	NewickTree resolve_multifurcations() const;

	// ========================================================================
	// Comparative methods
	// ========================================================================

	// Compute Felsenstein (1985) phylogenetic independent contrasts for one
	// numeric per-tip trait, keyed by tip name. Returns one contrast per internal
	// node (n-1 for a rooted bifurcating tree with n tips).
	//
	// Requires a strictly bifurcating tree (every internal node has exactly two
	// children) with unique tip names, a trait value for exactly the set of tips
	// (no missing tips, no extra names), and finite non-negative branch lengths on
	// every non-root edge; each contrast's variance (v_i + v_j) must be > 0
	// (so at most one of a node's two children may have zero effective length).
	// Zero-length internal edges are allowed. Throws std::runtime_error /
	// std::invalid_argument describing the first violation. The root's own branch
	// length is unused.
	//
	// `node_ids` (optional) maps this tree's dense node index to a caller-facing
	// identifier (e.g. the original node_index of the source table); when supplied,
	// error messages for unnamed nodes report that identifier instead of the
	// internal dense index. Pass nullptr to use the dense index.
	std::vector<IndependentContrast> independent_contrasts(const std::unordered_map<std::string, double> &trait_values,
	                                                       const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree: the trait-independent work
	// (structural + branch-length validation, variance-extended branch lengths, and
	// per-node contrast variances) is done once and reused across every trait.
	// Returns one contrast vector per input trait, in the same order. Per-trait
	// completeness is still validated; the same exceptions are thrown. `node_ids` is
	// as above (used only for error-message identifiers).
	std::vector<std::vector<IndependentContrast>>
	independent_contrasts(const std::vector<std::unordered_map<std::string, double>> &trait_values_list,
	                      const std::vector<int64_t> *node_ids = nullptr) const;

	// Reconstruct continuous ancestral states under a Brownian-motion model for one
	// numeric per-tip trait, keyed by tip name. Returns one AncestralStateBM per
	// internal node (including the root); tips are observed and not returned.
	//
	// Supports multifurcations — the precision-weighted Gaussian message passing is
	// arity-agnostic (unlike PIC's pairwise contrasts). Requires unique tip names, a
	// trait value for exactly the tip set, finite non-negative branch lengths on
	// every non-root edge, and a strictly positive branch length on every tip edge
	// (a zero-length tip edge would pin the estimate and is rejected here). The BM
	// rate is estimated by REML (the unbiased estimator: sum of squared contrasts /
	// (n_tips - 1)); this is the PIC-consistent choice and differs from full ML
	// (divisor n_tips, e.g. ape::ace method="ML") by a factor (n_tips - 1)/n_tips.
	// `variance` and the 95% CI are that rate times the structural variance. Throws
	// std::runtime_error on any violation.
	//
	// Every internal node must have >= 2 children: a unifurcation (1 child) is
	// rejected with a shear_tree hint (it carries no information). `node_ids` is as in
	// independent_contrasts (error-message identifiers only).
	std::vector<AncestralStateBM> ancestral_states_bm(const std::unordered_map<std::string, double> &trait_values,
	                                                  const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree: the trait-independent work
	// (structural + branch-length validation and the per-node structural variances)
	// is done once and reused across every trait. Returns one vector per input trait,
	// in the same order. Per-trait completeness + finiteness is still validated; the
	// same exceptions are thrown. `node_ids` is as above.
	std::vector<std::vector<AncestralStateBM>>
	ancestral_states_bm(const std::vector<std::unordered_map<std::string, double>> &trait_values_list,
	                    const std::vector<int64_t> *node_ids = nullptr) const;

	// Reconstruct discrete ancestral states under Sankoff (1975) parsimony for one
	// categorical per-tip trait, keyed by tip name. `tip_states` maps each tip name to
	// a state index in 0..k-1; `cost` is the k*k row-major substitution cost matrix
	// (cost[from*k + to], the cost of an edge whose parent end is `from` and child end
	// is `to`; the unit matrix — 0 on the diagonal, 1 elsewhere — gives Fitch 1971
	// parsimony). Returns one AncestralStateParsimony per (internal node, state): the
	// minimum total tree cost with that node fixed to that state, and whether that
	// state is in the node's most-parsimonious reconstruction (MPR) set. Tips are
	// observed and not returned.
	//
	// Branch lengths are ignored (topology-only trees are fine). Supports any arity:
	// multifurcations (the DP sums over an arbitrary child set) AND unifurcations (a
	// single-child internal node is a valid degree-2 ancestor here, unlike BM/PIC —
	// parsimony has no numerical obstruction to it). Requires unique tip names, a
	// state for exactly the tip set (no missing tips, no extra names), every tip state
	// index < k, and a k*k cost matrix with finite non-negative entries. Throws
	// std::runtime_error / std::invalid_argument on any violation. `node_ids` is as in
	// independent_contrasts (error-message identifiers only).
	std::vector<AncestralStateParsimony>
	ancestral_parsimony(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
	                    const std::vector<double> &cost, const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree: the trait-independent work
	// (structural validation and the post-/pre-order traversals) is done once and
	// reused across every trait. Each trait carries its own alphabet, so `k_list` and
	// `cost_list` are per-trait and parallel to `tip_states_list` (the unit-cost
	// default uses a per-trait alphabet of that trait's observed states; a shared cost
	// matrix repeats the same k/cost). Returns one vector per input trait, in order.
	// Per-trait completeness/validity is still checked; the same exceptions are thrown.
	std::vector<std::vector<AncestralStateParsimony>>
	ancestral_parsimony(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
	                    const std::vector<uint32_t> &k_list, const std::vector<std::vector<double>> &cost_list,
	                    const std::vector<int64_t> *node_ids = nullptr) const;

	// Reconstruct discrete ancestral states by maximum likelihood under the Mk model
	// with EQUAL rates (the ER model; Lewis 2001, a k-state generalization of Jukes-
	// Cantor 1969) for one categorical per-tip trait. `tip_states` maps each tip name
	// to a state index in 0..k-1. The single rate is fitted by maximum likelihood
	// (golden-section search on the log-likelihood) unless `fixed_rate` is given
	// (> 0), in which case that rate is used directly. Returns the per-(internal node,
	// state) marginal posteriors, the rate, and the log-likelihood. Tips are observed
	// and not returned.
	//
	// Uses Felsenstein's (1981) pruning algorithm with per-node rescaling for the
	// likelihood (scalers retained in logL), a uniform root prior (1/k, the ER
	// stationary distribution), and a marginal up-pass (Yang, Kumar & Nei 1995) for
	// posteriors. For ER, P(t) is closed-form (no eigendecomposition):
	// P(t)[i][j] = 1/k + (delta_ij - 1/k) * exp(-k * rate * t). Supports any arity
	// (multifurcations and unifurcations). Requires unique tip names, a state for
	// exactly the tip set (every index < k), finite non-negative branch lengths on
	// every non-root edge, and a strictly positive branch length on every tip edge (a
	// zero-length tip edge would put a hard zero in the leave-one-out and is rejected,
	// as in ancestral_states_bm). Zero-length internal edges are allowed. Throws
	// std::runtime_error on any violation. `node_ids` is as in independent_contrasts.
	AncestralMLResult ancestral_ml(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
	                               std::optional<double> fixed_rate = std::nullopt,
	                               const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree: the trait-independent work
	// (structural + branch-length validation and the traversals) is done once and
	// reused. Each trait carries its own alphabet, so `k_list` is per-trait and
	// parallel to `tip_states_list`. `fixed_rate`, if given, applies to every trait.
	// Returns one result per input trait, in order.
	std::vector<AncestralMLResult>
	ancestral_ml(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
	             const std::vector<uint32_t> &k_list, std::optional<double> fixed_rate = std::nullopt,
	             const std::vector<int64_t> *node_ids = nullptr) const;

	// Mk maximum-likelihood ASR under the SYM (symmetric-rates) model: the rate matrix Q
	// is symmetric (q_ij = q_ji), with k(k-1)/2 free off-diagonal rates, all fitted by
	// maximum likelihood (there is no single scalar rate to fix, so unlike the ER overload
	// there is no fixed-rate argument). Q symmetric => uniform stationary distribution
	// (matching the uniform root prior) and a real eigenspectrum, so P(t) = exp(Q*t) is
	// built from a SelfAdjointEigenSolver eigendecomposition. Reversible, so the
	// likelihood is invariant to root placement (pulley principle). Same pruning /
	// rescaling / marginal-posterior machinery and the same branch-length policy as the ER
	// overload. At k=2 there is exactly one free rate and SYM reduces to ER. The multi-
	// dimensional rate fit uses a Nelder-Mead simplex search (Nelder & Mead 1965),
	// warm-started from the ER fit. Returns `rate` = NaN and `rates` = the fitted
	// off-diagonal rates. Throws std::runtime_error on any input violation.
	AncestralMLResult ancestral_ml_sym(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
	                                   const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree (see the ER batch overload); the
	// trait-independent validation and traversals are computed once and reused.
	std::vector<AncestralMLResult>
	ancestral_ml_sym(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
	                 const std::vector<uint32_t> &k_list, const std::vector<int64_t> *node_ids = nullptr) const;

	// Mk maximum-likelihood ASR under the ARD (all-rates-different) model: the rate matrix Q
	// is fully general (all k*(k-1) off-diagonal rates free), so the model is NOT reversible
	// -- the reconstruction depends on root placement (no pulley principle) and the
	// stationary distribution is not uniform. A uniform root prior 1/k is used (matching the
	// ER/SYM convention and ape::ace). P(t) = exp(Q*t) is built by a general matrix
	// exponential (scaling-and-squaring Pade, Higham 2005). Same pruning / rescaling /
	// marginal-posterior machinery and branch-length policy as the other overloads. The
	// k*(k-1) rates are fitted by a Nelder-Mead simplex search (Nelder & Mead 1965),
	// warm-started from the ER fit. Returns `rate` = NaN and `rates` = the fitted off-
	// diagonal rates in row order (i,j), i!=j. ARD does not scale to large state counts
	// (k*(k-1) simplex parameters); high-k traits are rejected. Throws on any input
	// violation.
	AncestralMLResult ancestral_ml_ard(const std::unordered_map<std::string, uint32_t> &tip_states, uint32_t k,
	                                   const std::vector<int64_t> *node_ids = nullptr) const;

	// Batch overload for many traits over the same tree (see the ER batch overload).
	std::vector<AncestralMLResult>
	ancestral_ml_ard(const std::vector<std::unordered_map<std::string, uint32_t>> &tip_states_list,
	                 const std::vector<uint32_t> &k_list, const std::vector<int64_t> *node_ids = nullptr) const;

	// ========================================================================
	// Modification (for insert_fully_resolved)
	// ========================================================================

	// Add a new disconnected node, returns its index
	uint32_t add_node(const std::string &name, double branch_length, std::optional<int64_t> edge_id);

	// Set a node's parent (also updates parent's children list)
	// If node already has a parent, it is first removed from that parent
	void set_parent(uint32_t child, uint32_t parent);

	// Remove a child from a parent's children list
	// Does not update the child's parent field
	void remove_child(uint32_t parent, uint32_t child);

	// Set branch length
	void set_branch_length(uint32_t node, double length) {
		nodes_[node].branch_length = length;
	}

	// Set edge identifier
	void set_edge_id(uint32_t node, std::optional<int64_t> edge_id) {
		nodes_[node].edge_id = edge_id;
	}

	// Set node name
	void set_name(uint32_t node, const std::string &name) {
		nodes_[node].name = name;
	}

private:
	friend class NewickParser;

	std::vector<NewickNode> nodes_;
	uint32_t root_ = NO_PARENT;
};

} // namespace miint
