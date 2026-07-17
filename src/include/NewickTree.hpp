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
	std::vector<IndependentContrast>
	independent_contrasts(const std::unordered_map<std::string, double> &trait_values) const;

	// Batch overload for many traits over the same tree: the trait-independent work
	// (structural + branch-length validation, variance-extended branch lengths, and
	// per-node contrast variances) is done once and reused across every trait.
	// Returns one contrast vector per input trait, in the same order. Per-trait
	// completeness is still validated; the same exceptions are thrown.
	std::vector<std::vector<IndependentContrast>>
	independent_contrasts(const std::vector<std::unordered_map<std::string, double>> &trait_values_list) const;

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
