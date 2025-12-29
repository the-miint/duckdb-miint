#pragma once

#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <vector>

namespace miint {

// Forward declarations
class NewickTree;
class NewickParser;

// Node data stored in contiguous array for cache efficiency
struct NewickNode {
	std::string name;                    // Node label (empty for unlabeled internal nodes)
	double branch_length;                // Branch length (NaN if not specified)
	std::optional<int64_t> edge_id;      // Edge identifier from {n} syntax (jplace format)
	uint32_t parent;                     // Index of parent (UINT32_MAX for root)
	std::vector<uint32_t> children;      // Indices of children (empty for tips)

	NewickNode() : branch_length(std::numeric_limits<double>::quiet_NaN()), parent(UINT32_MAX) {
	}

	NewickNode(std::string name_, double branch_length_, std::optional<int64_t> edge_id_)
	    : name(std::move(name_)), branch_length(branch_length_), edge_id(edge_id_), parent(UINT32_MAX) {
	}
};

class NewickTree {
public:
	// Constants
	static constexpr uint32_t NO_PARENT = UINT32_MAX;

	// Parse a Newick string into a tree
	// Throws on parse errors (strict mode)
	static NewickTree parse(std::string_view newick);

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
