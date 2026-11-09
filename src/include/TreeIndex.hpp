#pragma once

#include <cstdint>
#include <optional>
#include <vector>

namespace miint {

// What TreeIndex::Build accepts as a branch length below the root.
enum class LengthPolicy {
	// Every length is known, finite and non-negative, so root distances grow away from the
	// root, as a search that prunes by distance needs.
	Searchable,
	// A NaN length (a NULL in a table) is unknown, and a distance whose path crosses it is
	// unknown too. Negative lengths are allowed: a sum of them is still a distance.
	UnknownAllowed,
};

// A read-only index over a rooted tree for ancestor, LCA and root-distance queries.
//
// Every subtree is a contiguous range of preorder positions, so "is a an ancestor of v"
// is one interval test, and the lowest ancestor of a node whose subtree contains some
// preorder position is found by climbing with skew-binary jump pointers (one pointer per
// node) in O(log depth). Memory is O(n): about 28 bytes per node, 32 when some length is
// unknown. On trees of tens of
// millions of nodes that matters: O(n log n) structures would take several times as much
// (binary lifting about 4 times, an Euler-tour sparse table about 7).
//
// Built iteratively, so arbitrarily deep trees cannot overflow the stack. The input's
// node numbering is not required to be in any order.
class TreeIndex {
public:
	static constexpr uint32_t NO_PARENT = UINT32_MAX;

	// Build from dense arrays: parent[v] is v's parent (NO_PARENT for the one root) and
	// branch_length[v] is the length of the edge above v. The root's length is ignored.
	//
	// Throws std::invalid_argument when the arrays differ in size or are empty, when there
	// is not exactly one root, when a parent is out of range, when some node is not
	// connected to the root (a cycle), when a non-root edge's length is infinite, when a root
	// distance is outside a double's range, or when a length breaks `policy` (see LengthPolicy).
	//
	// `node_ids` (optional, same size) maps each dense index to the caller's identifier,
	// used only to name nodes in error messages.
	static TreeIndex Build(const std::vector<uint32_t> &parent, const std::vector<double> &branch_length,
	                       const std::vector<int64_t> *node_ids = nullptr,
	                       LengthPolicy policy = LengthPolicy::Searchable);

	size_t size() const {
		return parent_.size();
	}
	uint32_t root() const {
		return root_;
	}
	// The policy the index was built under.
	LengthPolicy length_policy() const {
		return policy_;
	}
	uint32_t parent(uint32_t v) const {
		return parent_[v];
	}

	// Preorder position of v, and the last preorder position inside v's subtree.
	uint32_t pre(uint32_t v) const {
		return pre_[v];
	}
	uint32_t last(uint32_t v) const {
		return last_[v];
	}

	// Edges between v and the root.
	uint32_t depth(uint32_t v) const {
		return depth_[v];
	}

	// Sum of branch lengths from the root down to v (0 at the root), accumulated root first:
	// D[v] = D[parent] + branch_length[v], with an unknown length counted as 0. Where a
	// length may be unknown, distance_to_ancestor says whether a distance is known.
	double root_distance(uint32_t v) const {
		return dist_[v];
	}

	// True when preorder position x lies in a's subtree (a itself included).
	bool contains(uint32_t a, uint32_t x) const {
		return pre_[a] <= x && x <= last_[a];
	}

	// True when a is v or an ancestor of v.
	bool is_ancestor_or_self(uint32_t a, uint32_t v) const {
		return contains(a, pre_[v]);
	}

	// True when v has no children: its subtree is v alone.
	bool is_tip(uint32_t v) const {
		return pre_[v] == last_[v];
	}

	// The lowest ancestor-or-self of `from` whose subtree contains preorder position x.
	uint32_t lowest_ancestor_containing(uint32_t from, uint32_t x) const;

	// The lowest common ancestor of a and b.
	uint32_t lca(uint32_t a, uint32_t b) const {
		return lowest_ancestor_containing(a, pre_[b]);
	}

	// The length of the path from v up to a, which is v or an ancestor of v: D[v] - D[a].
	// nullopt when an edge on it has an unknown length. Throws std::out_of_range for a node
	// out of range, std::invalid_argument when a is neither v nor an ancestor of v, and
	// std::overflow_error when the length is outside a double's range (root distances at
	// opposite ends of it, which negative lengths allow).
	std::optional<double> distance_to_ancestor(uint32_t v, uint32_t a) const;

	// The length of the path between a and b, whose lowest common ancestor is lca:
	// (D[a] - D[lca]) + (D[b] - D[lca]), the expression the nearest-tip search scores with,
	// so both give the same double for the same pair, in either order. nullopt when an edge
	// on the path has an unknown length. Throws as distance_to_ancestor does, and
	// std::overflow_error when the sum of the two halves is outside a double's range.
	std::optional<double> path_length(uint32_t a, uint32_t b, uint32_t lca) const;

private:
	std::vector<uint32_t> parent_;
	std::vector<uint32_t> pre_;
	std::vector<uint32_t> last_;
	std::vector<uint32_t> depth_;
	std::vector<uint32_t> jump_;
	std::vector<double> dist_;
	// Edges of unknown length between the root and each node, held only when there are any.
	std::vector<uint32_t> unknown_above_;
	uint32_t root_ = NO_PARENT;
	LengthPolicy policy_ = LengthPolicy::Searchable;
};

// A tree's rows as a table holds them, numbered densely for TreeIndex::Build.
//
// Rows may come in any order and their ids may be any distinct values: dense index i is
// the row with the i-th smallest id, so the dense numbering keeps the caller's order.
struct DenseTree {
	static constexpr uint32_t NOT_FOUND = UINT32_MAX;

	std::vector<int64_t> node_ids; // each dense index's id, ascending
	std::vector<uint32_t> parent;  // dense; TreeIndex::NO_PARENT where the row has no parent
	std::vector<double> branch_length;

	// One entry per row: its id, its parent's id (ignored where has_parent is false) and
	// the length of the edge above it. Throws std::invalid_argument when the columns differ
	// in size, when there are no rows, when an id repeats, or when a parent id is no row's id.
	static DenseTree FromRows(std::vector<int64_t> node_ids, std::vector<int64_t> parent_ids,
	                          std::vector<bool> has_parent, std::vector<double> branch_length);

	// The dense index of an id, or NOT_FOUND when no row has it. O(1) when the ids are
	// contiguous, as read_newick's are; a binary search otherwise.
	uint32_t find(int64_t node_id) const {
		return Find(node_ids, node_id);
	}
	// The same over a DenseTree's node_ids, kept after the rest of it is released.
	static uint32_t Find(const std::vector<int64_t> &node_ids, int64_t node_id);
};

} // namespace miint
