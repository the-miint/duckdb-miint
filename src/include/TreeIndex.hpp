#pragma once

#include <cstdint>
#include <vector>

namespace miint {

class NewickTree;

// A read-only index over a rooted tree for ancestor, LCA and root-distance queries.
//
// Every subtree is a contiguous range of preorder positions, so "is a an ancestor of v"
// is one interval test, and the lowest ancestor of a node whose subtree contains some
// preorder position is found by climbing with skew-binary jump pointers (one pointer per
// node) in O(log depth). Memory is O(n): about 28 bytes per node, which matters on trees
// of tens of millions of nodes where O(n log n) structures (binary lifting, an Euler-tour
// sparse table) would not fit.
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
	// connected to the root (a cycle), or when a non-root edge's length is NaN (NULL),
	// infinite or negative. Root distances are sums of branch lengths, and a search that
	// prunes by distance needs them to grow monotonically away from the root.
	//
	// `node_ids` (optional, same size) maps each dense index to the caller's identifier,
	// used only to name nodes in error messages.
	static TreeIndex Build(const std::vector<uint32_t> &parent, const std::vector<double> &branch_length,
	                       const std::vector<int64_t> *node_ids = nullptr);

	// Build from a parsed or built tree, with the same validation.
	static TreeIndex FromNewickTree(const NewickTree &tree, const std::vector<int64_t> *node_ids = nullptr);

	size_t size() const {
		return parent_.size();
	}
	uint32_t root() const {
		return root_;
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
	// D[v] = D[parent] + branch_length[v].
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

private:
	std::vector<uint32_t> parent_;
	std::vector<uint32_t> pre_;
	std::vector<uint32_t> last_;
	std::vector<uint32_t> depth_;
	std::vector<uint32_t> jump_;
	std::vector<double> dist_;
	uint32_t root_ = NO_PARENT;
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
	uint32_t find(int64_t node_id) const;
};

} // namespace miint
