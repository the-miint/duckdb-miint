#include "TreeIndex.hpp"

#include "NewickTree.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace miint {

static_assert(TreeIndex::NO_PARENT == NewickTree::NO_PARENT, "both trees mark the root the same way");

namespace {

std::string NodeLabel(uint32_t v, const std::vector<int64_t> *node_ids) {
	if (node_ids) {
		return "node_index " + std::to_string((*node_ids)[v]);
	}
	return "node " + std::to_string(v);
}

} // namespace

TreeIndex TreeIndex::Build(const std::vector<uint32_t> &parent, const std::vector<double> &branch_length,
                           const std::vector<int64_t> *node_ids) {
	const size_t n_size = parent.size();
	if (n_size == 0) {
		throw std::invalid_argument("TreeIndex: the tree is empty");
	}
	if (branch_length.size() != n_size || (node_ids && node_ids->size() != n_size)) {
		throw std::invalid_argument("TreeIndex: the parent, branch length and node id arrays differ in size");
	}
	if (n_size >= NO_PARENT) {
		throw std::invalid_argument("TreeIndex: too many nodes (" + std::to_string(n_size) + ")");
	}
	const uint32_t n = static_cast<uint32_t>(n_size);

	// Exactly one root, and every other parent in range.
	uint32_t root = NO_PARENT;
	size_t n_roots = 0;
	for (uint32_t v = 0; v < n; ++v) {
		if (parent[v] == NO_PARENT) {
			root = v;
			++n_roots;
		} else if (parent[v] >= n) {
			throw std::invalid_argument("TreeIndex: the parent of " + NodeLabel(v, node_ids) + " is out of range");
		}
	}
	if (n_roots == 0) {
		throw std::invalid_argument("TreeIndex: the tree has no root (every node has a parent)");
	}
	if (n_roots > 1) {
		throw std::invalid_argument("TreeIndex: the tree has " + std::to_string(n_roots) +
		                            " roots; a tree has exactly one");
	}

	// Children in compressed form, each node's children in ascending index.
	std::vector<uint32_t> child_start(static_cast<size_t>(n) + 1, 0);
	for (uint32_t v = 0; v < n; ++v) {
		if (parent[v] != NO_PARENT) {
			++child_start[parent[v] + 1];
		}
	}
	for (uint32_t v = 0; v < n; ++v) {
		child_start[v + 1] += child_start[v];
	}
	std::vector<uint32_t> children(n - 1);
	{
		std::vector<uint32_t> fill(child_start.begin(), child_start.end() - 1);
		for (uint32_t v = 0; v < n; ++v) {
			if (parent[v] != NO_PARENT) {
				children[fill[parent[v]]++] = v;
			}
		}
	}

	TreeIndex idx;
	idx.root_ = root;
	idx.parent_ = parent;
	idx.pre_.assign(n, NO_PARENT);
	idx.last_.assign(n, 0);
	idx.depth_.assign(n, 0);
	idx.jump_.assign(n, root);
	idx.dist_.assign(n, 0.0);

	// Iterative preorder from the root. Pushing children in reverse visits them in
	// ascending index. A node never reached hangs off a cycle.
	std::vector<uint32_t> order;
	order.reserve(n);
	{
		std::vector<uint32_t> stack;
		stack.push_back(root);
		while (!stack.empty()) {
			uint32_t v = stack.back();
			stack.pop_back();
			idx.pre_[v] = static_cast<uint32_t>(order.size());
			order.push_back(v);
			for (uint32_t c = child_start[v + 1]; c > child_start[v]; --c) {
				stack.push_back(children[c - 1]);
			}
		}
	}
	if (order.size() != n) {
		for (uint32_t v = 0; v < n; ++v) {
			if (idx.pre_[v] == NO_PARENT) {
				throw std::invalid_argument("TreeIndex: " + NodeLabel(v, node_ids) +
				                            " is not connected to the root (the parent links form a cycle)");
			}
		}
	}

	// Branch lengths below the root: present, finite and non-negative.
	for (uint32_t v = 0; v < n; ++v) {
		if (v == root) {
			continue;
		}
		const double bl = branch_length[v];
		if (std::isnan(bl)) {
			throw std::invalid_argument("TreeIndex: the branch length of " + NodeLabel(v, node_ids) + " is NULL");
		}
		if (!std::isfinite(bl)) {
			throw std::invalid_argument("TreeIndex: the branch length of " + NodeLabel(v, node_ids) + " is not finite");
		}
		if (bl < 0) {
			throw std::invalid_argument("TreeIndex: the branch length of " + NodeLabel(v, node_ids) + " is negative (" +
			                            std::to_string(bl) + ")");
		}
	}

	// Parents before children: depth, root distance and the skew-binary jump pointer.
	// jump[v] is either v's parent or the parent's jump's jump, chosen so that following
	// jumps and parents reaches any ancestor in O(log depth) steps (Myers 1983).
	for (uint32_t i = 0; i < n; ++i) {
		const uint32_t v = order[i];
		idx.last_[v] = i;
		if (v == root) {
			continue;
		}
		const uint32_t p = parent[v];
		idx.depth_[v] = idx.depth_[p] + 1;
		idx.dist_[v] = idx.dist_[p] + branch_length[v];
		const uint32_t j1 = idx.jump_[p];
		const uint32_t j2 = idx.jump_[j1];
		idx.jump_[v] = (idx.depth_[p] - idx.depth_[j1] == idx.depth_[j1] - idx.depth_[j2]) ? j2 : p;
	}
	// Children after parents, in reverse: a subtree's last preorder position is the
	// largest among its nodes.
	for (uint32_t i = n; i-- > 0;) {
		const uint32_t v = order[i];
		if (v != root && idx.last_[v] > idx.last_[parent[v]]) {
			idx.last_[parent[v]] = idx.last_[v];
		}
	}
	return idx;
}

TreeIndex TreeIndex::FromNewickTree(const NewickTree &tree, const std::vector<int64_t> *node_ids) {
	const size_t n = tree.num_nodes();
	std::vector<uint32_t> parent(n);
	std::vector<double> branch_length(n);
	for (size_t v = 0; v < n; ++v) {
		parent[v] = tree.parent(static_cast<uint32_t>(v));
		branch_length[v] = tree.branch_length(static_cast<uint32_t>(v));
	}
	return Build(parent, branch_length, node_ids);
}

uint32_t TreeIndex::lowest_ancestor_containing(uint32_t from, uint32_t x) const {
	if (from >= parent_.size() || x >= parent_.size()) {
		throw std::out_of_range("TreeIndex::lowest_ancestor_containing: node or position out of range");
	}
	// The root contains every position, so the climb always ends. A jump is taken only
	// when its target still does not contain x, so it never overshoots the answer.
	uint32_t v = from;
	while (!contains(v, x)) {
		const uint32_t j = jump_[v];
		v = contains(j, x) ? parent_[v] : j;
	}
	return v;
}

} // namespace miint
