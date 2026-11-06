#include "TreeIndex.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

namespace miint {

namespace {

std::string NodeLabel(uint32_t v, const std::vector<int64_t> *node_ids) {
	if (node_ids) {
		return "node_index " + std::to_string((*node_ids)[v]);
	}
	return "node " + std::to_string(v);
}

template <class T>
std::vector<T> Permuted(const std::vector<T> &values, const std::vector<uint32_t> &order) {
	std::vector<T> out;
	out.reserve(order.size());
	for (uint32_t i : order) {
		out.push_back(values[i]);
	}
	return out;
}

} // namespace

TreeIndex TreeIndex::Build(const std::vector<uint32_t> &parent, const std::vector<double> &branch_length,
                           const std::vector<int64_t> *node_ids, LengthPolicy policy) {
	const size_t n_size = parent.size();
	if (n_size == 0) {
		throw std::invalid_argument("the tree is empty");
	}
	if (branch_length.size() != n_size || (node_ids && node_ids->size() != n_size)) {
		throw std::invalid_argument("the parent, branch length and node id arrays differ in size");
	}
	if (n_size >= NO_PARENT) {
		throw std::invalid_argument("too many nodes (" + std::to_string(n_size) + ")");
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
			throw std::invalid_argument("the parent of " + NodeLabel(v, node_ids) + " is out of range");
		}
	}
	if (n_roots == 0) {
		throw std::invalid_argument("the tree has no root (every node has a parent)");
	}
	if (n_roots > 1) {
		throw std::invalid_argument("the tree has " + std::to_string(n_roots) + " roots; a tree has exactly one");
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
	idx.policy_ = policy;
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
				throw std::invalid_argument(NodeLabel(v, node_ids) +
				                            " is not connected to the root (the parent links form a cycle)");
			}
		}
	}

	// Branch lengths below the root: never infinite; known and non-negative for a search.
	const bool searchable = policy == LengthPolicy::Searchable;
	bool any_unknown = false;
	for (uint32_t v = 0; v < n; ++v) {
		if (v == root) {
			continue;
		}
		const double bl = branch_length[v];
		if (std::isnan(bl)) {
			if (searchable) {
				throw std::invalid_argument("the branch length of " + NodeLabel(v, node_ids) + " is NULL or NaN");
			}
			any_unknown = true;
			continue;
		}
		if (!std::isfinite(bl)) {
			throw std::invalid_argument("the branch length of " + NodeLabel(v, node_ids) + " is not finite");
		}
		if (searchable && bl < 0) {
			throw std::invalid_argument("the branch length of " + NodeLabel(v, node_ids) + " is negative (" +
			                            std::to_string(bl) + ")");
		}
	}
	if (any_unknown) {
		idx.unknown_above_.assign(n, 0);
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
		const double bl = branch_length[v];
		if (std::isnan(bl)) {
			idx.dist_[v] = idx.dist_[p];
			idx.unknown_above_[v] = idx.unknown_above_[p] + 1;
		} else {
			idx.dist_[v] = idx.dist_[p] + bl;
			if (any_unknown) {
				idx.unknown_above_[v] = idx.unknown_above_[p];
			}
		}
		if (!std::isfinite(idx.dist_[v])) {
			throw std::invalid_argument("the root distance of " + NodeLabel(v, node_ids) +
			                            " is not finite (the branch lengths above it sum past the largest double)");
		}
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

DenseTree DenseTree::FromRows(std::vector<int64_t> node_ids, std::vector<int64_t> parent_ids,
                              std::vector<bool> has_parent, std::vector<double> branch_length) {
	const size_t n = node_ids.size();
	if (parent_ids.size() != n || has_parent.size() != n || branch_length.size() != n) {
		throw std::invalid_argument("the node, parent and branch length columns differ in size");
	}
	if (n == 0) {
		throw std::invalid_argument("the tree is empty");
	}
	if (n >= TreeIndex::NO_PARENT) {
		throw std::invalid_argument("too many nodes (" + std::to_string(n) + ")");
	}

	// Rows in ascending id. They usually are already (read_newick numbers in row order),
	// and then nothing moves.
	if (!std::is_sorted(node_ids.begin(), node_ids.end())) {
		std::vector<uint32_t> order(n);
		std::iota(order.begin(), order.end(), 0u);
		std::sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) { return node_ids[a] < node_ids[b]; });
		node_ids = Permuted(node_ids, order);
		parent_ids = Permuted(parent_ids, order);
		has_parent = Permuted(has_parent, order);
		branch_length = Permuted(branch_length, order);
	}
	for (size_t i = 1; i < n; ++i) {
		if (node_ids[i] == node_ids[i - 1]) {
			throw std::invalid_argument("node_index " + std::to_string(node_ids[i]) + " appears more than once");
		}
	}

	DenseTree out;
	out.node_ids = std::move(node_ids);
	out.branch_length = std::move(branch_length);
	out.parent.resize(n);
	for (size_t i = 0; i < n; ++i) {
		if (!has_parent[i]) {
			out.parent[i] = TreeIndex::NO_PARENT;
			continue;
		}
		const uint32_t p = out.find(parent_ids[i]);
		if (p == NOT_FOUND) {
			throw std::invalid_argument("the parent_index of node_index " + std::to_string(out.node_ids[i]) + " (" +
			                            std::to_string(parent_ids[i]) + ") is not the node_index of any row");
		}
		out.parent[i] = p;
	}
	return out;
}

uint32_t DenseTree::Find(const std::vector<int64_t> &node_ids, int64_t node_id) {
	if (node_ids.empty()) {
		return NOT_FOUND;
	}
	// Offsets in unsigned arithmetic, which is defined across the whole int64 range. The ids
	// are sorted and distinct, so they are contiguous exactly when they span n - 1.
	const uint64_t first = static_cast<uint64_t>(node_ids.front());
	if (static_cast<uint64_t>(node_ids.back()) - first == node_ids.size() - 1) {
		const uint64_t offset = static_cast<uint64_t>(node_id) - first;
		return offset < node_ids.size() ? static_cast<uint32_t>(offset) : NOT_FOUND;
	}
	auto it = std::lower_bound(node_ids.begin(), node_ids.end(), node_id);
	return (it != node_ids.end() && *it == node_id) ? static_cast<uint32_t>(it - node_ids.begin()) : NOT_FOUND;
}

std::optional<double> TreeIndex::distance_to_ancestor(uint32_t v, uint32_t a) const {
	if (v >= parent_.size() || a >= parent_.size()) {
		throw std::out_of_range("TreeIndex::distance_to_ancestor: node out of range");
	}
	if (!is_ancestor_or_self(a, v)) {
		throw std::invalid_argument("TreeIndex::distance_to_ancestor: node " + std::to_string(a) +
		                            " is not an ancestor of node " + std::to_string(v));
	}
	if (!unknown_above_.empty() && unknown_above_[v] != unknown_above_[a]) {
		return std::nullopt;
	}
	return dist_[v] - dist_[a];
}

std::optional<double> TreeIndex::path_length(uint32_t a, uint32_t b, uint32_t lca) const {
	const auto up = distance_to_ancestor(a, lca);
	const auto down = distance_to_ancestor(b, lca);
	if (!up || !down) {
		return std::nullopt;
	}
	return *up + *down;
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
