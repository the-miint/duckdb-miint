#pragma once

#include "NewickTree.hpp"
#include "TreeIndex.hpp"

#include <algorithm>
#include <cstdint>
#include <limits>
#include <numeric>
#include <random>
#include <stdexcept>
#include <vector>

namespace miint::test {

// A tree as plain arrays, the input TreeIndex::Build takes: parent[v] (NO_PARENT for the
// root) and the length of the edge above v.
struct TreeArrays {
	std::vector<uint32_t> parent;
	std::vector<double> length;
};

constexpr uint32_t kNoParent = TreeIndex::NO_PARENT;

static_assert(TreeIndex::NO_PARENT == NewickTree::NO_PARENT, "both trees mark the root the same way");

// The index of a parsed Newick tree, through the same validation as any other input.
inline TreeIndex index_of(const NewickTree &tree) {
	const uint32_t n = static_cast<uint32_t>(tree.num_nodes());
	std::vector<uint32_t> parent(n);
	std::vector<double> length(n);
	for (uint32_t v = 0; v < n; ++v) {
		parent[v] = tree.parent(v);
		length[v] = tree.branch_length(v);
	}
	return TreeIndex::Build(parent, length);
}

// ---------------------------------------------------------------------------
// A naive model of the tree, by walking parents: what the index and the searches built
// on it are checked against.
// ---------------------------------------------------------------------------

inline bool naive_is_ancestor_or_self(const TreeArrays &t, uint32_t a, uint32_t v) {
	for (uint32_t cur = v; cur != kNoParent; cur = t.parent[cur]) {
		if (cur == a) {
			return true;
		}
	}
	return false;
}

inline uint32_t naive_lca(const TreeArrays &t, uint32_t a, uint32_t b) {
	std::vector<bool> above_a(t.parent.size(), false);
	for (uint32_t cur = a; cur != kNoParent; cur = t.parent[cur]) {
		above_a[cur] = true;
	}
	for (uint32_t cur = b; cur != kNoParent; cur = t.parent[cur]) {
		if (above_a[cur]) {
			return cur;
		}
	}
	throw std::logic_error("naive_lca: no common ancestor");
}

inline uint32_t naive_depth(const TreeArrays &t, uint32_t v) {
	uint32_t d = 0;
	for (uint32_t cur = v; t.parent[cur] != kNoParent; cur = t.parent[cur]) {
		++d;
	}
	return d;
}

// The root distance, summed root first (the order TreeIndex promises).
inline double naive_root_distance(const TreeArrays &t, uint32_t v) {
	std::vector<uint32_t> path; // v up to, not including, the root
	for (uint32_t cur = v; t.parent[cur] != kNoParent; cur = t.parent[cur]) {
		path.push_back(cur);
	}
	double d = 0.0;
	for (auto it = path.rbegin(); it != path.rend(); ++it) {
		d += t.length[*it];
	}
	return d;
}

inline bool naive_is_tip(const TreeArrays &t, uint32_t v) {
	return std::find(t.parent.begin(), t.parent.end(), v) == t.parent.end();
}

// ---------------------------------------------------------------------------
// Random trees. Structure first (each node's parent is an earlier node, and often the
// immediately preceding one, so long chains appear), then the labels are shuffled so the
// numbering is in no particular order. Dyadic lengths are exact in binary and include
// many zeros, which produce ties; otherwise lengths are uniform in [0, 1).
// ---------------------------------------------------------------------------

inline TreeArrays random_tree(std::mt19937 &rng, uint32_t n, bool dyadic) {
	std::vector<uint32_t> structural_parent(n, kNoParent);
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	for (uint32_t k = 1; k < n; ++k) {
		if (unit(rng) < 0.3) {
			structural_parent[k] = k - 1;
		} else {
			structural_parent[k] = std::uniform_int_distribution<uint32_t>(0, k - 1)(rng);
		}
	}
	std::vector<uint32_t> label(n);
	std::iota(label.begin(), label.end(), 0u);
	std::shuffle(label.begin(), label.end(), rng);

	const double dyadic_lengths[] = {0.0, 0.0, 0.125, 0.25, 0.5, 1.0, 2.0};
	TreeArrays t;
	t.parent.assign(n, kNoParent);
	t.length.assign(n, 0.0);
	for (uint32_t k = 0; k < n; ++k) {
		if (structural_parent[k] != kNoParent) {
			t.parent[label[k]] = label[structural_parent[k]];
		}
		t.length[label[k]] = dyadic ? dyadic_lengths[std::uniform_int_distribution<int>(0, 6)(rng)] : unit(rng);
	}
	// The root's own length is ignored, whatever it holds.
	t.length[label[0]] = std::numeric_limits<double>::quiet_NaN();
	return t;
}

// A caterpillar `levels` deep: internal node k (index k) is the parent of internal node
// k+1 and of tip k (index levels + k); the deepest internal node also holds the last tip
// (index 2 * levels). Every edge has length 1, so root distances are exact integers:
// internal node k is at k, tip k at k + 1, and the last tip at levels.
inline TreeArrays caterpillar(uint32_t levels) {
	const uint32_t n = 2 * levels + 1;
	TreeArrays t;
	t.parent.assign(n, kNoParent);
	t.length.assign(n, 1.0);
	for (uint32_t k = 1; k < levels; ++k) {
		t.parent[k] = k - 1;
	}
	for (uint32_t k = 0; k < levels; ++k) {
		t.parent[levels + k] = k;
	}
	t.parent[2 * levels] = levels - 1;
	return t;
}

} // namespace miint::test
