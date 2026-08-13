#pragma once
//
// Pure (DuckDB-free) UPGMA hierarchical clustering (average linkage) over a
// dense symmetric distance matrix, used to score cluster recovery in the
// Kuczynski-2010 reproduction (UPGMA of community distances -> are the true
// clusters monophyletic clades?).
//
// "UPGMA" = Unweighted Pair Group Method with Arithmetic mean (Sokal & Michener
// 1958): each merge joins the two closest clusters and the merged-cluster
// distance to every other cluster is the arithmetic mean of member-pair
// distances, i.e. weighted by cluster size -- the same average-linkage
// definition as scipy.cluster.hierarchy.linkage(method='average'), which the
// parity test uses as its oracle (see THIRD_PARTY_LICENSES.md).
//
// Header/impl split with NO DuckDB dependency so the Catch2 unit-test target
// links it without libduckdb. The DuckDB wrapper (cluster_upgma_function.cpp)
// reads a condensed (sample_a, sample_b, distance) table into a dense matrix
// (reusing ReadDistanceTable), calls in here, and emits a read_newick-shaped
// tree table (node_index, name, branch_length, edge_id, parent_index, is_tip)
// so the result plugs into every miint tree consumer.

#include <cstdint>
#include <vector>

namespace miint {

//! One node of the UPGMA tree. Nodes are returned in a fixed order: tips first
//! at positions 0..n-1 (position == leaf index into the caller's sample list),
//! then the n-1 internal nodes in merge order at positions n..2n-2; the root is
//! the last node (position 2n-2). A node's array position is its node_index.
struct UpgmaNode {
	int32_t leaf_index;   //!< index into the caller's sample list for a tip; -1 for an internal node
	int32_t parent;       //!< node_index of the parent; -1 for the root
	double branch_length; //!< length of the edge to the parent (root: 0)
	double height;        //!< ultrametric node height (tips: 0; internal: merge_distance / 2)
};

//! UPGMA (average linkage) on a dense row-major n x n symmetric distance matrix
//! (zero diagonal). Returns 2n-1 nodes (see UpgmaNode ordering). Ties in the
//! minimum-distance search break toward the lowest (i, j) index pair, making the
//! merge order deterministic.
//!
//! Branch lengths are height differences; UPGMA on non-ultrametric input can
//! produce a "reversal" (parent height below a child's), yielding a negative
//! branch length -- this matches scipy's average linkage and is left as-is.
//!
//! Throws std::invalid_argument on n < 2 or dist.size() != n*n.
std::vector<UpgmaNode> UpgmaAverageLinkage(const std::vector<double> &dist, uint32_t n);

} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the cluster_upgma table function into the extension catalog.
void RegisterClusterUpgma(ExtensionLoader &loader);
} // namespace duckdb
