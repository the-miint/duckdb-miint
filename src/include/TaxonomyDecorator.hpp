#pragma once

#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace miint {

// The scoring rule that decides where a taxonomy name is placed on a tree: tax2tree's F1
// measure. For a candidate node, c is the number of informative tips below it carrying the
// name, m the number of informative tips below it, and total the number of informative tips
// in the tree carrying the name. Precision is c/m and recall c/total, so
// F1 = 2pr/(p+r) = 2c/(total+m), and candidates are compared by exact integer
// cross-multiplication, never in floating point.
//
// Contract. CladeCountIndex, and the placement built on it, are valid only while both hold:
//   1. A node is a candidate for a name exactly when its count c reaches a threshold that
//      depends only on (total, min_count). Counts never decrease toward the root, so the
//      candidates for a name are closed upward, and each name has at most two lowest
//      candidates, which is what makes placement a search over two root paths.
//   2. The score depends on (c, m, total) only.
// A rule that breaks (1), for example one that makes an ancestor ineligible while a
// descendant is eligible, invalidates that search. Placement would then need per-node count
// maps merged up the tree (small-to-large), which is much slower; nothing here supports it.
//
// The F-beta family generalizes this as (1+b^2)c/(b^2 total+m). For beta = 2 or 1/2 its cross
// products reach 5L^2 for L tips, against F1's 2L^2, which overflows uint64 at this cap: it
// would need a cap of 2^30 tips, or wider products.
struct F1Score {
	// The most tips a tree may have. Every count is at most this, so the cross products in
	// Compare, at most kMaxTips * 2 * kMaxTips = 2^63, fit in uint64.
	static constexpr uint32_t kMaxTips = uint32_t(1) << 31;

	// The count a node needs to be a candidate for a name: max(min_count, ceil(total/2)), so
	// recall is at least one half (2c >= total) and c is at least min_count. Throws
	// std::invalid_argument when min_count is 0, which would let a node with no informative
	// tips be a candidate, with precision 0/0.
	static uint32_t Threshold(uint32_t total, uint32_t min_count);

	// Orders candidate 1 against candidate 2 of one name by F1: greater when candidate 1 scores
	// higher. Requires 1 <= c <= m <= kMaxTips and c <= total <= kMaxTips for each.
	static std::strong_ordering Compare(uint32_t c1, uint32_t m1, uint32_t c2, uint32_t m2, uint32_t total);

	// F1 as a double, for reporting only; candidates are ordered by Compare.
	static double ScoreForReport(uint32_t c, uint32_t m, uint32_t total);
};

// Counts of taxonomy names over a tree's subtrees: for a node, a rank and a name, how many
// informative tips below the node carry the name, in two binary searches.
//
// Each subtree is a contiguous range of preorder positions (TreeIndex), so the tips carrying
// a name, listed by preorder position, hold the ones below a node as one contiguous run.
// Memory is O(n) per-node arrays plus one entry per (informative tip, named rank).
//
// A tip is informative when it carries a taxonomy row that kept a lineage after cleaning
// (TaxonomyStore::IsInformative). Other tips, unmapped or not, count in LeafCount only.
//
// Accessors take nodes, ranks (kLineageRanks indexes) and NameIds in range for the tree and
// store the index was built from.
class CladeCountIndex {
public:
	// `row_of_node` has one entry per node of `tree`: the store row a tip carries, or nullopt.
	// Throws std::invalid_argument when its size is not the tree's, when an internal node
	// carries a row, when a row is out of range for `store`, or when the tree has more than
	// F1Score::kMaxTips tips.
	CladeCountIndex(TreeIndex tree, const TaxonomyStore &store, const std::vector<std::optional<size_t>> &row_of_node);

	const TreeIndex &Tree() const;

	// Informative tips below `node` (itself included) carrying `name` at `rank`.
	uint32_t CladeCount(uint32_t node, size_t rank, NameId name) const;
	// Informative tips below `node`: the denominator of precision.
	uint32_t NumTips(uint32_t node) const;
	// All tips below `node`, informative or not.
	uint32_t LeafCount(uint32_t node) const;
	// Informative tips in the tree carrying `name` at `rank`: the denominator of recall.
	uint32_t Total(size_t rank, NameId name) const;

private:
	TreeIndex tree_;
	// Over preorder positions: entry p counts the tips (informative tips) at positions below p.
	std::vector<uint32_t> tip_prefix_;
	std::vector<uint32_t> informative_prefix_;
	// Per rank, the preorder positions of the informative tips carrying each name, ascending:
	// name k's are occurrences_[rank][offsets_[rank][k] .. offsets_[rank][k + 1]).
	std::array<std::vector<uint32_t>, kLineageRanks.size()> offsets_;
	std::array<std::vector<uint32_t>, kLineageRanks.size()> occurrences_;
};

} // namespace miint
