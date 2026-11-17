#pragma once

#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"

#include <array>
#include <compare>
#include <cstddef>
#include <cstdint>
#include <optional>
#include <span>
#include <vector>

namespace miint {

// The scoring rule that decides which node of a tree a taxonomy name is assigned to:
// tax2tree's F1 measure. For a node, c is the number of informative tips below it carrying
// the name, m the number of informative tips below it, and total the number of informative
// tips in the tree carrying the name. Precision is c/m and recall c/total, so
// F1 = 2pr/(p+r) = 2c/(total+m), and nodes are compared by exact integer cross-multiplication.
// tax2tree compares F1 as doubles, so two nodes whose F1 is the same rational can come out
// unequal there: c=2, m=2 against c=3, m=5 at total=4 gives 0.6666666666666666 against
// 0.6666666666666665. Here they tie.
//
// Contract. Decoration finds each name's lowest eligible nodes by searching at most two root
// paths, and that search relies on both of these:
//   1. A node is eligible for a name exactly when its count c reaches a threshold that depends
//      only on (total, min_count) and is at least 1 and at least ceil(total/2). Counts never
//      decrease toward the root, so a name's eligible nodes are closed upward. A threshold of
//      at least 1 keeps every node without the name ineligible, and one of at least
//      ceil(total/2) means two disjoint eligible nodes hold at least total between them, so
//      there are at most two lowest ones, and two only when each holds exactly half of an
//      even total.
//   2. The score depends on (c, m, total) only.
// A rule that breaks (1), for example one under which an ancestor is ineligible while a
// descendant is eligible, breaks that search. Decoration would then need per-node count maps
// merged up the tree, small-to-large: O(n log n) map operations to build, against O(n) work to
// build this index and two binary searches per count. Nothing here supports that.
//
// The F-beta family generalizes F1 as (1+b^2)c/(b^2 total+m). For beta = 2 or 1/2 its cross
// products reach 5L^2 for L tips, against F1's 2L^2, which overflows uint64 at kMaxTips: it
// would need a cap of 2^30 tips, or wider products.
struct F1Score {
	// The most tips a tree may have. Every count is at most this, so the cross products in
	// Compare, at most kMaxTips * 2 * kMaxTips = 2^63, fit in uint64.
	static constexpr uint32_t kMaxTips = uint32_t(1) << 31;

	// The count a node needs to be eligible for a name: max(min_count, ceil(total/2)), so
	// recall is at least one half (2c >= total) and c is at least min_count. Throws
	// std::invalid_argument when min_count is 0, which would let a node with no informative
	// tips be eligible, with precision 0/0.
	static uint32_t Threshold(uint32_t total, uint32_t min_count);

	// Orders two nodes' scores for one name by F1: greater when (c1, m1) scores higher than
	// (c2, m2). Requires 1 <= c <= m <= kMaxTips and c <= total <= kMaxTips for each.
	static std::strong_ordering Compare(uint32_t c1, uint32_t m1, uint32_t c2, uint32_t m2, uint32_t total);

	// F1 as a double, for reporting only; nodes are ordered by Compare.
	static double ScoreForReport(uint32_t c, uint32_t m, uint32_t total);
};

// Counts of taxonomy names over a tree's subtrees: for a node, a rank and a name, how many
// informative tips below the node carry the name, in two binary searches.
//
// Each subtree is a contiguous range of preorder positions (TreeIndex), so the tips carrying
// a name, listed by preorder position, hold the ones below a node as one contiguous run.
// Memory is O(n) per-node arrays, one entry per (informative tip, named rank), and one offset
// per name in the store.
//
// A tip is informative when it carries a taxonomy row that kept a lineage after cleaning
// (TaxonomyStore::IsInformative). A tip with no row, or with a row cleaning emptied, counts
// in TipCount only.
//
// Accessors take nodes, ranks (kLineageRanks indexes) and NameIds in range for the tree and
// store the index was built from; kMissingName is not a NameId.
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
	// Informative tips below `node`: the denominator of precision. tax2tree calls it NumTips.
	uint32_t InformativeTipCount(uint32_t node) const;
	// All tips below `node`, informative or not: what tax2tree breaks ties between equally
	// scored nodes with.
	uint32_t TipCount(uint32_t node) const;
	// Informative tips in the tree carrying `name` at `rank`: the denominator of recall.
	uint32_t Total(size_t rank, NameId name) const;

	// The number of names at `rank` in the store the index was built from; NameIds at that
	// rank are 0..NameCount-1.
	size_t NameCount(size_t rank) const;
	// The preorder positions of the informative tips carrying `name` at `rank`, ascending.
	std::span<const uint32_t> Occurrences(size_t rank, NameId name) const;
	// The node at preorder position `p`.
	uint32_t NodeAtPreorder(uint32_t p) const;

private:
	TreeIndex tree_;
	// Over preorder positions: entry p counts the tips (informative tips) at positions less than p.
	std::vector<uint32_t> tip_prefix_;
	std::vector<uint32_t> informative_prefix_;
	// Per rank, the preorder positions of the informative tips carrying each name, ascending:
	// name k's are occurrences_[rank][offsets_[rank][k] .. offsets_[rank][k + 1]).
	std::array<std::vector<uint32_t>, kLineageRanks.size()> offsets_;
	std::array<std::vector<uint32_t>, kLineageRanks.size()> occurrences_;
};

// A name's lowest eligible nodes at `rank`: the nodes whose count reaches `threshold` while no
// child's does, in ascending preorder. There are none when even the root falls short of
// `threshold`, and two only when each holds exactly half of an even total. Requires a threshold
// of at least 1 and at least ceil(Total/2), as F1Score::Threshold gives.
std::vector<uint32_t> LowestEligibleNodes(const CladeCountIndex &index, size_t rank, NameId name, uint32_t threshold);

} // namespace miint
