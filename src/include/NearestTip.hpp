#pragma once

#include "TreeIndex.hpp"

#include <cstdint>
#include <vector>

namespace miint {

enum class NearestTipScope {
	// The nearest eligible target by path length, anywhere in the tree.
	Path,
	// The lowest ancestor-or-self of the query whose subtree holds an eligible target,
	// then the nearest eligible target by path length within that subtree. This can be
	// farther than the path-nearest target, which may lie outside the smallest clade.
	Clade,
};

struct NearestTipOptions {
	NearestTipScope scope = NearestTipScope::Path;
	// A query that is itself a target ignores itself.
	bool exclude_self = true;
};

struct NearestTipResult {
	static constexpr uint32_t NO_TARGET = UINT32_MAX;

	uint32_t target = NO_TARGET; // tree node of the chosen target; NO_TARGET when none is eligible
	double distance = 0.0;       // query-to-target path length; meaningful only with a target
	uint64_t n_tied = 0;         // eligible targets at exactly `distance`; 0 when there is none
	uint32_t lca = NO_TARGET;    // lowest common ancestor of the query and the chosen target
};

// One set of eligible targets (one group), prepared once and searched by many queries.
//
// The distance between query q and target t, whose lowest common ancestor is l, is
// (D[q] - D[l]) + (D[t] - D[l]) with D the root distances of the TreeIndex. A pair's
// distance therefore depends only on the pair: not on the order of children and not on
// which other targets exist. The nearest target is the one with the smallest distance;
// a tie (exactly equal doubles) goes to the smallest name rank, and every target at that
// distance is counted.
//
// Search climbs from the query through only those ancestors under which a new target
// appears, scoring each level's targets in increasing root distance (at a fixed level,
// distance is monotone in it), and stops once no higher level can be closer, or, for the
// clade scope, after the first level that holds an eligible target.
//
// Groups are separate instances over one TreeIndex; Search is const and may be called
// from several threads at once.
class NearestTipTargets {
public:
	// `nodes[i]` is a target's tree node and `name_ranks[i]` its place in name order (a
	// tie goes to the smallest). Throws std::invalid_argument when the arrays differ in
	// size, a node is out of range, or a node appears twice. `index` must outlive this.
	NearestTipTargets(const TreeIndex &index, const std::vector<uint32_t> &nodes,
	                  const std::vector<uint32_t> &name_ranks);

	size_t size() const {
		return pos_.size();
	}

	NearestTipResult Search(uint32_t query, const NearestTipOptions &options = {}) const;

private:
	const TreeIndex *index_;
	// Targets in ascending preorder position: position, tree node and name rank.
	std::vector<uint32_t> pos_;
	std::vector<uint32_t> node_;
	std::vector<uint32_t> rank_;
	// A bottom-up segment tree of the targets' root distances in that order: leaves at
	// [size, 2 * size), each internal node the minimum of its two children.
	std::vector<double> seg_;
};

} // namespace miint
