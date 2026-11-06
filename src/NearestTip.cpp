#include "NearestTip.hpp"

#include <algorithm>
#include <functional>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace miint {

NearestTipTargets::NearestTipTargets(const TreeIndex &index, const std::vector<uint32_t> &nodes,
                                     const std::vector<uint32_t> &name_ranks)
    : index_(&index) {
	if (nodes.size() != name_ranks.size()) {
		throw std::invalid_argument("NearestTipTargets: the node and name rank arrays differ in size");
	}
	const size_t k = nodes.size();
	for (uint32_t node : nodes) {
		if (node >= index.size()) {
			throw std::invalid_argument("NearestTipTargets: target node " + std::to_string(node) + " is out of range");
		}
	}

	std::vector<size_t> order(k);
	std::iota(order.begin(), order.end(), size_t(0));
	std::sort(order.begin(), order.end(),
	          [&](size_t a, size_t b) { return index.pre(nodes[a]) < index.pre(nodes[b]); });

	pos_.resize(k);
	node_.resize(k);
	rank_.resize(k);
	seg_.assign(2 * k, std::numeric_limits<double>::infinity());
	for (size_t j = 0; j < k; ++j) {
		const uint32_t node = nodes[order[j]];
		pos_[j] = index.pre(node);
		if (j > 0 && pos_[j] == pos_[j - 1]) {
			throw std::invalid_argument("NearestTipTargets: target node " + std::to_string(node) + " appears twice");
		}
		node_[j] = node;
		rank_[j] = name_ranks[order[j]];
		seg_[k + j] = index.root_distance(node);
	}
	for (size_t i = k; i-- > 1;) {
		seg_[i] = std::min(seg_[2 * i], seg_[2 * i + 1]);
	}
}

NearestTipResult NearestTipTargets::Search(uint32_t query, const NearestTipOptions &options) const {
	const TreeIndex &tree = *index_;
	if (query >= tree.size()) {
		throw std::out_of_range("NearestTipTargets::Search: query node " + std::to_string(query) + " is out of range");
	}
	const size_t k = pos_.size();
	const double dq = tree.root_distance(query);

	NearestTipResult best;
	uint32_t best_rank = std::numeric_limits<uint32_t>::max();
	auto offer = [&](size_t slot, double d, uint32_t lca) {
		if (best.target == NearestTipResult::NO_TARGET || d < best.distance) {
			best = NearestTipResult {node_[slot], d, 1, lca};
			best_rank = rank_[slot];
		} else if (d == best.distance) {
			++best.n_tied;
			if (rank_[slot] < best_rank) {
				best.target = node_[slot];
				best.lca = lca;
				best_rank = rank_[slot];
			}
		}
	};
	// True when nothing at distance d or more can still win or tie.
	auto beaten = [&](double d) {
		return best.target != NearestTipResult::NO_TARGET && d > best.distance;
	};

	// Score the targets in slots [s, e), all of which meet the query at `lca`, in
	// increasing root distance. With up = D[q] - D[lca] fixed, a target's distance
	// up + (D[t] - D[lca]) never decreases as D[t] grows, so the enumeration stops at the
	// first one beyond the best; an internal segment-tree node's minimum bounds every
	// target under it the same way.
	using HeapItem = std::pair<double, size_t>;
	std::vector<HeapItem> heap;
	auto score_slice = [&](size_t s, size_t e, double up, double dl, uint32_t lca) {
		heap.clear();
		for (size_t l = s + k, r = e + k; l < r; l >>= 1, r >>= 1) {
			if (l & 1) {
				heap.emplace_back(seg_[l], l);
				++l;
			}
			if (r & 1) {
				--r;
				heap.emplace_back(seg_[r], r);
			}
		}
		std::make_heap(heap.begin(), heap.end(), std::greater<HeapItem>());
		while (!heap.empty()) {
			std::pop_heap(heap.begin(), heap.end(), std::greater<HeapItem>());
			const auto [d_min, seg_node] = heap.back();
			heap.pop_back();
			const double d = up + (d_min - dl);
			if (beaten(d)) {
				return;
			}
			if (seg_node >= k) {
				offer(seg_node - k, d, lca);
			} else {
				for (size_t child : {2 * seg_node, 2 * seg_node + 1}) {
					heap.emplace_back(seg_[child], child);
					std::push_heap(heap.begin(), heap.end(), std::greater<HeapItem>());
				}
			}
		}
	};

	// The covered part of the tree is a preorder interval, at first just the query. Slots
	// [0, below) lie before it and [above, k) after it. A query that is a target sits in
	// [below, above): scored only when it may be its own nearest target.
	const uint32_t pq = tree.pre(query);
	size_t below = std::lower_bound(pos_.begin(), pos_.end(), pq) - pos_.begin();
	size_t above = std::upper_bound(pos_.begin(), pos_.end(), pq) - pos_.begin();
	// For the clade scope, the lowest ancestor-or-self of the query holding an eligible
	// target, once known: the query itself when it is one.
	uint32_t clade = TreeIndex::NO_PARENT;
	if (above > below && !options.exclude_self) {
		offer(below, 0.0, query);
		clade = query;
	}

	uint32_t cur = query;
	while (below > 0 || above < k) {
		// The next level is the lowest ancestor that brings in an uncovered target. An
		// ancestor's subtree contains the covered interval, so it holds an uncovered
		// target before the interval exactly when it contains the nearest such target,
		// and likewise after it. Both candidates lie on the query's root path; the
		// deeper one is the next level.
		uint32_t level = TreeIndex::NO_PARENT;
		if (below > 0) {
			level = tree.lowest_ancestor_containing(cur, pos_[below - 1]);
		}
		if (above < k) {
			const uint32_t other = tree.lowest_ancestor_containing(cur, pos_[above]);
			if (level == TreeIndex::NO_PARENT || tree.depth(other) > tree.depth(level)) {
				level = other;
			}
		}

		// A clade search competes only within the first level holding an eligible target
		// (the query's own descendants count when the query is a target itself).
		if (options.scope == NearestTipScope::Clade && clade != TreeIndex::NO_PARENT && level != clade) {
			break;
		}

		// Every target from here up is at least D[q] - D[level] away, and that grows as
		// the levels rise (lengths are non-negative), so once it is beyond the best
		// nothing further can win or tie.
		const double dl = tree.root_distance(level);
		const double up = dq - dl;
		if (beaten(up)) {
			break;
		}

		const size_t new_below = std::lower_bound(pos_.begin(), pos_.begin() + below, tree.pre(level)) - pos_.begin();
		const size_t new_above = std::upper_bound(pos_.begin() + above, pos_.end(), tree.last(level)) - pos_.begin();
		score_slice(new_below, below, up, dl, level);
		score_slice(above, new_above, up, dl, level);
		// Every level brings in at least one target, so the first scored level is the
		// lowest one holding an eligible target.
		if (clade == TreeIndex::NO_PARENT) {
			clade = level;
		}

		below = new_below;
		above = new_above;
		cur = level;
	}
	return best;
}

} // namespace miint
