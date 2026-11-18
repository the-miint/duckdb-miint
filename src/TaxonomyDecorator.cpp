#include "TaxonomyDecorator.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace miint {

namespace {
constexpr size_t kRanks = kLineageRanks.size();
static_assert(kRanks <= 8, "a rank mask is one byte");

// min_count 0 would let a node with no informative tips be eligible, with precision 0/0.
void RequireMinCount(const char *caller, uint32_t min_count) {
	if (min_count == 0) {
		throw std::invalid_argument(std::string(caller) + ": min_count must be at least 1");
	}
}
} // namespace

uint32_t F1Score::Threshold(uint32_t total, uint32_t min_count) {
	RequireMinCount("F1Score", min_count);
	return std::max(min_count, total / 2 + total % 2);
}

std::strong_ordering F1Score::Compare(uint32_t c1, uint32_t m1, uint32_t c2, uint32_t m2, uint32_t total) {
	// 2c1/(total+m1) against 2c2/(total+m2), cross-multiplied; kMaxTips keeps the products in
	// uint64.
	const uint64_t lhs = uint64_t(c1) * (uint64_t(total) + m2);
	const uint64_t rhs = uint64_t(c2) * (uint64_t(total) + m1);
	return lhs <=> rhs;
}

double F1Score::ScoreForReport(uint32_t c, uint32_t m, uint32_t total) {
	return 2.0 * c / (double(total) + m);
}

CladeCountIndex::CladeCountIndex(TreeIndex tree, const TaxonomyStore &store,
                                 const std::vector<std::optional<size_t>> &row_of_node)
    : tree_(std::move(tree)) {
	const uint32_t n = static_cast<uint32_t>(tree_.size());
	if (row_of_node.size() != n) {
		throw std::invalid_argument("CladeCountIndex: row_of_node has " + std::to_string(row_of_node.size()) +
		                            " entries for a tree of " + std::to_string(n) + " nodes");
	}
	// Tips are visited in preorder below, through node_at_, so each name's occurrences come out
	// ascending.
	node_at_.resize(n);
	size_t tips = 0;
	for (uint32_t v = 0; v < n; v++) {
		node_at_[tree_.pre(v)] = v;
		tips += tree_.is_tip(v);
		const auto &row = row_of_node[v];
		if (!row) {
			continue;
		}
		if (!tree_.is_tip(v)) {
			throw std::invalid_argument("CladeCountIndex: node " + std::to_string(v) +
			                            " is not a tip but carries row " + std::to_string(*row));
		}
		if (*row >= store.RowCount()) {
			throw std::invalid_argument("CladeCountIndex: node " + std::to_string(v) + " carries row " +
			                            std::to_string(*row) + ", but the taxonomy has " +
			                            std::to_string(store.RowCount()) + " rows");
		}
	}
	if (tips > F1Score::kMaxTips) {
		throw std::invalid_argument("CladeCountIndex: the tree has " + std::to_string(tips) + " tips, more than the " +
		                            std::to_string(F1Score::kMaxTips) + " allowed");
	}

	auto informative = [&](uint32_t v) {
		return row_of_node[v] && store.IsInformative(*row_of_node[v]);
	};
	// Calls f(rank, name) for each named rank of v's row, if v is an informative tip.
	auto for_each_name = [&](uint32_t v, auto &&f) {
		if (!informative(v)) {
			return;
		}
		for (size_t r = 0; r < kRanks; r++) {
			const NameId name = store.Name(*row_of_node[v], r);
			if (name != kMissingName) {
				f(r, name);
			}
		}
	};

	// First pass, in preorder: the prefix sums, and how many tips carry each name.
	tip_prefix_.assign(n + 1, 0);
	informative_prefix_.assign(n + 1, 0);
	for (size_t r = 0; r < kRanks; r++) {
		offsets_[r].assign(store.NameCount(r) + 1, 0);
	}
	for (uint32_t p = 0; p < n; p++) {
		const uint32_t v = node_at_[p];
		tip_prefix_[p + 1] = tip_prefix_[p] + tree_.is_tip(v);
		informative_prefix_[p + 1] = informative_prefix_[p] + informative(v);
		for_each_name(v, [&](size_t r, NameId name) { offsets_[r][name + 1]++; });
	}

	// Second pass, in preorder: each tip's position into its names' runs.
	std::array<std::vector<uint32_t>, kRanks> next;
	for (size_t r = 0; r < kRanks; r++) {
		std::partial_sum(offsets_[r].begin(), offsets_[r].end(), offsets_[r].begin());
		occurrences_[r].resize(offsets_[r].back());
		next[r].assign(offsets_[r].begin(), offsets_[r].end() - 1);
	}
	for (uint32_t p = 0; p < n; p++) {
		for_each_name(node_at_[p], [&](size_t r, NameId name) { occurrences_[r][next[r][name]++] = p; });
	}
}

const TreeIndex &CladeCountIndex::Tree() const {
	return tree_;
}

uint32_t CladeCountIndex::CladeCount(uint32_t node, size_t rank, NameId name) const {
	const auto occ = Occurrences(rank, name);
	const auto lo = std::lower_bound(occ.begin(), occ.end(), tree_.pre(node));
	const auto hi = std::upper_bound(lo, occ.end(), tree_.last(node));
	return static_cast<uint32_t>(hi - lo);
}

uint32_t CladeCountIndex::InformativeTipCount(uint32_t node) const {
	return informative_prefix_[tree_.last(node) + 1] - informative_prefix_[tree_.pre(node)];
}

uint32_t CladeCountIndex::TipCount(uint32_t node) const {
	return tip_prefix_[tree_.last(node) + 1] - tip_prefix_[tree_.pre(node)];
}

uint32_t CladeCountIndex::Total(size_t rank, NameId name) const {
	return offsets_[rank][name + 1] - offsets_[rank][name];
}

size_t CladeCountIndex::NameCount(size_t rank) const {
	return offsets_[rank].size() - 1;
}

std::span<const uint32_t> CladeCountIndex::Occurrences(size_t rank, NameId name) const {
	return std::span<const uint32_t>(occurrences_[rank]).subspan(offsets_[rank][name], Total(rank, name));
}

uint32_t CladeCountIndex::NodeAtPreorder(uint32_t p) const {
	return node_at_[p];
}

std::vector<uint32_t> LowestEligibleNodes(const CladeCountIndex &index, size_t rank, NameId name, uint32_t threshold) {
	// An eligible node holds at least threshold >= ceil(t/2) of the name's t tips, and they are a
	// contiguous run of its occurrences, so the run contains the middle one (odd t) or one of the
	// middle two (even t). Every eligible node is therefore an ancestor-or-self of a middle tip,
	// and since eligibility is closed upward, the lowest are among the deepest eligible ancestors
	// of the middle tips. Another tip's deepest eligible ancestor need not be lowest: it can hold
	// an eligible node that does not contain that tip.
	const TreeIndex &tree = index.Tree();
	auto deepest_eligible = [&](uint32_t position) -> std::optional<uint32_t> {
		for (uint32_t v = index.NodeAtPreorder(position);; v = tree.parent(v)) {
			if (index.CladeCount(v, rank, name) >= threshold) {
				return v;
			}
			if (v == tree.root()) {
				return std::nullopt;
			}
		}
	};
	const auto occ = index.Occurrences(rank, name);
	const size_t t = occ.size();
	if (t == 0) {
		return {};
	}
	const auto a = deepest_eligible(occ[(t - 1) / 2]);
	if (!a) {
		return {};
	}
	if (t % 2 == 1) {
		return {*a};
	}
	// b exists, since a does only when the root is eligible. a's tip precedes b's in preorder.
	// Nested results keep the deeper; disjoint ones are the two halves of a 50/50 split, and a's
	// subtree, holding the earlier tip, comes first.
	const auto b = deepest_eligible(occ[t / 2]);
	if (tree.is_ancestor_or_self(*a, *b)) {
		return {*b};
	}
	if (tree.is_ancestor_or_self(*b, *a)) {
		return {*a};
	}
	return {*a, *b};
}

uint32_t MinCountAt(size_t rank, uint32_t min_count) {
	return rank == kRanks - 1 ? 1 : min_count;
}

namespace {

// Per rank, each name's lowest eligible nodes, by NameId.
using LowestByName = std::array<std::vector<std::vector<uint32_t>>, kRanks>;

LowestByName FindLowest(const CladeCountIndex &index, uint32_t min_count) {
	LowestByName lowest;
	for (size_t r = 0; r < kRanks; r++) {
		const uint32_t rank_min_count = MinCountAt(r, min_count);
		lowest[r].resize(index.NameCount(r));
		for (NameId name = 0; name < index.NameCount(r); name++) {
			const uint32_t threshold = F1Score::Threshold(index.Total(r, name), rank_min_count);
			lowest[r][name] = LowestEligibleNodes(index, r, name, threshold);
		}
	}
	return lowest;
}

std::vector<uint8_t> SafeMasks(const CladeCountIndex &index, const LowestByName &lowest) {
	const TreeIndex &tree = index.Tree();
	const uint32_t n = static_cast<uint32_t>(tree.size());
	std::vector<uint8_t> masks(n, 0);
	// Per preorder position, the number of names eligible at that node. Eligibility is closed
	// upward, so that is the number of names with a lowest eligible node at or below it: each
	// lowest node adds 1, and a name with two takes 1 back at their LCA, so it counts once from
	// there up. Before the accumulation a count can wrap below 0; unsigned arithmetic is modular,
	// and every accumulated count, at most the tree's tips, is exact. Counts are not capped at 2:
	// a sum capped before an LCA takes its 1 back would give 1 where two names are eligible, and
	// call the rank safe.
	std::vector<uint32_t> eligible(n);
	for (size_t r = 0; r < kRanks; r++) {
		std::fill(eligible.begin(), eligible.end(), 0);
		for (const auto &nodes : lowest[r]) {
			for (uint32_t v : nodes) {
				eligible[tree.pre(v)]++;
			}
			if (nodes.size() == 2) {
				eligible[tree.pre(tree.lca(nodes[0], nodes[1]))]--;
			}
		}
		// A node follows its parent in preorder, so in reverse preorder each count is complete
		// before it is added to the parent's. Position 0 is the root.
		for (uint32_t p = n; p-- > 1;) {
			eligible[tree.pre(tree.parent(index.NodeAtPreorder(p)))] += eligible[p];
		}
		for (uint32_t p = 0; p < n; p++) {
			const uint32_t v = index.NodeAtPreorder(p);
			if (eligible[p] == 1 && !tree.is_tip(v)) {
				masks[v] |= uint8_t(1) << r;
			}
		}
	}
	return masks;
}

} // namespace

std::vector<uint8_t> RankSafeMasks(const CladeCountIndex &index, uint32_t min_count) {
	RequireMinCount("RankSafeMasks", min_count);
	return SafeMasks(index, FindLowest(index, min_count));
}

uint8_t PickedRanks(uint8_t safe_mask) {
	// Adding the lowest set bit carries through the run of set bits that starts there, clearing
	// it, so the bits that m has and m + low lacks are that run.
	const uint32_t m = safe_mask;
	const uint32_t low = m & (~m + 1);
	return static_cast<uint8_t>(m & ~(m + low));
}

NameAssignment::NameAssignment(std::vector<NameScore> scores, std::vector<NotPlacedName> not_placed,
                               uint32_t node_count)
    : scores_(std::move(scores)), not_placed_(std::move(not_placed)), node_offsets_(node_count + 1, 0) {
	// scores_ is by (rank, name), so a stable counting sort by node keeps each node's names by
	// ascending rank.
	for (const auto &score : scores_) {
		node_offsets_[score.node + 1]++;
	}
	std::partial_sum(node_offsets_.begin(), node_offsets_.end(), node_offsets_.begin());
	by_node_.resize(scores_.size());
	std::vector<uint32_t> next(node_offsets_.begin(), node_offsets_.end() - 1);
	for (const auto &score : scores_) {
		by_node_[next[score.node]++] = score;
	}
}

std::span<const NameScore> NameAssignment::NamesAt(uint32_t node) const {
	return std::span<const NameScore>(by_node_).subspan(node_offsets_[node],
	                                                    node_offsets_[node + 1] - node_offsets_[node]);
}

const NameScore *NameAssignment::Preliminary(uint32_t node) const {
	const auto names = NamesAt(node);
	return names.empty() ? nullptr : &names.back();
}

const std::vector<NameScore> &NameAssignment::Scores() const {
	return scores_;
}

const std::vector<NotPlacedName> &NameAssignment::NotPlaced() const {
	return not_placed_;
}

NameAssignment AssignNames(const CladeCountIndex &index, uint32_t min_count) {
	RequireMinCount("AssignNames", min_count);
	const TreeIndex &tree = index.Tree();
	const auto lowest = FindLowest(index, min_count);
	const auto masks = SafeMasks(index, lowest);

	// Whether candidate a beats candidate b for one name: higher F1, then fewer tips, then later
	// in postorder. An ancestor comes after its descendants in postorder, and of two disjoint
	// nodes the one later in preorder comes later in postorder too.
	auto beats = [&](const NameScore &a, const NameScore &b) {
		if (const auto f1 = F1Score::Compare(a.c, a.m, b.c, b.m, a.total); f1 != 0) {
			return f1 > 0;
		}
		const uint32_t tips_a = index.TipCount(a.node), tips_b = index.TipCount(b.node);
		if (tips_a != tips_b) {
			return tips_a < tips_b;
		}
		if (tree.is_ancestor_or_self(a.node, b.node) || tree.is_ancestor_or_self(b.node, a.node)) {
			return tree.is_ancestor_or_self(a.node, b.node);
		}
		return tree.pre(a.node) > tree.pre(b.node);
	};

	std::vector<NameScore> scores;
	std::vector<NotPlacedName> not_placed;
	for (size_t r = 0; r < kRanks; r++) {
		const uint8_t bit = uint8_t(1) << r;
		for (NameId name = 0; name < index.NameCount(r); name++) {
			const uint32_t total = index.Total(r, name);
			if (total == 0) {
				continue;
			}
			const auto &nodes = lowest[r][name];
			if (nodes.empty()) {
				not_placed.push_back({r, name, total, NotPlacedReason::kNeverEligible});
				continue;
			}
			// The nodes where the name is eligible are the lowest ones and their ancestors. At such
			// a node with rank r safe, this name is the one eligible there, so its candidates are
			// the ones that also pick rank r. The second walk stops below the LCA, which the first
			// walk reached.
			std::optional<NameScore> best;
			bool ever_safe = false;
			auto walk = [&](uint32_t from, std::optional<uint32_t> stop) {
				for (uint32_t v = from;; v = tree.parent(v)) {
					if (stop && v == *stop) {
						return;
					}
					if (masks[v] & bit) {
						ever_safe = true;
						if (PickedRanks(masks[v]) & bit) {
							const NameScore candidate {
							    r, name, v, index.CladeCount(v, r, name), index.InformativeTipCount(v), total};
							if (!best || beats(candidate, *best)) {
								best = candidate;
							}
						}
					}
					if (v == tree.root()) {
						return;
					}
				}
			};
			walk(nodes[0], std::nullopt);
			if (nodes.size() == 2) {
				walk(nodes[1], tree.lca(nodes[0], nodes[1]));
			}
			if (best) {
				scores.push_back(*best);
			} else {
				not_placed.push_back(
				    {r, name, total,
				     ever_safe ? NotPlacedReason::kBlockedByContiguity : NotPlacedReason::kNeverRankSafe});
			}
		}
	}
	return NameAssignment(std::move(scores), std::move(not_placed), static_cast<uint32_t>(tree.size()));
}

} // namespace miint
