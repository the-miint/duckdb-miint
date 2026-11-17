#include "TaxonomyDecorator.hpp"

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>

namespace miint {

namespace {
constexpr size_t kRanks = kLineageRanks.size();
} // namespace

uint32_t F1Score::Threshold(uint32_t total, uint32_t min_count) {
	if (min_count == 0) {
		throw std::invalid_argument("F1Score: min_count must be at least 1");
	}
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
	const auto first = occurrences_[rank].begin() + offsets_[rank][name];
	const auto last = occurrences_[rank].begin() + offsets_[rank][name + 1];
	const auto lo = std::lower_bound(first, last, tree_.pre(node));
	const auto hi = std::upper_bound(lo, last, tree_.last(node));
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

} // namespace miint
