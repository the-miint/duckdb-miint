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

// See F1Score::Threshold for why min_count 0 is rejected.
void RequireMinCount(const char *caller, uint32_t min_count) {
	if (min_count == 0) {
		throw std::invalid_argument(std::string(caller) + ": min_count must be at least 1");
	}
}

// ceil(total / 2): the fewest of a name's tips that give recall of at least one half.
uint32_t CeilHalf(uint32_t total) {
	return total / 2 + total % 2;
}
} // namespace

uint32_t F1Score::Threshold(uint32_t total, uint32_t min_count) {
	RequireMinCount("F1Score", min_count);
	return std::max(min_count, CeilHalf(total));
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
	const uint32_t total = index.Total(rank, name);
	if (threshold == 0 || threshold < CeilHalf(total)) {
		throw std::invalid_argument("LowestEligibleNodes: threshold " + std::to_string(threshold) +
		                            " is below 1 or below half of the name's " + std::to_string(total) + " tips");
	}
	// The root holds every tip carrying the name.
	if (total < threshold) {
		return {};
	}
	// An eligible node holds at least threshold >= ceil(t/2) of the name's t tips, and they are a
	// contiguous run of its occurrences, so the run contains the middle one (odd t) or one of the
	// middle two (even t). Every eligible node is therefore an ancestor-or-self of a middle tip,
	// and since eligibility is closed upward, the lowest are among the deepest eligible ancestors
	// of the middle tips. Another tip's deepest eligible ancestor need not be lowest: it can hold
	// an eligible node that does not contain that tip.
	const TreeIndex &tree = index.Tree();
	// The root is eligible, so every walk ends.
	auto deepest_eligible = [&](uint32_t position) {
		uint32_t v = index.NodeAtPreorder(position);
		while (index.CladeCount(v, rank, name) < threshold) {
			v = tree.parent(v);
		}
		return v;
	};
	const auto occ = index.Occurrences(rank, name);
	const size_t t = occ.size();
	const uint32_t a = deepest_eligible(occ[(t - 1) / 2]);
	if (t % 2 == 1) {
		return {a};
	}
	// a's tip precedes b's in preorder. Nested results keep the deeper; disjoint ones are the two
	// halves of a 50/50 split, and a's subtree, holding the earlier tip, comes first.
	const uint32_t b = deepest_eligible(occ[t / 2]);
	if (tree.is_ancestor_or_self(a, b)) {
		return {b};
	}
	if (tree.is_ancestor_or_self(b, a)) {
		return {a};
	}
	return {a, b};
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

	// v's position in postorder. The nodes before it are its descendants, last(v) - pre(v) of
	// them, and the nodes before it in preorder that are not its ancestors, pre(v) - depth(v).
	auto postorder = [&](uint32_t v) {
		return tree.last(v) - tree.depth(v);
	};
	// Whether candidate a beats candidate b for one name: higher F1, then fewer tips, then later
	// in postorder.
	auto beats = [&](const NameScore &a, const NameScore &b) {
		if (const auto f1 = F1Score::Compare(a.c, a.m, b.c, b.m, a.total); f1 != 0) {
			return f1 > 0;
		}
		const uint32_t tips_a = index.TipCount(a.node), tips_b = index.TipCount(b.node);
		if (tips_a != tips_b) {
			return tips_a < tips_b;
		}
		return postorder(a.node) > postorder(b.node);
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
			// the ones that also pick rank r. At an internal one with rank r unsafe, another name is
			// eligible too, and so at every ancestor, so the walk stops. The second walk stops below
			// the LCA too: the LCA and its ancestors are on the first walk's path, which either went
			// through them or stopped below them at a node where every ancestor is unsafe.
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
					} else if (!tree.is_tip(v)) {
						return;
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

Backfill BackfillNames(const NameAssignment &assignment, const TaxonomyStore &store, const TreeIndex &tree) {
	const uint32_t n = static_cast<uint32_t>(tree.size());
	constexpr uint32_t kNone = TreeIndex::NO_PARENT;

	// Per node, itself if it has a preliminary name, else its nearest ancestor that has one, else
	// kNone; filled in preorder, so a node's parent is filled before it.
	std::vector<uint32_t> preorder(n);
	for (uint32_t v = 0; v < n; v++) {
		preorder[tree.pre(v)] = v;
	}
	std::vector<uint32_t> named(n, kNone);
	for (const uint32_t v : preorder) {
		const uint32_t parent = tree.parent(v);
		named[v] = assignment.Preliminary(v) ? v : parent == kNone ? kNone : named[parent];
	}

	Backfill out;
	out.labels.resize(n);
	for (uint32_t v = 0; v < n; v++) {
		const NameScore *own = assignment.Preliminary(v);
		if (own == nullptr) {
			continue;
		}
		// The shallowest rank the label reaches: 0 with no named ancestor, else the rank below the
		// ancestor's, and the node's own when that is deeper than the node's.
		size_t from = 0;
		const uint32_t parent = tree.parent(v);
		if (const uint32_t ancestor = parent == kNone ? kNone : named[parent]; ancestor != kNone) {
			const size_t ancestor_rank = assignment.Preliminary(ancestor)->rank;
			if (own->rank < ancestor_rank) {
				out.inverted.push_back({v, own->rank, own->name, ancestor, ancestor_rank});
			}
			from = std::min(ancestor_rank + 1, own->rank);
		}
		// The names at ranks from..own->rank up the hierarchy from the preliminary name. Every name
		// below rank 0 has a parent one rank shallower (TaxonomyStore::Parent), so the walk never
		// meets a gap.
		auto &label = out.labels[v];
		label.resize(own->rank - from + 1);
		NameId name = own->name;
		for (size_t r = own->rank;; r--) {
			label[r - from] = {r, name};
			if (r == from) {
				break;
			}
			name = store.Parent(r, name);
		}
	}
	return out;
}

std::string_view PolyphyleticBase(std::string_view name) {
	const size_t underscore = name.rfind('_');
	if (underscore == std::string_view::npos || underscore == 0 || underscore + 1 == name.size()) {
		return name;
	}
	const std::string_view suffix = name.substr(underscore + 1);
	const bool capitals = std::all_of(suffix.begin(), suffix.end(), [](char c) { return c >= 'A' && c <= 'Z'; });
	return capitals ? name.substr(0, underscore) : name;
}

namespace {
// The name the first names of a node's labeled descendants agree on: their one distinct name, or
// of two distinct names at one rank, the one whose PolyphyleticBase is the other.
std::optional<RankName> AgreedName(const std::vector<uint32_t> &descendants, const std::vector<NodeLabel> &labels,
                                   const TaxonomyStore &store) {
	std::optional<RankName> first, second;
	for (const uint32_t d : descendants) {
		const RankName &name = labels[d].front();
		if (!first || name == *first) {
			first = name;
		} else if (!second || name == *second) {
			second = name;
		} else {
			return std::nullopt;
		}
	}
	if (!second) {
		return first;
	}
	if (first->rank != second->rank) {
		return std::nullopt;
	}
	const std::string &a = store.NameText(first->rank, first->name);
	const std::string &b = store.NameText(second->rank, second->name);
	if (PolyphyleticBase(b) == a) {
		return second;
	}
	if (PolyphyleticBase(a) == b) {
		return first;
	}
	return std::nullopt;
}
} // namespace

std::vector<NodeLabel> PromoteCommonNames(std::vector<NodeLabel> labels, const TaxonomyStore &store,
                                          const TreeIndex &tree) {
	const uint32_t n = static_cast<uint32_t>(tree.size());
	// Each node's internal children: those of v are children[offsets[v] .. offsets[v + 1]).
	std::vector<uint32_t> offsets(n + 1, 0);
	for (uint32_t v = 0; v < n; v++) {
		if (v != tree.root() && !tree.is_tip(v)) {
			offsets[tree.parent(v) + 1]++;
		}
	}
	std::partial_sum(offsets.begin(), offsets.end(), offsets.begin());
	std::vector<uint32_t> children(offsets[n]);
	std::vector<uint32_t> next(offsets.begin(), offsets.end() - 1);
	std::vector<uint32_t> preorder(n);
	for (uint32_t v = 0; v < n; v++) {
		preorder[tree.pre(v)] = v;
		if (v != tree.root() && !tree.is_tip(v)) {
			children[next[tree.parent(v)]++] = v;
		}
	}

	std::vector<uint32_t> stack, descendants;
	for (const uint32_t v : preorder) {
		// v's nearest internal descendants with a non-empty label, found below the empty ones.
		descendants.clear();
		stack.assign(children.begin() + offsets[v], children.begin() + offsets[v + 1]);
		while (!stack.empty()) {
			const uint32_t d = stack.back();
			stack.pop_back();
			if (labels[d].empty()) {
				stack.insert(stack.end(), children.begin() + offsets[d], children.begin() + offsets[d + 1]);
			} else {
				descendants.push_back(d);
			}
		}
		if (descendants.empty()) {
			continue;
		}
		for (bool more = true; more;) {
			const auto name = AgreedName(descendants, labels, store);
			if (!name || (!labels[v].empty() && name->rank <= labels[v].back().rank)) {
				break;
			}
			labels[v].push_back(*name);
			for (const uint32_t d : descendants) {
				labels[d].erase(labels[d].begin());
				more = more && labels[d].size() > 1;
			}
		}
	}
	return labels;
}

std::string FormatLabel(const NodeLabel &label, const TaxonomyStore &store) {
	return {};
}

DecorateResult Decorate(const CladeCountIndex &index, const TaxonomyStore &store, uint32_t min_count) {
	return {AssignNames(index, min_count), {}, {}};
}

std::vector<LabelPart> ParseNodeLabel(std::string_view label) {
	return {};
}

EmittedLineages EmitLineages(const TreeIndex &tree, const std::vector<std::vector<LabelPart>> &labels) {
	return {};
}

} // namespace miint
