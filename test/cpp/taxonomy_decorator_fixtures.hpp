#pragma once

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "tree_test_util.hpp"

#include <algorithm>
#include <array>
#include <compare>
#include <cstdint>
#include <optional>
#include <ostream>
#include <random>
#include <string>
#include <vector>

// Random trees with random taxonomies on their tips, and the counts CladeCountIndex should
// give for them, found by walking parents.

namespace miint {

// So Catch prints assignments it compares.
inline std::ostream &operator<<(std::ostream &os, const NameScore &s) {
	return os << "{rank " << s.rank << ", name " << s.name << ", node " << s.node << ", c " << s.c << ", m " << s.m
	          << ", total " << s.total << "}";
}

inline std::ostream &operator<<(std::ostream &os, const NotPlacedName &s) {
	return os << "{rank " << s.rank << ", name " << s.name << ", total " << s.total << ", reason "
	          << static_cast<int>(s.reason) << "}";
}

} // namespace miint

namespace miint::test {

constexpr size_t kRankCount = kLineageRanks.size();

// A tree, the taxonomy rows, and the row each node carries (tips only; nullopt for none).
struct TaxonomyCase {
	TreeArrays tree;
	std::vector<TaxonomyRow> rows;
	std::vector<std::optional<size_t>> row_of_node;
};

// A random tree of n nodes (random_tree: chains and multifurcations, numbered in no order),
// with a taxonomy drawn from a random name hierarchy, so no lineage has a gap or a name under
// two parents, as TaxonomyStore requires. Lineages stop at random depths (NULL tails, up to
// strain); some are domain-only or all NULL, some carry a bad name, so cleaning truncates or
// wipes them; some tips carry no row; and a few rows sit on no tip, so the store can hold
// names the tree does not.
inline TaxonomyCase random_taxonomy_case(std::mt19937 &rng, uint32_t n) {
	TaxonomyCase c;
	c.tree = random_tree(rng, n, true);
	c.row_of_node.assign(n, std::nullopt);

	// parent_name[r][k]: the index, at rank r - 1, of name k at rank r's parent.
	std::array<std::vector<size_t>, kRankCount> parent_name;
	for (size_t r = 0; r < kRankCount; r++) {
		const size_t count = 1 + std::uniform_int_distribution<size_t>(0, r + 1)(rng);
		for (size_t k = 0; k < count; k++) {
			parent_name[r].push_back(
			    r == 0 ? 0 : std::uniform_int_distribution<size_t>(0, parent_name[r - 1].size() - 1)(rng));
		}
	}
	auto name_text = [](size_t r, size_t k) {
		return "N" + std::to_string(r) + "_" + std::to_string(k);
	};
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	auto random_row = [&](std::string id) {
		TaxonomyRow row {std::move(id), {}};
		const double u = unit(rng);
		const size_t depth = u < 0.05 ? 0 : u < 0.15 ? 1 : std::uniform_int_distribution<size_t>(2, kRankCount)(rng);
		if (depth > 0) {
			size_t k = std::uniform_int_distribution<size_t>(0, parent_name[depth - 1].size() - 1)(rng);
			for (size_t r = depth; r-- > 0;) {
				row.ranks[r] = name_text(r, k);
				k = parent_name[r][k];
			}
			if (unit(rng) < 0.1) {
				const size_t bad = std::uniform_int_distribution<size_t>(0, depth - 1)(rng);
				row.ranks[bad] = "uncultured " + *row.ranks[bad];
			}
		}
		return row;
	};

	const auto internal = has_child(c.tree);
	for (uint32_t v = 0; v < n; v++) {
		if (!internal[v] && unit(rng) >= 0.15) {
			c.row_of_node[v] = c.rows.size();
			c.rows.push_back(random_row("T" + std::to_string(v)));
		}
	}
	const size_t extra = std::uniform_int_distribution<size_t>(0, 2)(rng);
	for (size_t k = 0; k < extra; k++) {
		c.rows.push_back(random_row("X" + std::to_string(k)));
	}
	return c;
}

inline CladeCountIndex build_index(const TaxonomyCase &c, const TaxonomyStore &store) {
	return CladeCountIndex(TreeIndex::Build(c.tree.parent, c.tree.length), store, c.row_of_node);
}

// The counts, found by walking parent pointers up from every tip. Preorder positions are
// never used, so this is independent of how the index finds a subtree.
struct NaiveCounts {
	// clade[rank][name][node]
	std::array<std::vector<std::vector<uint32_t>>, kRankCount> clade;
	std::vector<uint32_t> informative_tip_count;
	std::vector<uint32_t> tip_count;

	uint32_t Total(size_t rank, NameId name, uint32_t root) const {
		return clade[rank][name][root];
	}
};

inline NaiveCounts naive_counts(const TaxonomyCase &c, const TaxonomyStore &store) {
	const uint32_t n = static_cast<uint32_t>(c.tree.parent.size());
	NaiveCounts out;
	for (size_t r = 0; r < kRankCount; r++) {
		out.clade[r].assign(store.NameCount(r), std::vector<uint32_t>(n, 0));
	}
	out.informative_tip_count.assign(n, 0);
	out.tip_count.assign(n, 0);
	const auto internal = has_child(c.tree);
	for (uint32_t tip = 0; tip < n; tip++) {
		if (internal[tip]) {
			continue;
		}
		const auto &row = c.row_of_node[tip];
		const bool informative = row && store.IsInformative(*row);
		for (uint32_t v = tip; v != kNoParent; v = c.tree.parent[v]) {
			out.tip_count[v]++;
			if (!informative) {
				continue;
			}
			out.informative_tip_count[v]++;
			for (size_t r = 0; r < kRankCount; r++) {
				const NameId name = store.Name(*row, r);
				if (name != kMissingName) {
					out.clade[r][name][v]++;
				}
			}
		}
	}
	return out;
}

// A name's lowest eligible nodes by brute force: the nodes whose naive count reaches
// `threshold` while no child's does, in node order.
inline std::vector<uint32_t> naive_lowest_eligible(const TaxonomyCase &c, const NaiveCounts &naive, size_t rank,
                                                   NameId name, uint32_t threshold) {
	const auto &count = naive.clade[rank][name];
	const uint32_t n = static_cast<uint32_t>(c.tree.parent.size());
	std::vector<bool> has_eligible_child(n, false);
	for (uint32_t v = 0; v < n; v++) {
		if (count[v] >= threshold && c.tree.parent[v] != kNoParent) {
			has_eligible_child[c.tree.parent[v]] = true;
		}
	}
	std::vector<uint32_t> out;
	for (uint32_t v = 0; v < n; v++) {
		if (count[v] >= threshold && !has_eligible_child[v]) {
			out.push_back(v);
		}
	}
	return out;
}

// What AssignNames should give, by brute force over naive counts: every rule applied at every
// node directly, with no lowest-node search, no counting up the tree and no path walk.
struct NaiveAssignment {
	// By (rank, name).
	std::vector<NameScore> scores;
	std::vector<NotPlacedName> not_placed;
	// How many names' nodes were decided by each tie-break, so a test can require that both occur.
	size_t won_on_tip_count = 0;
	size_t won_on_postorder = 0;
};

inline NaiveAssignment naive_assign(const TaxonomyCase &c, const TaxonomyStore &store, const NaiveCounts &naive,
                                    uint32_t min_count) {
	const uint32_t n = static_cast<uint32_t>(c.tree.parent.size());
	const uint32_t root = root_of(c.tree);
	const auto internal = has_child(c.tree);

	// Each node's position in postorder, children visited in ascending index.
	std::vector<std::vector<uint32_t>> children(n);
	for (uint32_t v = 0; v < n; v++) {
		if (c.tree.parent[v] != kNoParent) {
			children[c.tree.parent[v]].push_back(v);
		}
	}
	std::vector<uint32_t> post(n);
	uint32_t next = 0;
	auto visit = [&](auto &&self, uint32_t v) -> void {
		for (uint32_t child : children[v]) {
			self(self, child);
		}
		post[v] = next++;
	};
	visit(visit, root);

	// safe_name[r][v]: the one name eligible at internal node v and rank r, else kMissingName.
	std::array<std::vector<NameId>, kRankCount> safe_name;
	std::array<std::vector<bool>, kRankCount> ever_eligible;
	std::vector<uint32_t> picked(n, 0);
	for (size_t r = 0; r < kRankCount; r++) {
		const uint32_t rank_min_count = r == kRankCount - 1 ? 1 : min_count;
		std::vector<uint32_t> eligible_names(n, 0);
		std::vector<NameId> eligible_name(n, kMissingName);
		ever_eligible[r].assign(store.NameCount(r), false);
		for (NameId name = 0; name < store.NameCount(r); name++) {
			const uint32_t threshold = F1Score::Threshold(naive.Total(r, name, root), rank_min_count);
			for (uint32_t v = 0; v < n; v++) {
				if (naive.clade[r][name][v] >= threshold) {
					eligible_names[v]++;
					eligible_name[v] = name;
					ever_eligible[r][name] = true;
				}
			}
		}
		safe_name[r].assign(n, kMissingName);
		for (uint32_t v = 0; v < n; v++) {
			if (internal[v] && eligible_names[v] == 1) {
				safe_name[r][v] = eligible_name[v];
			}
		}
	}
	// tax2tree's pick_names: from domain down, take each safe rank, and stop at the first unsafe
	// rank after one was taken.
	for (uint32_t v = 0; v < n; v++) {
		bool taken = false;
		for (size_t r = 0; r < kRankCount; r++) {
			if (safe_name[r][v] != kMissingName) {
				picked[v] |= 1u << r;
				taken = true;
			} else if (taken) {
				break;
			}
		}
	}

	NaiveAssignment out;
	for (size_t r = 0; r < kRankCount; r++) {
		for (NameId name = 0; name < store.NameCount(r); name++) {
			const uint32_t total = naive.Total(r, name, root);
			if (total == 0) {
				continue;
			}
			std::vector<NameScore> candidates;
			bool ever_safe = false;
			for (uint32_t v = 0; v < n; v++) {
				if (safe_name[r][v] != name) {
					continue;
				}
				ever_safe = true;
				if (picked[v] >> r & 1) {
					candidates.push_back({r, name, v, naive.clade[r][name][v], naive.informative_tip_count[v], total});
				}
			}
			if (candidates.empty()) {
				const auto reason = !ever_eligible[r][name] ? NotPlacedReason::kNeverEligible
				                    : !ever_safe            ? NotPlacedReason::kNeverRankSafe
				                                            : NotPlacedReason::kBlockedByContiguity;
				out.not_placed.push_back({r, name, total, reason});
				continue;
			}
			// F1 = 2c/(total+m), compared by cross-multiplying; then fewer tips; then later in
			// postorder.
			auto f1_cmp = [&](const NameScore &a, const NameScore &b) {
				return uint64_t(a.c) * (uint64_t(total) + b.m) <=> uint64_t(b.c) * (uint64_t(total) + a.m);
			};
			auto better = [&](const NameScore &a, const NameScore &b) {
				if (const auto f = f1_cmp(a, b); f != 0) {
					return f > 0;
				}
				if (naive.tip_count[a.node] != naive.tip_count[b.node]) {
					return naive.tip_count[a.node] < naive.tip_count[b.node];
				}
				return post[a.node] > post[b.node];
			};
			std::sort(candidates.begin(), candidates.end(), better);
			if (candidates.size() > 1 && f1_cmp(candidates[0], candidates[1]) == 0) {
				if (naive.tip_count[candidates[0].node] != naive.tip_count[candidates[1].node]) {
					out.won_on_tip_count++;
				} else {
					out.won_on_postorder++;
				}
			}
			out.scores.push_back(candidates[0]);
		}
	}
	return out;
}

} // namespace miint::test
