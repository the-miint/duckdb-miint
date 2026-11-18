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
#include <set>
#include <string>
#include <utility>
#include <vector>

// Random trees with random taxonomies on their tips; the counts CladeCountIndex should give
// for them, found by walking parents; and, from those counts, brute-force versions of each
// stage built on the index, which the tests check the stages against.

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
constexpr size_t kStrainRank = kRankCount - 1;
constexpr size_t kSpeciesRank = kRankCount - 2;

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

// The min_count a rank uses, written out from the rule rather than taken from MinCountAt: 1 at
// strain.
inline uint32_t naive_min_count(size_t rank, uint32_t min_count) {
	return rank == kStrainRank ? 1 : min_count;
}

// The count a node needs at `rank` to be eligible for a name `total` tips carry.
inline uint32_t naive_threshold(size_t rank, uint32_t total, uint32_t min_count) {
	return F1Score::Threshold(total, naive_min_count(rank, min_count));
}

// Per rank and node, how many names are eligible there, and the last one found, which is the
// only one where the count is 1.
struct NaiveEligible {
	std::array<std::vector<uint32_t>, kRankCount> count;
	std::array<std::vector<NameId>, kRankCount> name;
};

inline NaiveEligible naive_eligible(const TaxonomyCase &c, const TaxonomyStore &store, const NaiveCounts &naive,
                                    uint32_t min_count) {
	const uint32_t n = static_cast<uint32_t>(c.tree.parent.size());
	const uint32_t root = root_of(c.tree);
	NaiveEligible out;
	for (size_t r = 0; r < kRankCount; r++) {
		out.count[r].assign(n, 0);
		out.name[r].assign(n, kMissingName);
		for (NameId name = 0; name < store.NameCount(r); name++) {
			const uint32_t threshold = naive_threshold(r, naive.Total(r, name, root), min_count);
			for (uint32_t v = 0; v < n; v++) {
				if (naive.clade[r][name][v] >= threshold) {
					out.count[r][v]++;
					out.name[r][v] = name;
				}
			}
		}
	}
	return out;
}

// tax2tree's pick_names over a safe mask: from domain down, take each safe rank, and stop at the
// first unsafe rank after one was taken.
inline unsigned naive_picked(unsigned safe_mask) {
	unsigned picked = 0;
	bool taken = false;
	for (size_t r = 0; r < kRankCount; r++) {
		if (safe_mask >> r & 1) {
			picked |= 1u << r;
			taken = true;
		} else if (taken) {
			break;
		}
	}
	return picked;
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
	const auto eligible = naive_eligible(c, store, naive, min_count);
	std::array<std::vector<NameId>, kRankCount> safe_name;
	std::vector<unsigned> safe_mask(n, 0);
	for (size_t r = 0; r < kRankCount; r++) {
		safe_name[r].assign(n, kMissingName);
		for (uint32_t v = 0; v < n; v++) {
			if (internal[v] && eligible.count[r][v] == 1) {
				safe_name[r][v] = eligible.name[r][v];
				safe_mask[v] |= 1u << r;
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
				if (naive_picked(safe_mask[v]) >> r & 1) {
					candidates.push_back({r, name, v, naive.clade[r][name][v], naive.informative_tip_count[v], total});
				}
			}
			if (candidates.empty()) {
				const auto reason = total < naive_min_count(r, min_count) ? NotPlacedReason::kNeverEligible
				                    : !ever_safe                          ? NotPlacedReason::kNeverRankSafe
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

// Over many naive_assign results, whether each tie-break decided some name and each reason left
// some name unassigned. A test over random cases requires that all did, or its agreement shows
// nothing about the ones that didn't.
struct NaiveCoverage {
	size_t won_on_tip_count = 0;
	size_t won_on_postorder = 0;
	std::set<NotPlacedReason> reasons;

	void Add(const NaiveAssignment &assignment) {
		won_on_tip_count += assignment.won_on_tip_count;
		won_on_postorder += assignment.won_on_postorder;
		for (const auto &name : assignment.not_placed) {
			reasons.insert(name.reason);
		}
	}

	// What never occurred; empty when everything did.
	std::vector<std::string> Missing() const {
		std::vector<std::string> missing;
		if (won_on_tip_count == 0) {
			missing.push_back("a name decided by fewer tips");
		}
		if (won_on_postorder == 0) {
			missing.push_back("a name decided by postorder");
		}
		const std::pair<NotPlacedReason, const char *> all_reasons[] = {
		    {NotPlacedReason::kNeverEligible, "kNeverEligible"},
		    {NotPlacedReason::kNeverRankSafe, "kNeverRankSafe"},
		    {NotPlacedReason::kBlockedByContiguity, "kBlockedByContiguity"}};
		for (const auto &[reason, label] : all_reasons) {
			if (!reasons.contains(reason)) {
				missing.push_back(std::string("a name not placed for ") + label);
			}
		}
		return missing;
	}
};

} // namespace miint::test
