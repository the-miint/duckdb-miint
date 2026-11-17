#pragma once

#include "TaxonomyDecorator.hpp"
#include "TaxonomyStore.hpp"
#include "TreeIndex.hpp"
#include "tree_test_util.hpp"

#include <array>
#include <cstdint>
#include <optional>
#include <random>
#include <string>
#include <vector>

// Random trees with random taxonomies on their tips, and the counts CladeCountIndex should
// give for them, found by walking parents.

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

} // namespace miint::test
