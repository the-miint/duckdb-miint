#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint::unifrac {

struct MetadataRow {
	std::string sample_id;
	std::string variable;
	std::string value;
};

struct NamedGrouping {
	std::string variable;
	std::vector<uint32_t> labels; // size == ordered_sample_ids.size()
	uint32_t n_groups;
};

// Factorize long-form metadata into per-variable grouping arrays aligned to
// `ordered_sample_ids` (the canonical sample order of the distance matrix).
//
// Contract:
//   - For each variable, each sample in `ordered_sample_ids` must have a row
//     in `rows`. Missing → throws std::invalid_argument naming the sample and
//     variable.
//   - Labels are *value-driven, not order-driven*: identical partitions
//     across variables produce identical label arrays. Encoded by
//     first-appearance of value in `ordered_sample_ids` order — independent
//     of input row order, so seeded PERMANOVA reproducibility holds across
//     metadata shuffles.
//   - `variables` empty: returns one grouping per distinct variable found in
//     `rows`, sorted lexicographically by variable name.
//   - `variables` lists a name not present in any row → throws naming the
//     variable.
//   - Rows whose sample_id is not in `ordered_sample_ids` are silently
//     ignored (metadata tables routinely oversample relative to the analysis
//     set).
//   - Duplicate (variable, sample_id) pairs with conflicting values throw
//     naming both keys (botched join / cartesian product). Duplicates with
//     identical values are accepted as a no-op.
//   - `ordered_sample_ids` empty → throws std::invalid_argument.
std::vector<NamedGrouping> BuildGroupings(const std::vector<MetadataRow> &rows,
                                          const std::vector<std::string> &ordered_sample_ids,
                                          const std::vector<std::string> &variables);

} // namespace miint::unifrac
