#pragma once

#include "taxonomy_lineage.hpp"

#include <array>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <optional>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

namespace miint {

// A rank name, interned per rank: the same NameId at two ranks names two different taxa.
using NameId = uint32_t;
// A missing rank (SQL NULL in the input, or emptied by cleaning).
inline constexpr NameId kMissingName = std::numeric_limits<NameId>::max();

// One input row of the shared lineage schema: a feature id, as text, and its eight rank
// values in kLineageRanks order, nullopt for NULL.
struct TaxonomyRow {
	std::string feature_id;
	LineageRanks ranks;
};

// tax2tree's bad-name tokens. A rank value containing one, case-insensitively, is a
// placeholder ("uncultured bacterium"), not a name.
const std::vector<std::string> &DefaultBadNameTokens();

// A lineage cut short by a bad-name token: `rank` and every deeper rank were emptied.
struct LineageTruncation {
	size_t row;
	size_t rank;
	std::string token;
};

// A taxonomy in the shared lineage schema, validated, cleaned and interned for decoration.
//
// Construction checks every row and throws std::invalid_argument giving the number of
// problems and the first 10, in input order, when any of these holds:
//   - a feature_id is empty, or appears more than once;
//   - a rank value is empty, starts or ends with whitespace, starts with a rank prefix
//     ("p__Bacillota"), or contains ';';
//   - a rank is NULL while a deeper rank is not (a gap);
//   - after cleaning, a name has more than one parent name (one rank shallower).
// Cleaning, per row: the first rank whose value contains a bad-name token, and every
// rank below it, become missing; then a lineage missing its domain or phylum becomes
// missing entirely (the row is non-informative).
class TaxonomyStore {
public:
	explicit TaxonomyStore(const std::vector<TaxonomyRow> &rows,
	                       const std::vector<std::string> &bad_name_tokens = DefaultBadNameTokens());

	// Not copyable: row_of_id_ views the strings in feature_ids_, and a copy's views would
	// still point into the original. A move keeps the vector's buffer, so its views stay valid.
	TaxonomyStore(const TaxonomyStore &) = delete;
	TaxonomyStore &operator=(const TaxonomyStore &) = delete;
	TaxonomyStore(TaxonomyStore &&) = default;
	TaxonomyStore &operator=(TaxonomyStore &&) = default;

	size_t RowCount() const;
	const std::string &FeatureId(size_t row) const;
	// The input row with this feature_id, if any.
	std::optional<size_t> FindRow(std::string_view feature_id) const;

	// The row's cleaned name at `rank` (a kLineageRanks index), or kMissingName.
	NameId Name(size_t row, size_t rank) const;
	// Whether the row kept a lineage after cleaning.
	bool IsInformative(size_t row) const;

	// The number of distinct names at `rank`; NameIds at that rank are 0..NameCount-1.
	size_t NameCount(size_t rank) const;
	const std::string &NameText(size_t rank, NameId name) const;
	// The name's parent, one rank shallower; kMissingName at rank 0.
	NameId Parent(size_t rank, NameId name) const;

	// Lineages cut at a bad name, in input order. A cut that left no domain and phylum also
	// puts the row in WipedRows: the cut says why the lineage was emptied.
	const std::vector<LineageTruncation> &Truncations() const;
	// Rows whose lineage the domain/phylum rule emptied, in input order. A row that was
	// all NULL from the start is non-informative but not listed: nothing was removed.
	const std::vector<size_t> &WipedRows() const;

private:
	std::vector<std::string> feature_ids_;
	std::unordered_map<std::string_view, size_t> row_of_id_;
	// kLineageRanks.size() NameIds per row.
	std::vector<NameId> names_;
	// Per rank: the text of each NameId, and its parent's NameId.
	std::array<std::vector<std::string>, kLineageRanks.size()> name_text_;
	std::array<std::vector<NameId>, kLineageRanks.size()> parent_;
	std::vector<LineageTruncation> truncations_;
	std::vector<size_t> wiped_rows_;
};

} // namespace miint
