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
// A missing rank: NULL in the input row, or emptied by cleaning.
inline constexpr NameId kMissingName = std::numeric_limits<NameId>::max();

// One input row: a feature id, as text, and its values for the shared lineage schema's eight
// rank columns, in kLineageRanks order, nullopt for NULL.
struct TaxonomyRow {
	std::string feature_id;
	LineageRanks ranks;
};

// tax2tree's bad-name tokens, then "incertae sedis" and "incertae_sedis", which tax2tree
// lacks. A rank value containing one, case-insensitively, is a placeholder ("uncultured
// bacterium", "Clostridiales Family XVII. Incertae Sedis"), not a name.
const std::vector<std::string> &DefaultBadNameTokens();

// The kinds of problem a single row can have. taxonomy_violations returns them as rows;
// TaxonomyStore throws on them.
enum class TaxonomyViolationKind {
	// The feature_id is empty.
	InvalidId,
	// A rank value is not a bare name; RowViolation::reason says why.
	InvalidValue,
	// The feature_id is on more than one row; every one of those rows is reported.
	DuplicateId,
	// A rank is NULL while a deeper rank is not.
	Gap,
};

// Each kind's name in taxonomy_violations' `kind` column, the SQL output's contract.
const char *TaxonomyViolationKindName(TaxonomyViolationKind kind);
// The `kind` of a name under more than one parent (MultipleParents).
inline constexpr const char *kMultipleParentsKindName = "multiple_parents";

// One problem with one input row.
struct RowViolation {
	TaxonomyViolationKind kind;
	// The input row.
	size_t row;
	// InvalidValue: the value's rank. Gap: the first NULL rank. 0 for the id kinds.
	size_t rank;
	// InvalidValue: why the value is not a name ("is empty", ...). nullptr otherwise.
	const char *reason;
	// DuplicateId: how many rows carry the id. 0 otherwise.
	size_t count;
	// Gap: the first rank below `rank` that is set. 0 otherwise.
	size_t deeper;
};

// Every row problem, each row read as given, in input order; within a row, its id first, then
// its values by rank, then its gap.
std::vector<RowViolation> FindRowViolations(const std::vector<TaxonomyRow> &rows);

// A label, a rank and a whole name, under more than one parent. A label's parent is the value
// at the rank directly above it, and a NULL there is a parent too.
struct MultipleParents {
	// The label's rank: a kLineageRanks index of at least 1, since a domain has no parent.
	size_t rank;
	// One entry per distinct parent, in the order each was first seen: every row that puts the
	// label under that parent, in input order. The text is read from the rows: the label is
	// rows[rows_by_parent[p][0]].ranks[rank], and the parent is .ranks[rank - 1].
	std::vector<std::vector<size_t>> rows_by_parent;
};

// The labels under more than one parent, in the order in which each one's second parent first
// appears. `depths`, if given, holds one value per row: that row's ranks from depths[row] down
// are read as NULL (TaxonomyStore passes the depths of its cleaned lineages). Without it, every
// row is read as given, gaps included.
std::vector<MultipleParents> FindMultipleParents(const std::vector<TaxonomyRow> &rows,
                                                 const std::vector<size_t> *depths = nullptr);

// A lineage cut short by a bad-name token: `rank` and every deeper rank were emptied.
struct LineageTruncation {
	size_t row;
	size_t rank;
	std::string token;
};

// A taxonomy's lineages, validated, cleaned and interned per rank, with each name's parent.
// Every name has one parent, so walking Parent() from a name gives its one lineage.
//
// Construction throws std::invalid_argument from two passes, each giving the number of
// problems and the first 10 in input order:
//   1. Rows, all checked before anything is cleaned: a feature_id is empty or repeated; a
//      rank value is empty, starts or ends with whitespace, starts with a rank prefix
//      ("p__Bacillota") or contains ';'; a rank is NULL while a deeper rank is not (a gap).
//      Any problem stops construction here.
//   2. Names, over the cleaned lineages: a name is under more than one name at the rank
//      above it (its parent).
// It also throws, without a count, for an empty bad-name token, and for more distinct names
// at one rank than NameId can number.
//
// Cleaning, per row: the first rank whose value contains a bad-name token, and every rank
// below it, become missing; then a lineage missing its domain or its phylum becomes missing
// entirely (the row is non-informative). The bad-name cut runs first so that a lineage cut back to
// its domain is emptied too. tax2tree differs: it empties only the bad rank, which leaves a
// gap; its min-informative rule tests the phylum alone, before the bad-name filter; and its
// tokens do not include incertae sedis.
//
// Accessors take rows, ranks (kLineageRanks indexes) and NameIds in range; kMissingName is
// not a NameId.
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

	// Lineages cut at a bad name, in input order. A cut that left no domain or no phylum
	// also puts the row in WipedRows: the cut says why the lineage was emptied.
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
