#pragma once

#include "TreeIndex.hpp"
#include "catalog_utils.hpp"

#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/validity_mask.hpp"
#include "duckdb/main/client_context.hpp"

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Shared by the table functions that take a tree, and the relations that go with it, by name.

namespace duckdb {

// The type of `column` among a relation's columns, matched case-insensitively as DuckDB
// matches column names.
bool FindColumnType(const TableOrViewColumns &columns, const std::string &column, LogicalType &type);

// The type of the tree's `name` column, the type tip names are matched in. Throws
// BinderException, naming `fn`, when there is no such column or it is not an id type
// (VARCHAR, BIGINT or UUID).
LogicalType TipNameType(const char *fn, const std::string &tree_table, const TableOrViewColumns &tree_columns);

// The tree is read in a projection of its own, so passing it again as another relation
// argument would read it twice. Throws BinderException "<fn>: '<relation>' is the tree
// table; <remedy>" when `relation` names `tree_table` (identifiers are case-insensitive).
void RejectTreeAs(const char *fn, const std::string &tree_table, const std::string &relation, const char *remedy);

// Read `columns` (each an expression and the type it was bound with) from a relation, once,
// at execution, streamed over a helper connection (docs/internals/reading-tables-views.md),
// and hand each non-empty chunk to `consume`. `fn` names the function in errors. The
// expressions are spliced into the SQL as they are, and an error names a column by its
// expression, so pass constants of your own (a column name, or a cast of one), never
// anything a user supplied; the relation's name is quoted.
//
// Bind fixed the output types. A relation can change type without a catalog change (a
// view over a file that is rewritten, under a prepared statement), so a column read with
// another type is refused rather than emitted under the old one. The check uses the types
// of the read itself (docs/internals/reading-tables-views.md says why not the catalog's).
void ReadRelation(ClientContext &context, const char *fn, const std::string &relation,
                  const std::vector<std::pair<std::string, LogicalType>> &columns,
                  const std::function<void(DataChunk &)> &consume);

// "'a', 'b', ... (n more)": the first ten names, quoted.
std::string QuoteSome(const std::vector<std::string> &names);

// The tip a name names, and how many tips carry it (0: none; more than 1: ambiguous).
struct TipMatch {
	uint32_t node = 0;
	uint32_t count = 0;
};
// Node-based: a reference to an entry stays valid while the map grows, which callers that
// keep pointers to entries rely on.
using TipMatches = std::unordered_map<std::string, TipMatch>;

struct IndexedTree {
	miint::TreeIndex index;
	// The caller's node_index of each node of the index: ascending and distinct, so
	// miint::DenseTree::Find maps an id back to its node.
	std::vector<int64_t> node_ids;
};

// Read the tree once (ReadTreeArrays), keeping only the names that are keys of `tips`, and
// index it under `policy`; then fill in each name's TipMatch. Internal nodes are never
// matched. `name_type` is the type Bind saw for the tree's name column, which the read
// checks. A tree that cannot be numbered or indexed, or whose name column has changed
// type, is an InvalidInputException, "<fn>: the tree '<tree_table>' cannot be <use>: <why>".
IndexedTree ReadIndexedTree(ClientContext &context, const char *fn, const std::string &tree_table, TipMatches &tips,
                            const LogicalType &name_type, miint::LengthPolicy policy, const char *use);

// Throws InvalidInputException when a name names more than one tip, ending the message
// "<who> must name exactly one tip".
void RequireUnambiguousTips(const char *fn, const std::string &tree_table, const TipMatches &tips, const char *who);

// How many threads to run when `items` are claimed `per_claim` at a time: one per claim,
// at most the scheduler's thread count, and at least one.
idx_t ThreadsForClaims(ClientContext &context, idx_t items, idx_t per_claim);

// Write a distance into row `row` of a DOUBLE output vector's data and validity: NULL when
// it is unknown.
inline void SetDistance(double *data, ValidityMask &validity, idx_t row, const std::optional<double> &distance) {
	if (distance) {
		data[row] = *distance;
		validity.SetValid(row);
	} else {
		validity.SetInvalid(row);
	}
}

} // namespace duckdb
