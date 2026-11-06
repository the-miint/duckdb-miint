#pragma once

#include "TreeIndex.hpp"
#include "catalog_utils.hpp"

#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/client_context.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// What the tree table functions that take relations by name share: finding a column's
// type at bind, reading a relation once at execution, and reading and indexing the tree
// with the names its tips are matched by.

namespace duckdb {

// The type of `column` among a relation's columns, matched case-insensitively as DuckDB
// matches column names.
bool FindColumnType(const TableOrViewColumns &columns, const std::string &column, LogicalType &type);

// Read `columns` (each an expression and the type it was bound with) from a relation, once,
// at execution, over a helper connection (docs/internals/reading-tables-views.md), and hand
// each non-empty chunk to `consume`. `fn` names the function in errors.
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
using TipMatches = std::unordered_map<std::string, TipMatch>;

struct IndexedTree {
	miint::TreeIndex index;
	// The caller's node_index of each node of the index.
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

} // namespace duckdb
