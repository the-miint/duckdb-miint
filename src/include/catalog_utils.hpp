#pragma once

#include "duckdb/common/types.hpp"
#include "duckdb/main/client_context.hpp"
#include <string>
#include <vector>

namespace duckdb {

struct TableOrViewColumns {
	vector<string> names;
	vector<LogicalType> types;
	bool is_physical_table; // true for tables (have rowid), false for views
};

// Look up a table or view by name and return its column names and types.
// entity_type is used in error messages (e.g., "Placement table", "Tree table").
// Throws BinderException if the table/view does not exist.
TableOrViewColumns GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                         const std::string &entity_type = "Table");

// Case-insensitive check for whether `columns` contains a column named `col`.
bool HasColumn(const TableOrViewColumns &columns, const std::string &col);

// Guard for functions whose argument is a LITERAL identifier fetched from a remote
// service (e.g. an NCBI accession), NOT the name of a relation. Throws if `literal`
// resolves to a table or view in the catalog; returns silently otherwise.
//
// Passing a table of accessions by name is the obvious mistake and used to fail
// silently: the name was sent upstream as if it were an accession and came back
// parsing to zero rows (#179). Keying off catalog residency rather than the string's
// shape is deliberate — accession formats are open-ended (bare GenBank ids such as
// J02459 have no distinguishing prefix, see read_ncbi_fasta.cpp), so any pattern
// test would reject valid input.
//
// `function_name` prefixes the message and is quoted back as the remedy, matching
// the convention at read_ncbi.cpp's assembly-accession error.
void RejectRelationNameAsLiteral(ClientContext &context, const std::string &function_name, const std::string &literal);

} // namespace duckdb
