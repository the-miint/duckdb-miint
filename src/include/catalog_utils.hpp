#pragma once

#include "duckdb/common/types.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
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

// Make the caller's TEMP tables and views resolvable on a helper connection (#193).
//
// Every function that takes a relation *by name* reads it through a second
// Connection, because context.Query() deadlocks during bind/execution — see
// docs/internals/reading-tables-views.md. DuckDB constructs
// ClientData::temporary_objects per ClientContext, so that second connection starts
// with its own empty TEMP catalog and cannot see anything the caller created with
// CREATE TEMP TABLE / CREATE TEMP VIEW. The failure is not a clean "unsupported":
// bind-time validation runs Catalog::GetEntry() against the caller's own context,
// which *does* see TEMP objects, so bind passes and only execution fails — with a
// message claiming a table the user can plainly see does not exist.
//
// Assigning the caller's shared_ptr makes both contexts resolve against one
// AttachedDatabase. Sharing the pointer (rather than the object) also keeps the
// catalog alive for as long as either context holds it.
//
// ONLY for connections that read user-named relations and create nothing of their
// own. A CREATE TEMP on an inheriting connection lands in the *caller's* catalog,
// where two hazards bite: parallel workers sharing a fixed temp name collide there
// (silent wrong answers, e.g. the per-sample __<tool>_per_sample views), and the
// object outlives the helper connection instead of dying with it. Sites that create
// TEMP objects must either stay isolated, or uniquify their names AND drop them via
// HelperTempRelation / DropHelperTempRelation below — both halves, not either one.
void InheritTempObjects(ClientContext &context, Connection &conn);

// A helper Connection that already inherits the caller's TEMP catalog, replacing the
// two-line `DatabaseInstance::GetDatabase(context)` + `Connection conn(db)` idiom.
// Carries the same read-only constraint as InheritTempObjects.
Connection MakeReadOnlyHelperConnection(ClientContext &context);

// Drop a TEMP relation that a helper connection created. `quoted_name` must already
// be SQL-quoted; an empty name is a no-op (the relation was never created). `conn`
// must still be live — a moved-from Connection has a null context and is not
// checked for, the same precondition every other Connection& helper here carries.
//
// Safe to call from a destructor: never throws. Unlike a plain fire-and-forget DROP,
// a failure is reported through miint::EmitWarning instead of being discarded,
// because on an inheriting connection the relation lives in the *caller's* session
// catalog — a failed drop leaks an internal relation into it, where the user will
// see it in SHOW TABLES, rather than disappearing with the connection.
void DropHelperTempRelation(Connection &conn, const std::string &quoted_name);

// RAII form of DropHelperTempRelation, for a relation created and finished with
// inside a single function. States that create it and hand it to a longer-lived
// consumer should call DropHelperTempRelation from their own destructor instead
// (see the rype_* global states).
class HelperTempRelation {
public:
	HelperTempRelation(Connection &conn, std::string quoted_name) : conn(conn), quoted_name(std::move(quoted_name)) {
	}
	HelperTempRelation(const HelperTempRelation &) = delete;
	HelperTempRelation &operator=(const HelperTempRelation &) = delete;
	~HelperTempRelation() {
		DropHelperTempRelation(conn, quoted_name);
	}

private:
	Connection &conn;
	std::string quoted_name;
};

} // namespace duckdb
