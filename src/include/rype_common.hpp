#pragma once

#include "rype.h"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/arrow/result_arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <string>
#include <vector>

namespace duckdb {

// ============================================================================
// Shared Arrow local state for all RYpe table functions.
// Manages per-column ArrowArrayScanState for zero-copy Arrow→DuckDB conversion.
// ============================================================================
struct RypeArrowLocalState : public LocalTableFunctionState {
	std::unordered_map<idx_t, unique_ptr<ArrowArrayScanState>> array_states;
	ClientContext &context;

	explicit RypeArrowLocalState(ClientContext &ctx) : context(ctx) {
	}

	~RypeArrowLocalState() {
		array_states.clear();
	}

	ArrowArrayScanState &GetState(idx_t col_idx) {
		auto it = array_states.find(col_idx);
		if (it == array_states.end()) {
			auto state = make_uniq<ArrowArrayScanState>(context);
			auto &ref = *state;
			array_states.emplace(col_idx, std::move(state));
			return ref;
		}
		return *it->second;
	}

	void ResetStates() {
		for (auto &state : array_states) {
			state.second->Reset();
		}
	}
};

//! Sample average read length from the first 1000 sequences in a table.
//! Returns fallback (default 300) if the query fails or the table is empty.
inline size_t SampleAvgReadLength(Connection &conn, const std::string &table_quoted, size_t fallback = 300) {
	std::string query =
	    "SELECT AVG(LENGTH(sequence1))::BIGINT FROM (SELECT sequence1 FROM " + table_quoted + " LIMIT 1000)";
	auto result = conn.Query(query);
	if (!result->HasError()) {
		auto &materialized = result->Cast<MaterializedQueryResult>();
		auto chunk = materialized.Fetch();
		if (chunk && chunk->size() > 0 && !chunk->data[0].GetValue(0).IsNull()) {
			auto val = chunk->data[0].GetValue(0).GetValue<int64_t>();
			if (val > 0) {
				return static_cast<size_t>(val);
			}
		}
	}
	return fallback;
}

//! Look up a table or view by name and return its column names, lowercased.
//! Throws BinderException if the entry does not exist or is not a table/view.
//! `role` names the entry in the "does not exist" message (e.g. "Table or view",
//! "chunk table"), so callers can produce role-specific diagnostics.
inline vector<string> GetTableColumnNamesLower(ClientContext &context, const std::string &table_name,
                                               const std::string &role = "Table or view") {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("%s '%s' does not exist", role, table_name);
	}

	vector<string> col_names;
	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			col_names.push_back(StringUtil::Lower(columns.GetColumn(LogicalIndex(i)).Name()));
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		for (const auto &name : col_info->names) {
			col_names.push_back(StringUtil::Lower(name));
		}
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}
	return col_names;
}

//! Validate that a table/view exists and has the required columns for RYpe functions.
//! Returns true if the optional "sequence2" column is present (used by rype_classify
//! for paired-end reads). All RYpe functions require id_column and "sequence1".
inline bool ValidateSequenceTable(ClientContext &context, const std::string &table_name, const std::string &id_column) {
	auto col_names = GetTableColumnNamesLower(context, table_name);

	auto id_col_lower = StringUtil::Lower(id_column);
	bool has_id = std::find(col_names.begin(), col_names.end(), id_col_lower) != col_names.end();
	bool has_seq1 = std::find(col_names.begin(), col_names.end(), "sequence1") != col_names.end();
	bool has_seq2 = std::find(col_names.begin(), col_names.end(), "sequence2") != col_names.end();

	if (!has_id) {
		throw BinderException("Table '%s' missing required column '%s'", table_name, id_column);
	}
	if (!has_seq1) {
		throw BinderException("Table '%s' missing required column 'sequence1'", table_name);
	}

	return has_seq2;
}

//! Validate that a table/view exists and contains every column in
//! `required_columns` (case-insensitive). `role` names the table in error
//! messages (e.g. "chunk table", "mapping table"). Throws BinderException on a
//! missing table or a missing required column.
inline void ValidateTableHasColumns(ClientContext &context, const std::string &table_name,
                                    const std::vector<std::string> &required_columns, const std::string &role) {
	auto col_names = GetTableColumnNamesLower(context, table_name, role);
	for (const auto &req : required_columns) {
		if (std::find(col_names.begin(), col_names.end(), StringUtil::Lower(req)) == col_names.end()) {
			throw BinderException("%s '%s' is missing required column '%s'", role, table_name, req);
		}
	}
}

// ============================================================================
// Shared input pipeline for RYpe table functions (the row-index ↔ read_id fix).
// ============================================================================
//
// Background: rype_classify, rype_extract_*, and rype_log_ratio all need to
// build (a) a vector of read_ids indexed by a synthetic id and (b) an Arrow
// stream of (id, sequence, [pair_sequence]) for RYpe to consume. The naive
// approach of issuing two independent SELECTs against the user's
// sequence_table corrupts the correspondence whenever the two scans see rows
// in different orders (multi-threaded scans, views, parquet sources,
// preserve_insertion_order=false). The helpers below materialize the source
// once into a per-call TEMP table with an explicit id column, then read
// read_ids and stream sequences from that same table ORDER BY id — the id
// column is now a stable attribute of the data, not an emergent property of
// independent scans.
//
// Usage pattern (in InitGlobal, on a per-GlobalState sub-Connection):
//
//   gstate->tmp_table_name = MaterializeRypeInputTempTable(
//       conn, table_quoted, id_col_quoted, source_name_for_errors,
//       has_sequence2, "_rype_classify_", gstate->read_ids,
//       avg_read_length);
//   // ... compute batch_size from avg_read_length ...
//   auto wrapper = BuildRypeArrowInput(conn, gstate->tmp_table_name,
//                                      /*include_pair_column=*/true,
//                                      batch_size);
//   ArrowArrayStream *input_stream = &wrapper->stream;
//   // hand input_stream to rype_*_arrow(); on success, wrapper.release().
//
// The destructor must call DropRypeTempTable(*input_connection,
// gstate->tmp_table_name) AFTER releasing the RYpe output stream and BEFORE
// resetting input_connection.

//! Materialize the user's sequence_table into a per-call TEMP table on `conn`,
//! populate `out_read_ids` from it ordered by the synthetic id, and sample the
//! average read length for batch-size estimation. Returns the name of the
//! created TEMP table; the caller must store this and drop it later via
//! DropRypeTempTable.
//!
//! `source_name_for_errors` is the user-facing source-table name (unquoted);
//! used only in error messages.
//! `name_prefix` is the per-function debug prefix (e.g. "_rype_classify_").
inline std::string MaterializeRypeInputTempTable(Connection &conn, const std::string &table_quoted,
                                                 const std::string &id_col_quoted,
                                                 const std::string &source_name_for_errors, bool has_sequence2,
                                                 const std::string &name_prefix, std::vector<std::string> &out_read_ids,
                                                 size_t &out_avg_read_length) {
	std::string tmp_table_name = name_prefix + StringUtil::Replace(UUID::ToString(UUID::GenerateRandomUUID()), "-", "");
	std::string tmp_quoted = KeywordHelper::WriteOptionallyQuoted(tmp_table_name);
	std::string seq2_proj = has_sequence2 ? "sequence2" : "NULL::BLOB AS sequence2";
	std::string create_sql = "CREATE TEMP TABLE " + tmp_quoted +
	                         " AS SELECT (row_number() OVER () - 1)::BIGINT AS id, " + id_col_quoted +
	                         " AS read_id, sequence1, " + seq2_proj + " FROM " + table_quoted;
	auto create_result = conn.Query(create_sql);
	if (create_result->HasError()) {
		throw InvalidInputException("Failed to materialize sequence table '%s': %s", source_name_for_errors,
		                            create_result->GetError());
	}

	std::string id_query = "SELECT read_id FROM " + tmp_quoted + " ORDER BY id";
	auto id_result = conn.Query(id_query);
	if (id_result->HasError()) {
		throw InvalidInputException("Failed to read read_ids from temp table: %s", id_result->GetError());
	}

	out_read_ids.reserve(id_result->RowCount());
	auto &id_materialized = id_result->Cast<MaterializedQueryResult>();
	while (auto chunk = id_materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			out_read_ids.push_back(chunk->data[0].GetValue(i).ToString());
		}
	}

	out_avg_read_length = SampleAvgReadLength(conn, tmp_quoted);
	return tmp_table_name;
}

//! Build the Arrow input stream RYpe will consume. Reads (id, sequence1, [sequence2])
//! from the named TEMP table ordered by id. `include_pair_column` exposes a
//! pair_sequence column (true for classify/log_ratio, false for extract).
//!
//! Caller transfers ownership of the returned wrapper to RYpe by calling
//! .release() AFTER rype_*_arrow() succeeds; on failure, the unique_ptr's
//! destructor cleans up the wrapper.
inline unique_ptr<ResultArrowArrayStreamWrapper>
BuildRypeArrowInput(Connection &conn, const std::string &tmp_table_name, bool include_pair_column, size_t batch_size) {
	std::string tmp_quoted = KeywordHelper::WriteOptionallyQuoted(tmp_table_name);
	std::string select_cols = include_pair_column
	                              ? std::string("id, sequence1::BLOB AS sequence, sequence2::BLOB AS pair_sequence")
	                              : std::string("id, sequence1::BLOB AS sequence");
	std::string query = "SELECT " + select_cols + " FROM " + tmp_quoted + " ORDER BY id";
	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from temp table: %s", query_result->GetError());
	}
	return make_uniq<ResultArrowArrayStreamWrapper>(std::move(query_result), batch_size);
}

//! Drop the per-call TEMP table on `conn`. Safe with empty name (no-op) and
//! with a name that doesn't exist (uses IF EXISTS). Errors are silently
//! ignored — this runs in destructors where we cannot usefully propagate
//! failures, and the connection's catalog cleanup will reap the table on
//! teardown anyway.
inline void DropRypeTempTable(Connection &conn, const std::string &tmp_table_name) {
	if (tmp_table_name.empty()) {
		return;
	}
	conn.Query("DROP TABLE IF EXISTS " + KeywordHelper::WriteOptionallyQuoted(tmp_table_name));
}

} // namespace duckdb
