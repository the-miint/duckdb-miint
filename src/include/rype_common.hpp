#pragma once

#include "rype.h"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "miint_log.hpp"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/enums/arrow_format_version.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/storage/buffer_manager.hpp"

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

//! Floor for the memory budget handed to RYpe. Below this the recommender just
//! returns its own 1000-row minimum, so a smaller number buys nothing.
constexpr size_t RYPE_MIN_MEMORY_BUDGET = 256ULL * 1024 * 1024;

//! Resolve the `max_memory` byte budget to pass to RYpe's batch sizing.
//!
//! `requested` is the caller's `max_memory` named parameter; > 0 is used
//! verbatim, on the principle that someone who names a number knows their
//! deployment better than this heuristic does.
//!
//! Otherwise: RYpe's own auto-detection (`max_memory = 0`) resolves to the
//! cgroups/SLURM limit for the whole process
//! (`ext/rype/src/memory.rs` `detect_available_memory`). Inside DuckDB that
//! double-counts whatever DuckDB is already holding — two budgets drawn against
//! one allocation, neither aware of the other, which is #204.
//!
//! Policy: subtract what DuckDB is *actually* holding, plus a reserve for it to
//! grow into. Deliberately not `memory_limit`: that is a ceiling on the buffer
//! pool rather than a reservation, it defaults to ~80% of RAM, and subtracting
//! it would cut the budget roughly fourfold on an unconfigured host — which,
//! since batch size scales with the budget and each batch costs a full pass over
//! the index, trades an occasional OOM for a reliable severalfold slowdown. A
//! deployment that needs disjoint budgets guaranteed rather than estimated
//! should say so with `max_memory`, which is the point of the parameter.
//!
//! The reserve matches what `rype_extract` already withheld before this existed
//! (a tenth, floor 256 MiB), so no call site becomes more permissive.
//!
//! Returns 0 only when detection fails, which leaves RYpe to its own devices —
//! the same behavior as before this existed, and better than inventing a number
//! from a failed measurement.
inline size_t ResolveRypeMemoryBudget(ClientContext &context, int64_t requested) {
	if (requested > 0) {
		return static_cast<size_t>(requested);
	}
	const size_t available = rype_detect_available_memory();
	if (available == 0) {
		return 0;
	}
	const size_t duckdb_held = BufferManager::GetBufferManager(context).GetUsedMemory();
	const size_t reserve = MaxValue<size_t>(available / 10, RYPE_MIN_MEMORY_BUDGET);
	const size_t claimed = duckdb_held + reserve;
	if (claimed >= available) {
		return RYPE_MIN_MEMORY_BUDGET;
	}
	return MaxValue<size_t>(available - claimed, RYPE_MIN_MEMORY_BUDGET);
}

//! Ask RYpe to size a classification batch, and report the numbers when asked.
//!
//! Wraps rype_calculate_batch_config rather than rype_recommend_batch_size
//! because it returns the same batch_size plus the memory estimates behind it.
//! Those estimates were the missing signal in #204/#199: shard loading is ~99.9%
//! of a classify, so batch count is very nearly a direct multiplier on runtime,
//! and nothing reached the caller. `debug` surfaces them through
//! miint::EmitWarning, i.e. stderr and miint_warnings().
//!
//! `label` names the calling function in the message. `sizing_index` is the
//! index whose shard sizes bound the estimate — for log-ratio, the larger of the
//! two. Returns STANDARD_VECTOR_SIZE if RYpe cannot size a batch, warning either
//! way since that fallback is not a considered choice.
inline size_t ResolveRypeBatchSize(ClientContext &context, const char *label, const RypeIndex *sizing_index,
                                   size_t avg_read_length, int is_paired, size_t memory_budget, bool debug) {
	// is_large_binary=1 tells RYpe to skip its 2 GiB batch cap. That is only sound
	// because ConfigureRypeArrowExport pinned the connection to Arrow LargeBinary
	// (i64 offsets) — the two must be changed together (#222).
	const RypeBatchConfig config =
	    rype_calculate_batch_config(sizing_index, avg_read_length, is_paired, memory_budget, 1);

	if (config.batch_size == 0) {
		const char *err = rype_get_last_error();
		miint::EmitWarning(context, "%s: could not size a classification batch (%s); falling back to %llu reads", label,
		                   err ? err : "unknown error", (unsigned long long)STANDARD_VECTOR_SIZE);
		return STANDARD_VECTOR_SIZE;
	}

	if (debug) {
		// batch_count is documented as "always 1, reserved for
		// forward-compatibility", so it is deliberately not reported.
		miint::EmitWarning(context,
		                   "%s debug: classification batch = %llu reads; memory budget %.2f GB, "
		                   "estimated %.2f GB per batch and %.2f GB peak; avg read length %llu, is_paired %d",
		                   label, (unsigned long long)config.batch_size, double(memory_budget) / 1e9,
		                   double(config.per_batch_memory) / 1e9, double(config.peak_memory) / 1e9,
		                   (unsigned long long)avg_read_length, is_paired);
	}
	return config.batch_size;
}

//! Sample average read length from the first 1000 sequences in a table or view.
//! Returns fallback (default 300) if the query fails or the relation is empty.
//!
//! Runs against the caller's relation. The LIMIT pushes down, so for a view over
//! a join this touches 1000 rows rather than the whole corpus.
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

//! Whether the caller's relation actually carries paired reads: true iff at least
//! one row has a non-NULL sequence2.
//!
//! Callers must use this rather than the mere presence of a `sequence2` column.
//! read_fastx always emits that column (see read_fastx.hpp), so single-end reads
//! loaded the obvious way arrive with it all-NULL. Treating that as paired-end
//! doubles RYpe's per-read memory estimate, which roughly halves the batch size,
//! which doubles the number of full index loads — measured at ~1.8 h of waste on a
//! 4 h job (#199).
//!
//! Callers must ALSO skip this entirely when the relation has no `sequence2`
//! column at all: miint then projects `NULL::BLOB AS sequence2` itself, so the
//! answer is known to be false without touching the data. That is the common
//! single-end shape, and the scan below is the expensive part of it.
//!
//! Deliberately a FULL scan, unlike SampleAvgReadLength's LIMIT 1000: a paired row
//! beyond the sample would be missed, and a false negative under-budgets memory,
//! which is the direction that OOMs. A false positive only costs time. Only the
//! `sequence2` column is projected, so column pruning keeps it off the sequence
//! bytes.
//!
//! Returns false for an empty relation: bool_or over zero rows is NULL, and zero
//! rows genuinely means zero paired rows, so false is the correct answer there
//! rather than a guess.
//!
//! Throws on query failure — it does NOT fall back to false. Falling back would
//! produce exactly the false negative described above, i.e. the OOM direction, and
//! it would do so precisely when memory pressure is the likeliest cause of the
//! failure. Rule 10: surface it.
inline bool TableHasPairedContent(Connection &conn, const std::string &table_quoted) {
	auto result = conn.Query("SELECT bool_or(sequence2 IS NOT NULL) FROM " + table_quoted);
	if (result->HasError()) {
		throw InvalidInputException("Failed to detect paired content in '%s': %s", table_quoted, result->GetError());
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	auto chunk = materialized.Fetch();
	if (!chunk || chunk->size() == 0) {
		return false;
	}
	auto val = chunk->data[0].GetValue(0);
	return !val.IsNull() && val.GetValue<bool>();
}

//! Column names (lowercased) and their storage types for a table/view,
//! index-aligned. The single backing query for both GetTableColumnNamesLower
//! and the id-type capture in ValidateSequenceTable.
struct LowercasedColumns {
	vector<string> names;
	vector<LogicalType> types;
};

//! Look up a table or view by name and return its column names (lowercased) and
//! types, index-aligned. Throws BinderException if the entry does not exist or
//! is not a table/view. `role` names the entry in the "does not exist" message
//! (e.g. "Table or view", "chunk table"), so callers can produce role-specific
//! diagnostics.
inline LowercasedColumns GetTableColumnsLower(ClientContext &context, const std::string &table_name,
                                              const std::string &role = "Table or view") {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("%s '%s' does not exist", role, table_name);
	}

	LowercasedColumns cols;
	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			cols.names.push_back(StringUtil::Lower(col.Name()));
			cols.types.push_back(col.Type());
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		for (idx_t i = 0; i < col_info->names.size(); i++) {
			cols.names.push_back(StringUtil::Lower(col_info->names[i]));
			cols.types.push_back(col_info->types[i]);
		}
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}
	return cols;
}

//! Look up a table or view by name and return its column names, lowercased.
//! Throws BinderException if the entry does not exist or is not a table/view.
//! `role` names the entry in the "does not exist" message (e.g. "Table or view",
//! "chunk table"), so callers can produce role-specific diagnostics.
inline vector<string> GetTableColumnNamesLower(ClientContext &context, const std::string &table_name,
                                               const std::string &role = "Table or view") {
	return GetTableColumnsLower(context, table_name, role).names;
}

//! Result of validating a RYpe sequence table: whether the optional "sequence2"
//! column is present (paired-end, used by rype_classify/rype_log_ratio) and the
//! storage type of the id column.
struct RypeSequenceTableInfo {
	bool has_sequence2 = false;
	LogicalType id_type = LogicalType(LogicalTypeId::INVALID);
};

//! Validate that a table/view exists and has the required columns for RYpe
//! functions. Reports whether the optional "sequence2" column is present and the
//! id column's storage type. All RYpe functions require id_column and
//! "sequence1". The id column must be VARCHAR, BIGINT, or UUID — RYpe carries
//! ids as strings internally and emits them back in their SQL storage type, so
//! only the codec-supported set round-trips losslessly; anything else is
//! rejected loud at bind, matching align_minimap2 and the rest of the id-aware
//! tree.
inline RypeSequenceTableInfo ValidateSequenceTable(ClientContext &context, const std::string &table_name,
                                                   const std::string &id_column) {
	auto cols = GetTableColumnsLower(context, table_name);
	auto &col_names = cols.names;

	auto id_col_lower = StringUtil::Lower(id_column);
	auto id_it = std::find(col_names.begin(), col_names.end(), id_col_lower);
	bool has_seq1 = std::find(col_names.begin(), col_names.end(), "sequence1") != col_names.end();
	bool has_seq2 = std::find(col_names.begin(), col_names.end(), "sequence2") != col_names.end();

	if (id_it == col_names.end()) {
		throw BinderException("Table '%s' missing required column '%s'", table_name, id_column);
	}
	if (!has_seq1) {
		throw BinderException("Table '%s' missing required column 'sequence1'", table_name);
	}

	RypeSequenceTableInfo info;
	info.has_sequence2 = has_seq2;
	info.id_type = cols.types[static_cast<size_t>(std::distance(col_names.begin(), id_it))];

	if (!IsAllowedIdType(info.id_type)) {
		throw BinderException("Column '%s' in table '%s' must be %s", id_column, table_name, AllowedIdTypeList());
	}

	return info;
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
// build (a) a map from a synthetic id to the caller's read_id and (b) an Arrow
// stream of (id, sequence, [pair_sequence]) for RYpe to consume. The naive
// approach of issuing two independent SELECTs against the user's
// sequence_table corrupts the correspondence whenever the two scans see rows
// in different orders (multi-threaded scans, views, parquet sources,
// preserve_insertion_order=false). This used to be fixed by materializing the
// source into a per-call TEMP table; RypeInputStream now does it with a single
// streaming scan that carries the identifier and the sequence in the same row,
// which is the same guarantee without the second copy of the corpus. See
// rype_input_stream.hpp.
//
// Usage pattern (in InitGlobal, on a per-GlobalState sub-Connection):
//
//   ConfigureRypeArrowExport(conn);
//   size_t avg_read_length = SampleAvgReadLength(conn, table_quoted);
//   int is_paired = (has_sequence2 && TableHasPairedContent(conn, table_quoted)) ? 1 : 0;
//   // ... compute batch_size from avg_read_length and is_paired ...
//   auto input = BuildRypeInputStream(conn, gstate->id_map, std::move(opts));
//   ArrowArrayStream *input_stream = &input->stream;
//   // hand input_stream to rype_*_arrow(); on success, input.release().
//
// The destructor must release the RYpe output stream BEFORE resetting
// input_connection: releasing it releases the input stream, which owns the
// streaming QueryResult that holds a non-owning ClientContext pointer into that
// connection.

//! Pin `conn` to exporting variable-length columns with 64-bit offsets (Arrow
//! LargeBinary / LargeUtf8). Call once on the sub-connection, before building any
//! RYpe input stream.
//!
//! Do NOT express this as `arrow_output_version = '1.4'` (BinaryView). BinaryView
//! lifts the limit on the *number* of variadic data buffers, not the 32-bit offset
//! field inside each view, and DuckDB emits exactly one buffer with `buffer_index`
//! hardcoded to 0. That offset is written through `UnsafeNumericCast<int32_t>`,
//! which is an unchecked `static_cast` in release builds, so a record batch
//! carrying more than 2 GiB of sequence silently truncates it modulo 2^32. RYpe
//! reads the field as u32 and is handed a *different read's* bytes under the
//! correct read_id; a spec-conformant signed reader indexes before the buffer
//! entirely. Both boundaries were measured exactly (#222). Batch size is derived
//! from available memory, so a large corpus on a large machine reaches this
//! routinely rather than exceptionally.
//!
//! BOTH settings must be pinned, not just the offset size: for BLOB the BinaryView
//! branch wins over `arrow_large_buffer_size` in ArrowAppender's type dispatch, and
//! both settings are GLOBAL_DEFAULT scope — so a caller who has done
//! `SET arrow_output_version='1.4'` for their own reasons would otherwise have it
//! inherited by this fresh sub-connection and silently reintroduce the corruption.
//!
//! Throws rather than warns. Continuing on failure would fall back to 32-bit
//! offsets, which is precisely the corruption this exists to prevent.
inline void ConfigureRypeArrowExport(Connection &conn) {
	// SET SESSION, not bare SET. Both settings are GLOBAL_DEFAULT scope, meaning a
	// bare SET writes the *global* default — which is shared with the user's own
	// connection. The pre-#222 code used a bare SET and so silently rewrote the
	// caller's global arrow_output_version to '1.4' as a side effect of classifying,
	// changing the format of their own unrelated Arrow exports. Session scope keeps
	// this confined to the sub-connection that feeds RYpe.
	for (const char *stmt :
	     {"SET SESSION arrow_large_buffer_size = true", "SET SESSION arrow_output_version = '1.0'"}) {
		auto result = conn.Query(stmt);
		if (result->HasError()) {
			throw InvalidInputException("Failed to configure Arrow export for RYpe (%s): %s", stmt, result->GetError());
		}
	}

	// Verify what the settings actually resolved to rather than trusting that the
	// statements above mean what we think. This is the assertion that fails loudly
	// if the export format is ever changed out from under RYpe again.
	auto props = conn.context->GetClientProperties();
	if (props.arrow_offset_size != ArrowOffsetSize::LARGE || props.arrow_output_version >= ArrowFormatVersion::V1_4) {
		throw InvalidInputException(
		    "RYpe Arrow export must use 64-bit offsets: expected arrow_large_buffer_size=true and "
		    "arrow_output_version < 1.4, got offset_size=%s output_version=%d",
		    props.arrow_offset_size == ArrowOffsetSize::LARGE ? "LARGE" : "REGULAR",
		    static_cast<int>(props.arrow_output_version));
	}
}

} // namespace duckdb
