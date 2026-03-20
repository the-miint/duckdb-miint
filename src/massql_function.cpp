#include "massql_function.hpp"
#include "massql_parser.hpp"
#include "massql_transpiler.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <atomic>

namespace duckdb {

// ── massql() table function ──────────────────────────────────────────────────
// Usage: SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=220', 'spectra')
//        SELECT * FROM massql('QUERY ...', 'path/to/file.mzML')

struct MassQLData : public TableFunctionData {
	// Non-sample_id: pre-run result at Bind time; ownership moved to GlobalState at InitGlobal.
	unique_ptr<MaterializedQueryResult> non_sample_result;

	// Sample iteration state (sample_id path only)
	bool has_sample_id = false;
	string sample_id_col;
	LogicalType sample_id_type;
	vector<Value> sample_values;

	// Sample_id: sample[0] pre-run for schema inference; moved to GlobalState at InitGlobal.
	unique_ptr<MaterializedQueryResult> sample0_result;

	// Deferred execution state (kept alive for per-sample pipeline in Execute)
	miint::MassQLQuery parsed;
	string effective_source;
};

struct MassQLGlobalState : public GlobalTableFunctionState {
	// Sample_id: atomic counter starts at 1 (sample[0] pre-run at Bind time).
	atomic<idx_t> next_sample_idx {0};
	idx_t max_threads = 1;

	// Non-sample_id: pre-run result transferred from MassQLData at InitGlobal.
	// Single thread drains it from Execute; no synchronization needed.
	unique_ptr<MaterializedQueryResult> non_sample_result;

	// Sample_id: sample[0] pre-run result (schema probe); claimed by the first Execute thread.
	// Each per-thread Connection owns its own TEMP objects (__massql_base, __massql_per_sample,
	// __massql_ms1), so parallel threads never collide on those names despite the shared names.
	unique_ptr<MaterializedQueryResult> sample0_result;
	atomic<bool> sample0_claimed {false};

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

struct MassQLLocalState : public LocalTableFunctionState {
	unique_ptr<Connection> conn;                // sample_id: per-thread isolated connection
	unique_ptr<MaterializedQueryResult> result; // current sample's result
	// Keeps the fetched chunk alive while output.Reference() points into its buffers.
	// Overwritten only on the next Execute call, after the upstream operator has consumed
	// the previous output (DuckDB pull-based execution guarantee).
	unique_ptr<DataChunk> current_chunk;
};

// Materialize base peaks (and optionally MS1 peaks) into temp tables.
// Returns the ms1_table name (empty string if not needed).
static string MaterializePeaks(Connection &conn, const miint::MassQLQuery &parsed, const string &source) {
	auto mat_sql = miint::MassQLTranspiler::materialize_base_sql(parsed, source, "__massql_base");
	auto mat_result = conn.Query(mat_sql);
	if (mat_result->HasError()) {
		throw InvalidInputException("MassQL: failed to materialize peaks: %s", mat_result->GetError());
	}

	auto plan = miint::MassQLTranspiler::get_materialization_plan(parsed);
	string ms1_table;
	if (plan.needs_ms1) {
		ms1_table = "__massql_ms1";
		auto ms1_sql = miint::MassQLTranspiler::materialize_ms1_sql(source, ms1_table);
		auto ms1_result = conn.Query(ms1_sql);
		if (ms1_result->HasError()) {
			throw InvalidInputException("MassQL: failed to materialize MS1 peaks: %s", ms1_result->GetError());
		}
	}

	return ms1_table;
}

// Materialize peaks, transpile, and execute. Returns the materialized result.
static unique_ptr<MaterializedQueryResult> RunPipeline(Connection &conn, const miint::MassQLQuery &parsed,
                                                       const string &source) {
	auto ms1_table = MaterializePeaks(conn, parsed, source);

	string exec_sql;
	try {
		exec_sql = miint::MassQLTranspiler::to_sql_materialized(parsed, source, "__massql_base", ms1_table);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	auto result = conn.Query(exec_sql);
	if (result->HasError()) {
		throw InvalidInputException("MassQL query failed: %s\nGenerated SQL: %s", result->GetError(), exec_sql);
	}
	return result;
}

// Run the MassQL pipeline for a single sample value. Returns the materialized result.
// Each call creates its own TEMP objects (__massql_per_sample, __massql_base, __massql_ms1)
// scoped to `conn`. Since each thread has its own Connection, parallel calls are safe.
static unique_ptr<MaterializedQueryResult> RunSamplePipeline(Connection &conn, const miint::MassQLQuery &parsed,
                                                             const string &effective_source,
                                                             const string &sample_id_col, const Value &sample_value) {
	auto quoted_col = KeywordHelper::WriteOptionallyQuoted(sample_id_col);
	auto quoted_source = KeywordHelper::WriteOptionallyQuoted(effective_source);

	// Use ToSQLString() for safe SQL literal construction — handles all Value types
	// (integers, timestamps, intervals, etc.) without manual escaping.
	auto filter_val = sample_value.ToSQLString();
	auto view_sql = "CREATE OR REPLACE TEMP VIEW __massql_per_sample AS SELECT * FROM " + quoted_source +
	                " WHERE CAST(" + quoted_col + " AS VARCHAR) = CAST(" + filter_val + " AS VARCHAR)";
	auto view_result = conn.Query(view_sql);
	if (view_result->HasError()) {
		throw InvalidInputException("MassQL: failed to create per-sample view: %s", view_result->GetError());
	}

	// Materialize peaks from the filtered view
	auto ms1_table = MaterializePeaks(conn, parsed, "__massql_per_sample");

	// Generate transpiled SQL and wrap with sample_id column (single execution)
	string exec_sql;
	try {
		exec_sql =
		    miint::MassQLTranspiler::to_sql_materialized(parsed, "__massql_per_sample", "__massql_base", ms1_table);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	auto sample_literal = sample_value.ToSQLString();
	auto wrapped_sql = "SELECT " + sample_literal + " AS " + quoted_col + ", __q.* FROM (" + exec_sql + ") __q";

	auto result = conn.Query(wrapped_sql);
	if (result->HasError()) {
		throw InvalidInputException("MassQL query failed: %s\nGenerated SQL: %s", result->GetError(), wrapped_sql);
	}
	return result;
}

static void ExtractSchema(MaterializedQueryResult &result, vector<LogicalType> &return_types, vector<string> &names) {
	for (idx_t i = 0; i < result.ColumnCount(); i++) {
		names.push_back(result.ColumnName(i));
		return_types.push_back(result.types[i]);
	}
}

static unique_ptr<FunctionData> MassQLBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MassQLData>();

	auto query_str = input.inputs[0].GetValue<string>();
	auto source_str = input.inputs[1].GetValue<string>();

	// Read optional sample_id named parameter
	string sample_id_col;
	bool has_sample_id = false;
	auto it = input.named_parameters.find("sample_id");
	if (it != input.named_parameters.end()) {
		sample_id_col = it->second.GetValue<string>();
		if (sample_id_col.empty()) {
			throw InvalidInputException("sample_id column name must not be empty");
		}
		has_sample_id = true;
	}

	// Parse MassQL (convert parser exceptions to DuckDB exceptions)
	miint::MassQLQuery parsed;
	try {
		parsed = miint::MassQLParser::parse(query_str);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	auto &db = DatabaseInstance::GetDatabase(context);

	// Reject sample_id with file path source early
	auto effective_source = source_str;
	if (miint::MassQLTranspiler::is_file_path(source_str)) {
		if (has_sample_id) {
			throw InvalidInputException("sample_id cannot be used with a file path source");
		}
	}

	if (has_sample_id) {
		// ── sample_id path: use stack-local conn for bind-time validation only ──
		Connection conn(db);

		// Validate column exists and get its type
		auto quoted_col = KeywordHelper::WriteOptionallyQuoted(sample_id_col);
		auto quoted_source = KeywordHelper::WriteOptionallyQuoted(effective_source);
		auto probe = conn.Query("SELECT " + quoted_col + " FROM " + quoted_source + " LIMIT 0");
		if (probe->HasError()) {
			throw InvalidInputException("sample_id column '%s' not found", sample_id_col);
		}
		data->sample_id_type = probe->types[0];

		// Reject NULLs early before collecting distinct values
		auto null_check = conn.Query("SELECT COUNT(*) FROM " + quoted_source + " WHERE " + quoted_col + " IS NULL");
		if (null_check->HasError()) {
			throw InvalidInputException("MassQL: failed to check for NULLs: %s", null_check->GetError());
		}
		auto null_count = null_check->GetValue(0, 0).GetValue<int64_t>();
		if (null_count > 0) {
			throw InvalidInputException("NULL values in sample_id column '%s'", sample_id_col);
		}

		// Collect distinct values
		auto distinct_result =
		    conn.Query("SELECT DISTINCT " + quoted_col + " FROM " + quoted_source + " ORDER BY " + quoted_col);
		if (distinct_result->HasError()) {
			throw InvalidInputException("MassQL: failed to query sample_id values: %s", distinct_result->GetError());
		}
		auto &materialized = distinct_result->Cast<MaterializedQueryResult>();
		while (auto chunk = materialized.Fetch()) {
			for (idx_t i = 0; i < chunk->size(); i++) {
				data->sample_values.push_back(chunk->data[0].GetValue(i));
			}
		}

		if (data->sample_values.empty()) {
			throw InvalidInputException("sample_id column '%s' has no non-NULL values", sample_id_col);
		}

		// Store state for deferred execution
		data->has_sample_id = true;
		data->sample_id_col = sample_id_col;
		data->parsed = parsed;
		data->effective_source = effective_source;

		// Run sample[0] to determine the output schema. The result is stored here and
		// moved to GlobalState at InitGlobal so it is not re-run in Execute.
		data->sample0_result = RunSamplePipeline(conn, parsed, effective_source, sample_id_col, data->sample_values[0]);
		ExtractSchema(*data->sample0_result, return_types, names);
	} else {
		// ── non-sample_id path: existing behavior, conn is stack-local ──
		Connection conn(db);

		if (miint::MassQLTranspiler::is_file_path(source_str)) {
			effective_source = "__massql_src";
			auto escaped_path = StringUtil::Replace(source_str, "'", "''");
			auto view_result = conn.Query("CREATE OR REPLACE TEMP VIEW " + effective_source +
			                              " AS SELECT * FROM read_mzml('" + escaped_path + "')");
			if (view_result->HasError()) {
				throw InvalidInputException("MassQL: failed to read source file: %s", view_result->GetError());
			}
		}

		// Run pipeline at bind time. Result is fully materialized before conn goes out of scope,
		// so any TEMP views created above are no longer needed after this call.
		data->non_sample_result = RunPipeline(conn, parsed, effective_source);
		ExtractSchema(*data->non_sample_result, return_types, names);
	}

	return data;
}

static unique_ptr<GlobalTableFunctionState> MassQLInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	// CastNoConst: we move result ownership out of bind data into global state.
	// InitGlobal is called exactly once before any Execute thread starts, so this is safe.
	auto &data = input.bind_data->CastNoConst<MassQLData>();
	auto gstate = make_uniq<MassQLGlobalState>();
	if (data.has_sample_id) {
		idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		gstate->max_threads = std::max<idx_t>(1, std::min(db_threads, data.sample_values.size()));
		// Transfer sample[0]'s pre-run result; Execute threads start claiming from index 1.
		gstate->sample0_result = std::move(data.sample0_result);
		gstate->next_sample_idx = 1;
	} else {
		// Transfer the pre-run result; single Execute thread drains it from global state.
		gstate->non_sample_result = std::move(data.non_sample_result);
	}
	return gstate;
}

static unique_ptr<LocalTableFunctionState> MassQLInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                           GlobalTableFunctionState * /*global_state*/) {
	auto &data = input.bind_data->Cast<MassQLData>();
	auto lstate = make_uniq<MassQLLocalState>();
	if (data.has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context.client);
		lstate->conn = make_uniq<Connection>(db);
	}
	return lstate;
}

static void MassQLExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &data = input.bind_data->Cast<MassQLData>();
	auto &gstate = input.global_state->Cast<MassQLGlobalState>();
	auto &lstate = input.local_state->Cast<MassQLLocalState>();

	if (!data.has_sample_id) {
		// Single thread drains the pre-run result from global state.
		// current_chunk keeps the chunk alive while output references its buffers;
		// it is overwritten only on the next Execute call (DuckDB pull-based guarantee).
		lstate.current_chunk = gstate.non_sample_result->Fetch();
		if (lstate.current_chunk && lstate.current_chunk->size() > 0) {
			output.Reference(*lstate.current_chunk);
			return;
		}
		output.SetCardinality(0);
		return;
	}

	// Sample_id: each thread atomically claims sample indices and processes independently.
	while (true) {
		if (lstate.result) {
			lstate.current_chunk = lstate.result->Fetch();
			if (lstate.current_chunk && lstate.current_chunk->size() > 0) {
				output.Reference(*lstate.current_chunk);
				return;
			}
			lstate.result.reset();
		}
		// Claim sample[0]'s pre-run result before pulling new samples from the counter.
		// exchange returns the old value: false means this thread is the claimant.
		if (!gstate.sample0_claimed.exchange(true, std::memory_order_acq_rel)) {
			lstate.result = std::move(gstate.sample0_result);
			continue;
		}
		// Atomically claim the next sample index. relaxed: only atomicity needed, no
		// ordering with respect to other memory operations.
		idx_t sample_idx = gstate.next_sample_idx.fetch_add(1, std::memory_order_relaxed);
		if (sample_idx >= data.sample_values.size()) {
			output.SetCardinality(0);
			return;
		}
		lstate.result = RunSamplePipeline(*lstate.conn, data.parsed, data.effective_source, data.sample_id_col,
		                                  data.sample_values[sample_idx]);
	}
}

// ── massql_to_sql() scalar function ──────────────────────────────────────────
// Returns the generated SQL for debugging: SELECT massql_to_sql('QUERY ...', 'source')

static void MassQLToSQLFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	BinaryExecutor::Execute<string_t, string_t, string_t>(
	    args.data[0], args.data[1], result, args.size(), [&](string_t query_str, string_t source_str) {
		    try {
			    auto parsed = miint::MassQLParser::parse(query_str.GetString());
			    auto sql = miint::MassQLTranspiler::to_sql(parsed, source_str.GetString());
			    return StringVector::AddString(result, sql);
		    } catch (const std::exception &e) {
			    throw InvalidInputException(e.what());
		    }
	    });
}

// ── Registration ─────────────────────────────────────────────────────────────

void MassQLFunction::Register(ExtensionLoader &loader) {
	// massql(query, source) table function
	// order_preservation_type=NO_ORDER: parallel samples produce non-deterministic interleaving.
	TableFunction massql_func("massql", {LogicalType::VARCHAR, LogicalType::VARCHAR}, MassQLExecute, MassQLBind,
	                          MassQLInitGlobal, MassQLInitLocal);
	massql_func.named_parameters["sample_id"] = LogicalType::VARCHAR;
	massql_func.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(massql_func);

	// massql_to_sql(query, source) scalar function
	ScalarFunction to_sql_func("massql_to_sql", {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::VARCHAR,
	                           MassQLToSQLFunction);
	loader.RegisterFunction(to_sql_func);
}

} // namespace duckdb
