#include "massql_function.hpp"
#include "massql_parser.hpp"
#include "massql_transpiler.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector_operations/binary_executor.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

// ── massql() table function ──────────────────────────────────────────────────
// Usage: SELECT * FROM massql('QUERY scaninfo(MS2DATA) WHERE MS2PROD=220', 'spectra')
//        SELECT * FROM massql('QUERY ...', 'path/to/file.mzML')

struct MassQLData : public TableFunctionData {
	// Result iteration (both paths)
	unique_ptr<MaterializedQueryResult> result;
	unique_ptr<DataChunk> current_chunk; // keeps fetched chunk alive for Reference()

	// Sample iteration state (sample_id path only)
	bool has_sample_id = false;
	string sample_id_col;
	LogicalType sample_id_type;
	vector<Value> sample_values;
	idx_t current_sample_idx = 0;

	// Deferred execution state (kept alive for per-sample pipeline in Execute)
	miint::MassQLQuery parsed;
	string effective_source;
	unique_ptr<Connection> conn;
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

// Materialize peaks, transpile, execute, and store result in data.result.
static void RunPipeline(Connection &conn, const miint::MassQLQuery &parsed, const string &source,
                        MassQLData &data) {
	auto ms1_table = MaterializePeaks(conn, parsed, source);

	string exec_sql;
	try {
		exec_sql = miint::MassQLTranspiler::to_sql_materialized(parsed, source, "__massql_base", ms1_table);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	data.result = conn.Query(exec_sql);
	if (data.result->HasError()) {
		throw InvalidInputException("MassQL query failed: %s\nGenerated SQL: %s", data.result->GetError(), exec_sql);
	}
}

// Run the MassQL pipeline for a single sample value: create a filtered view,
// materialize peaks, and execute a single wrapped query with the sample_id column.
static void RunSamplePipeline(MassQLData &data, const Value &sample_value) {
	auto &conn = *data.conn;
	auto quoted_col = KeywordHelper::WriteOptionallyQuoted(data.sample_id_col);
	auto quoted_source = KeywordHelper::WriteOptionallyQuoted(data.effective_source);

	// Create filtered view for this sample
	auto escaped_val = StringUtil::Replace(sample_value.ToString(), "'", "''");
	auto view_sql = "CREATE OR REPLACE TEMP VIEW __massql_per_sample AS SELECT * FROM " + quoted_source +
	                " WHERE CAST(" + quoted_col + " AS VARCHAR) = '" + escaped_val + "'";
	auto view_result = conn.Query(view_sql);
	if (view_result->HasError()) {
		throw InvalidInputException("MassQL: failed to create per-sample view: %s", view_result->GetError());
	}

	// Materialize peaks from the filtered view
	auto ms1_table = MaterializePeaks(conn, data.parsed, "__massql_per_sample");

	// Generate transpiled SQL and wrap with sample_id column (single execution)
	string exec_sql;
	try {
		exec_sql =
		    miint::MassQLTranspiler::to_sql_materialized(data.parsed, "__massql_per_sample", "__massql_base", ms1_table);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	auto sample_literal = sample_value.ToSQLString();
	auto wrapped_sql = "SELECT " + sample_literal + " AS " + quoted_col + ", __q.* FROM (" + exec_sql + ") __q";

	data.result = conn.Query(wrapped_sql);
	if (data.result->HasError()) {
		throw InvalidInputException("MassQL query failed: %s\nGenerated SQL: %s", data.result->GetError(), wrapped_sql);
	}
}

static void ExtractSchema(MassQLData &data, vector<LogicalType> &return_types, vector<string> &names) {
	for (idx_t i = 0; i < data.result->ColumnCount(); i++) {
		names.push_back(data.result->ColumnName(i));
		return_types.push_back(data.result->types[i]);
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
		// ── sample_id path: store conn and state for per-sample iteration in Execute ──
		data->conn = make_uniq<Connection>(db);
		auto &conn = *data->conn;

		// Validate column exists and get its type
		auto quoted_col = KeywordHelper::WriteOptionallyQuoted(sample_id_col);
		auto quoted_source = KeywordHelper::WriteOptionallyQuoted(effective_source);
		auto probe = conn.Query("SELECT " + quoted_col + " FROM " + quoted_source + " LIMIT 0");
		if (probe->HasError()) {
			throw InvalidInputException("sample_id column '%s' not found", sample_id_col);
		}
		data->sample_id_type = probe->types[0];

		// Reject NULLs early before collecting distinct values
		auto null_check =
		    conn.Query("SELECT COUNT(*) FROM " + quoted_source + " WHERE " + quoted_col + " IS NULL");
		if (null_check->HasError()) {
			throw InvalidInputException("MassQL: failed to check for NULLs: %s", null_check->GetError());
		}
		auto null_count = null_check->GetValue(0, 0).GetValue<int64_t>();
		if (null_count > 0) {
			throw InvalidInputException("NULL values in sample_id column '%s'", sample_id_col);
		}

		// Collect distinct values
		auto distinct_result = conn.Query("SELECT DISTINCT " + quoted_col + " FROM " + quoted_source + " ORDER BY " +
		                                  quoted_col);
		if (distinct_result->HasError()) {
			throw InvalidInputException("MassQL: failed to query sample_id values: %s", distinct_result->GetError());
		}
		auto &materialized = distinct_result->Cast<MaterializedQueryResult>();
		while (auto chunk = materialized.Fetch()) {
			for (idx_t i = 0; i < chunk->size(); i++) {
				data->sample_values.push_back(chunk->data[0].GetValue(i));
			}
		}

		// Store state for deferred execution
		data->has_sample_id = true;
		data->sample_id_col = sample_id_col;
		data->parsed = parsed;
		data->effective_source = effective_source;

		if (data->sample_values.empty()) {
			throw InvalidInputException("sample_id column '%s' has no non-NULL values", sample_id_col);
		}

		// Run first sample to determine schema
		RunSamplePipeline(*data, data->sample_values[0]);
		data->current_sample_idx = 1;
		ExtractSchema(*data, return_types, names);
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

		RunPipeline(conn, parsed, effective_source, *data);
		ExtractSchema(*data, return_types, names);
	}

	return data;
}

static void MassQLExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &data = input.bind_data->CastNoConst<MassQLData>();

	while (true) {
		// Drain current result
		data.current_chunk = data.result->Fetch();
		if (data.current_chunk && data.current_chunk->size() > 0) {
			output.Reference(*data.current_chunk);
			return;
		}

		// For non-sample_id path, we're done
		if (!data.has_sample_id) {
			output.SetCardinality(0);
			return;
		}

		// All samples exhausted
		if (data.current_sample_idx >= data.sample_values.size()) {
			output.SetCardinality(0);
			return;
		}

		// Run pipeline for next sample
		RunSamplePipeline(data, data.sample_values[data.current_sample_idx]);
		data.current_sample_idx++;
		// Loop back to drain this sample's result
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
	TableFunction massql_func("massql", {LogicalType::VARCHAR, LogicalType::VARCHAR}, MassQLExecute, MassQLBind);
	massql_func.named_parameters["sample_id"] = LogicalType::VARCHAR;
	loader.RegisterFunction(massql_func);

	// massql_to_sql(query, source) scalar function
	ScalarFunction to_sql_func("massql_to_sql", {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::VARCHAR,
	                           MassQLToSQLFunction);
	loader.RegisterFunction(to_sql_func);
}

} // namespace duckdb
