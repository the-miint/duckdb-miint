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
	unique_ptr<MaterializedQueryResult> result;
	unique_ptr<DataChunk> current_chunk; // keeps fetched chunk alive for Reference()
};

static unique_ptr<FunctionData> MassQLBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MassQLData>();

	auto query_str = input.inputs[0].GetValue<string>();
	auto source_str = input.inputs[1].GetValue<string>();

	// Parse MassQL (convert parser exceptions to DuckDB exceptions)
	miint::MassQLQuery parsed;
	try {
		parsed = miint::MassQLParser::parse(query_str);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// If source is a file path, create a temp view so mzml_peaks() can reference it by name
	auto effective_source = source_str;
	if (miint::MassQLTranspiler::is_file_path(source_str)) {
		effective_source = "__massql_src";
		auto escaped_path = StringUtil::Replace(source_str, "'", "''");
		auto view_result = conn.Query("CREATE OR REPLACE TEMP VIEW " + effective_source +
		                              " AS SELECT * FROM read_mzml('" + escaped_path + "')");
		if (view_result->HasError()) {
			throw InvalidInputException("MassQL: failed to read source file: %s", view_result->GetError());
		}
	}

	// Materialize peaks into a temp table to avoid repeated mzml_peaks() UNNEST evaluation
	auto mat_sql = miint::MassQLTranspiler::materialize_base_sql(parsed, effective_source, "__massql_base");
	auto mat_result = conn.Query(mat_sql);
	if (mat_result->HasError()) {
		throw InvalidInputException("MassQL: failed to materialize peaks: %s", mat_result->GetError());
	}

	// For cross-level queries (MS1MZ=X AND MS2PREC=X), also materialize MS1 peaks.
	// get_materialization_plan() is the authoritative source for which tables are needed.
	auto plan = miint::MassQLTranspiler::get_materialization_plan(parsed);
	string ms1_table;
	if (plan.needs_ms1) {
		ms1_table = "__massql_ms1";
		auto ms1_sql = miint::MassQLTranspiler::materialize_ms1_sql(effective_source, ms1_table);
		auto ms1_result = conn.Query(ms1_sql);
		if (ms1_result->HasError()) {
			throw InvalidInputException("MassQL: failed to materialize MS1 peaks: %s", ms1_result->GetError());
		}
	}

	// Execute using pre-materialized tables for performance
	string exec_sql;
	try {
		exec_sql = miint::MassQLTranspiler::to_sql_materialized(parsed, effective_source, "__massql_base", ms1_table);
	} catch (const std::exception &e) {
		throw InvalidInputException(e.what());
	}

	data->result = conn.Query(exec_sql);

	if (data->result->HasError()) {
		throw InvalidInputException("MassQL query failed: %s\nGenerated SQL: %s", data->result->GetError(), exec_sql);
	}

	// Extract schema from result
	for (idx_t i = 0; i < data->result->ColumnCount(); i++) {
		names.push_back(data->result->ColumnName(i));
		return_types.push_back(data->result->types[i]);
	}

	return data;
}

static void MassQLExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &data = input.bind_data->CastNoConst<MassQLData>();

	data.current_chunk = data.result->Fetch();
	if (!data.current_chunk || data.current_chunk->size() == 0) {
		output.SetCardinality(0);
		return;
	}

	output.Reference(*data.current_chunk);
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
	loader.RegisterFunction(massql_func);

	// massql_to_sql(query, source) scalar function
	ScalarFunction to_sql_func("massql_to_sql", {LogicalType::VARCHAR, LogicalType::VARCHAR}, LogicalType::VARCHAR,
	                           MassQLToSQLFunction);
	loader.RegisterFunction(to_sql_func);
}

} // namespace duckdb
