#include "rype_index_create.hpp"
#include "rype_common.hpp"

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

// ============================================================================
// Bind
// ============================================================================
unique_ptr<FunctionData> RypeIndexCreateTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                            vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<Data>();

	// Required positional parameters: chunk_table, output_path
	if (input.inputs.size() < 2) {
		throw BinderException("rype_index_create requires chunk_table and output_path parameters");
	}
	data->chunk_table = input.inputs[0].ToString();
	data->output_path = input.inputs[1].ToString();

	// Optional: mapping_table (feature_idx, bucket_name). Empty => single bucket.
	auto mapping_param = input.named_parameters.find("mapping_table");
	if (mapping_param != input.named_parameters.end() && !mapping_param->second.IsNull()) {
		data->mapping_table = mapping_param->second.ToString();
	}

	auto k_param = input.named_parameters.find("k");
	if (k_param != input.named_parameters.end() && !k_param->second.IsNull()) {
		data->k = k_param->second.GetValue<int32_t>();
	}

	auto w_param = input.named_parameters.find("w");
	if (w_param != input.named_parameters.end() && !w_param->second.IsNull()) {
		data->w = w_param->second.GetValue<int32_t>();
	}

	auto salt_param = input.named_parameters.find("salt");
	if (salt_param != input.named_parameters.end() && !salt_param->second.IsNull()) {
		data->salt = salt_param->second.GetValue<uint64_t>();
	}

	auto orient_param = input.named_parameters.find("orient");
	if (orient_param != input.named_parameters.end() && !orient_param->second.IsNull()) {
		data->orient = orient_param->second.GetValue<bool>();
	}

	auto mem_param = input.named_parameters.find("max_memory");
	if (mem_param != input.named_parameters.end() && !mem_param->second.IsNull()) {
		data->max_memory = mem_param->second.GetValue<int64_t>();
	}

	// Validate at bind time (fail fast, before any build side effects). The chunk
	// table is checked first so a bad chunk table is not misattributed to the
	// synthesized mapping query (which is derived from it).
	ValidateTableHasColumns(context, data->chunk_table, {"feature_idx", "chunk_index", "chunk_data"}, "chunk table");
	if (!data->mapping_table.empty()) {
		ValidateTableHasColumns(context, data->mapping_table, {"feature_idx", "bucket_name"}, "mapping table");
	}

	// k must be one of RYpe's supported RY-space k-mer sizes. RYpe also validates
	// this, but checking here fails before the (non-atomic) build writes anything.
	if (data->k != 16 && data->k != 32 && data->k != 64) {
		throw BinderException("rype_index_create: k must be 16, 32, or 64 (got %d)", data->k);
	}

	// w is the minimizer window size and must be >= 1. RYpe does NOT validate this:
	// w <= 0 silently builds a minimizer-free, useless index that still reports
	// success, and a negative w would wrap to a huge size_t via the FFI cast.
	if (data->w < 1) {
		throw BinderException("rype_index_create: w must be >= 1 (got %d)", data->w);
	}

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ============================================================================
// InitGlobal — performs the build (a synchronous side effect)
// ============================================================================
unique_ptr<GlobalTableFunctionState> RypeIndexCreateTableFunction::InitGlobal(ClientContext &context,
                                                                              TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Sub-connection to read the input tables as Arrow (avoids the context.Query
	// re-entrancy deadlock). Held in GlobalState so it outlives RYpe's lazy
	// consumption of the streamed chunk cursor during the synchronous build.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->input_connection = make_uniq<Connection>(db);
	auto &conn = *gstate->input_connection;

	// Arrow v1.4 (Utf8View/BinaryView) for VARCHAR — no i32 offset 2 GiB cap. The
	// rype build path accepts Utf8/LargeUtf8/Binary/LargeBinary/Utf8View/BinaryView.
	conn.Query("SET arrow_output_version = '1.4'");

	std::string chunk_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.chunk_table);

	// Mapping stream: feature_idx Int64, bucket_name Utf8. Materialized — it is
	// small (one row per feature) and RYpe reads the whole mapping into memory
	// before processing chunks, so streaming it would buy nothing. Done BEFORE the
	// streaming cursor below so it does not occupy the connection's stream slot.
	// When no mapping_table is supplied, every feature goes into a single bucket
	// named 'unnamed-bucket', synthesized from the distinct feature_idx in the
	// chunk table.
	std::string mapping_sql;
	if (bind_data.mapping_table.empty()) {
		mapping_sql = "SELECT DISTINCT feature_idx::BIGINT AS feature_idx, "
		              "'unnamed-bucket'::VARCHAR AS bucket_name FROM " +
		              chunk_quoted;
	} else {
		std::string mapping_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.mapping_table);
		mapping_sql =
		    "SELECT feature_idx::BIGINT AS feature_idx, bucket_name::VARCHAR AS bucket_name FROM " + mapping_quoted;
	}
	auto mapping_result = conn.Query(mapping_sql);
	if (mapping_result->HasError()) {
		throw InvalidInputException("Failed to read mapping table '%s': %s",
		                            bind_data.mapping_table.empty() ? bind_data.chunk_table : bind_data.mapping_table,
		                            mapping_result->GetError());
	}
	auto mapping_wrapper = make_uniq<ResultArrowArrayStreamWrapper>(std::move(mapping_result), STANDARD_VECTOR_SIZE);

	// Chunk stream: feature_idx Int64, chunk_index Int32, chunk_data. A feature's
	// chunks must be contiguous, ascending, 0-based; ORDER BY guarantees this.
	// Streamed via SendQuery (lazy fetch) so the full sequence corpus is never
	// materialized in memory at once.
	std::string chunk_sql = "SELECT feature_idx::BIGINT AS feature_idx, chunk_index::INTEGER AS chunk_index, "
	                        "chunk_data FROM " +
	                        chunk_quoted + " ORDER BY feature_idx, chunk_index";
	auto chunk_result = conn.SendQuery(chunk_sql);
	if (chunk_result->HasError()) {
		throw InvalidInputException("Failed to read chunk table '%s': %s", bind_data.chunk_table,
		                            chunk_result->GetError());
	}
	auto chunk_wrapper = make_uniq<ResultArrowArrayStreamWrapper>(std::move(chunk_result), STANDARD_VECTOR_SIZE);

	int rc =
	    rype_index_build_from_arrow(bind_data.output_path.c_str(), &chunk_wrapper->stream, &mapping_wrapper->stream,
	                                static_cast<size_t>(bind_data.k), static_cast<size_t>(bind_data.w), bind_data.salt,
	                                bind_data.orient ? 1 : 0, static_cast<size_t>(bind_data.max_memory));

	// RYpe takes ownership of both (non-NULL) streams and invokes their release
	// callbacks during the build — release the unique_ptrs to avoid a double-free.
	(void)chunk_wrapper.release();
	(void)mapping_wrapper.release();

	if (rc != 0) {
		const char *err = rype_get_last_error();
		throw IOException("rype_index_create failed to build '%s': %s", bind_data.output_path,
		                  err ? err : "unknown error");
	}

	return gstate;
}

// ============================================================================
// InitLocal
// ============================================================================
unique_ptr<LocalTableFunctionState> RypeIndexCreateTableFunction::InitLocal(ExecutionContext &context,
                                                                            TableFunctionInitInput &input,
                                                                            GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ============================================================================
// Execute — emit the single status row
// ============================================================================
void RypeIndexCreateTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	output.data[0].SetValue(0, Value(bind_data.output_path));
	output.data[1].SetValue(0, Value::INTEGER(bind_data.k));
	output.data[2].SetValue(0, Value::INTEGER(bind_data.w));
	output.data[3].SetValue(0, Value("ok"));

	output.SetCardinality(1);
	gstate.done = true;
}

// ============================================================================
// Function registration
// ============================================================================
TableFunction RypeIndexCreateTableFunction::GetFunction() {
	TableFunction tf("rype_index_create", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                 InitLocal);

	tf.named_parameters["mapping_table"] = LogicalType::VARCHAR;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["w"] = LogicalType::INTEGER;
	tf.named_parameters["salt"] = LogicalType::UBIGINT;
	tf.named_parameters["orient"] = LogicalType::BOOLEAN;
	tf.named_parameters["max_memory"] = LogicalType::BIGINT;

	return tf;
}

void RypeIndexCreateTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
