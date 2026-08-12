#include "rype_log_ratio.hpp"
#include "catalog_utils.hpp"
#include "rype_common.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/main/config.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

// ============================================================================
// GlobalState destructor - cleanup RYpe resources
// ============================================================================
RypeLogRatioTableFunction::GlobalState::~GlobalState() {
	// Release shared_ptr to current batch. Any DuckDB Vectors that still
	// reference this batch via ArrowAuxiliaryData will keep it alive.
	current_chunk.reset();

	// IMPORTANT: Clear arrow_table BEFORE releasing output_schema.
	// ArrowTableSchema holds pointers into the schema data, so we must
	// clear it before freeing that memory.
	arrow_table = ArrowTableSchema();

	// Now safe to release output schema
	if (output_schema.release) {
		output_schema.release(&output_schema);
	}

	// Release output stream
	if (output_stream.release) {
		output_stream.release(&output_stream);
	}

	// Free RYpe indices (both)
	if (denominator_index) {
		rype_index_free(denominator_index);
	}
	if (numerator_index) {
		rype_index_free(numerator_index);
	}

	// Drop the materialized temp table BEFORE releasing input_connection — see
	// rype_classify.cpp destructor for rationale.
	if (input_connection) {
		DropRypeTempTable(*input_connection, tmp_table_name);
	}

	// Release sub-connection LAST — see rype_classify.cpp destructor for rationale.
	input_connection.reset();
}

// ============================================================================
// Bind
// ============================================================================
unique_ptr<FunctionData> RypeLogRatioTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<Data>();

	// Required: numerator_path (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("rype_log_ratio requires numerator index path parameter");
	}
	data->numerator_path = input.inputs[0].ToString();

	// Required: denominator_path (second positional parameter)
	if (input.inputs.size() < 2) {
		throw BinderException("rype_log_ratio requires denominator index path parameter");
	}
	data->denominator_path = input.inputs[1].ToString();

	// Required: sequence_table (third positional parameter)
	if (input.inputs.size() < 3) {
		throw BinderException("rype_log_ratio requires sequence_table parameter");
	}
	data->sequence_table = input.inputs[2].ToString();

	// Optional: id_column (defaults to 'read_id')
	auto id_col_param = input.named_parameters.find("id_column");
	if (id_col_param != input.named_parameters.end() && !id_col_param->second.IsNull()) {
		data->id_column = id_col_param->second.ToString();
	}

	// Optional: skip_threshold (defaults to 0.5)
	// Unlike rype_classify's threshold (which must be [0,1]), negative values are intentionally
	// allowed here: per the RYpe API, skip_threshold <= 0.0 disables the fast-path entirely,
	// forcing exact classification against both indices for every read.
	auto skip_param = input.named_parameters.find("skip_threshold");
	if (skip_param != input.named_parameters.end() && !skip_param->second.IsNull()) {
		data->skip_threshold = skip_param->second.GetValue<double>();
		if (std::isnan(data->skip_threshold)) {
			throw BinderException("skip_threshold must not be NaN");
		}
		if (data->skip_threshold > 1.0) {
			throw BinderException("skip_threshold must be <= 1.0");
		}
	}

	// Validate sequence table exists and has required columns. Capture the id
	// column's storage type so read_id mirrors it on output (BIGINT/UUID/VARCHAR)
	// instead of always VARCHAR.
	auto table_info = ValidateSequenceTable(context, data->sequence_table, data->id_column);
	data->has_sequence2 = table_info.has_sequence2;
	data->id_type = table_info.id_type;
	data->types[0] = data->id_type;

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
// InitGlobal
// ============================================================================
unique_ptr<GlobalTableFunctionState> RypeLogRatioTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Step 1: Load numerator index
	gstate->numerator_index = rype_index_load(bind_data.numerator_path.c_str());
	if (!gstate->numerator_index) {
		const char *err = rype_get_last_error();
		throw IOException("Failed to load numerator index '%s': %s", bind_data.numerator_path,
		                  err ? err : "unknown error");
	}

	// Step 2: Load denominator index
	gstate->denominator_index = rype_index_load(bind_data.denominator_path.c_str());
	if (!gstate->denominator_index) {
		const char *err = rype_get_last_error();
		throw IOException("Failed to load denominator index '%s': %s", bind_data.denominator_path,
		                  err ? err : "unknown error");
	}

	// Step 3: Validate indices are compatible for log-ratio
	int validate_result = rype_validate_log_ratio_indices(gstate->numerator_index, gstate->denominator_index);
	if (validate_result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("Log-ratio index validation failed: %s", err ? err : "unknown error");
	}

	// Step 4: Build read_id mapping and query sequence data for RYpe
	// Store connection in GlobalState — see rype_classify.cpp InitGlobal for rationale.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->input_connection = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->input_connection);
	auto &conn = *gstate->input_connection;

	// Export BLOB with 64-bit offsets — see ConfigureRypeArrowExport in rype_common.hpp (#222).
	ConfigureRypeArrowExport(conn);

	std::string id_col_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.id_column);
	std::string table_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.sequence_table);

	// Materialize the user's sequence_table into a per-call TEMP table with an
	// explicit id column, then read read_ids and the sequence stream from it.
	// See rype_common.hpp for the design rationale (row-index ↔ read_id fix).
	size_t avg_read_length = 0;
	gstate->tmp_table_name =
	    MaterializeRypeInputTempTable(conn, table_quoted, id_col_quoted, bind_data.sequence_table,
	                                  bind_data.has_sequence2, "_rype_log_ratio_", gstate->read_ids, avg_read_length);

	// Step 5: Estimate batch size.
	// Log-ratio loads shards from BOTH indices per batch, so use whichever index has larger
	// shards for a conservative memory estimate. rype_recommend_batch_size accounts for shard
	// size in its memory budget, so the index with larger shards yields a smaller batch size.
	//
	// is_paired follows sequence2 CONTENT, not the column's presence (#199) — see
	// rype_classify.cpp InitGlobal and TableHasPairedContent in rype_common.hpp for rationale.
	int is_paired = TableHasPairedContent(conn, KeywordHelper::WriteOptionallyQuoted(gstate->tmp_table_name)) ? 1 : 0;

	size_t num_shard_bytes = rype_index_largest_shard_bytes(gstate->numerator_index);
	size_t denom_shard_bytes = rype_index_largest_shard_bytes(gstate->denominator_index);
	const RypeIndex *sizing_index =
	    (denom_shard_bytes > num_shard_bytes) ? gstate->denominator_index : gstate->numerator_index;

	// is_large_binary=1 tells RYpe to skip its 2 GiB batch cap. That is only sound
	// because ConfigureRypeArrowExport pinned this connection to Arrow LargeBinary
	// (i64 offsets) above — the two must be changed together (#222).
	size_t batch_size = rype_recommend_batch_size(sizing_index, avg_read_length, is_paired, 0, 1);
	if (batch_size == 0) {
		// rype_recommend_batch_size returns 0 on error — log but use safe fallback
		const char *err = rype_get_last_error();
		Printer::Print(StringUtil::Format("Warning: rype_recommend_batch_size failed (%s), using default",
		                                  err ? err : "unknown error"));
		batch_size = STANDARD_VECTOR_SIZE;
	}

	// Step 6: Build the Arrow input stream from the same temp table, ordered by id.
	auto input_wrapper = BuildRypeArrowInput(conn, gstate->tmp_table_name, /*include_pair_column=*/true, batch_size);
	ArrowArrayStream *input_stream = &input_wrapper->stream;

	// Step 7: Call RYpe log-ratio classification
	int result = rype_classify_arrow_log_ratio(gstate->numerator_index, gstate->denominator_index, input_stream,
	                                           bind_data.skip_threshold, &gstate->output_stream);

	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe log-ratio classification failed: %s", err ? err : "unknown error");
	}

	// Success - RYpe now owns the stream
	(void)input_wrapper.release();

	// Step 8: Get output schema
	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	// Step 9: Parse Arrow schema for conversion
	ArrowTableFunction::PopulateArrowTableSchema(context, gstate->arrow_table, gstate->output_schema);

	// Verify RYpe's output schema matches expected columns (query_id, log_ratio, fast_path)
	if (gstate->arrow_table.GetColumns().size() != 3) {
		throw IOException("RYpe log-ratio returned %zu columns, expected 3", gstate->arrow_table.GetColumns().size());
	}

	return gstate;
}

// ============================================================================
// InitLocal
// ============================================================================
unique_ptr<LocalTableFunctionState> RypeLogRatioTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>(context.client);
}

// ============================================================================
// Execute
// ============================================================================
void RypeLogRatioTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();
	auto &bind_data = data_p.bind_data->Cast<Data>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	// Fetch next batch if needed
	auto batch_length = gstate.current_chunk ? gstate.current_chunk->arrow_array.length : 0;
	while (gstate.batch_offset >= static_cast<idx_t>(batch_length) || !gstate.current_chunk) {
		gstate.current_chunk.reset();
		gstate.batch_offset = 0;
		lstate.ResetStates();

		// Get next batch from RYpe
		auto wrapper = make_shared_ptr<ArrowArrayWrapper>();
		if (gstate.output_stream.get_next(&gstate.output_stream, &wrapper->arrow_array) != 0) {
			const char *err = gstate.output_stream.get_last_error(&gstate.output_stream);
			throw IOException("Error getting next batch from RYpe: %s", err ? err : "unknown error");
		}

		// Check if stream is exhausted
		if (!wrapper->arrow_array.release) {
			gstate.done = true;
			output.SetCardinality(0);
			return;
		}

		gstate.current_chunk = std::move(wrapper);
		batch_length = gstate.current_chunk->arrow_array.length;
	}

	auto &batch = gstate.current_chunk->arrow_array;
	idx_t remaining = static_cast<idx_t>(batch.length) - gstate.batch_offset;
	idx_t to_output = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	output.SetCardinality(to_output);

	// RYpe output schema: query_id (Int64), log_ratio (Float64), fast_path (Int32)
	// Our output schema: read_id (VARCHAR), log_ratio (DOUBLE), fast_path (INTEGER)

	// --- Column 0 (read_id): manual transformation ---
	auto &query_id_array = *batch.children[0];

	if (!query_id_array.buffers[1]) {
		throw IOException("Arrow array has null data buffer in RYpe log-ratio output");
	}

	auto query_ids = reinterpret_cast<const int64_t *>(query_id_array.buffers[1]);

	for (idx_t i = 0; i < to_output; i++) {
		idx_t array_idx = gstate.batch_offset + i + batch.offset;

		// Column 0: read_id (lookup from query_id which is a row index)
		int64_t query_id = query_ids[array_idx + query_id_array.offset];
		if (query_id < 0 || static_cast<size_t>(query_id) >= gstate.read_ids.size()) {
			throw IOException("RYpe returned invalid query_id %lld (expected 0-%zu)", static_cast<long long>(query_id),
			                  gstate.read_ids.size() - 1);
		}
		EmitIdCell(output.data[0], i, gstate.read_ids[query_id], bind_data.id_type);
	}

	// --- Column 1 (log_ratio) and Column 2 (fast_path): zero-copy via Arrow conversion ---
	auto &arrow_columns = gstate.arrow_table.GetColumns();

	// log_ratio: Arrow col 1 → DuckDB col 1
	{
		auto &array = *batch.children[1];
		auto &arrow_type = *arrow_columns.at(1);
		auto &array_state = lstate.GetState(1);
		array_state.owned_data = gstate.current_chunk;
		ArrowToDuckDBConversion::SetValidityMask(output.data[1], array, gstate.batch_offset, to_output, batch.offset,
		                                         -1);
		ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[1], array, gstate.batch_offset, array_state, to_output,
		                                             arrow_type);
	}

	// fast_path: Arrow col 2 → DuckDB col 2
	{
		auto &array = *batch.children[2];
		auto &arrow_type = *arrow_columns.at(2);
		auto &array_state = lstate.GetState(2);
		array_state.owned_data = gstate.current_chunk;
		ArrowToDuckDBConversion::SetValidityMask(output.data[2], array, gstate.batch_offset, to_output, batch.offset,
		                                         -1);
		ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[2], array, gstate.batch_offset, array_state, to_output,
		                                             arrow_type);
	}

	gstate.batch_offset += to_output;
}

// ============================================================================
// GetFunction
// ============================================================================
TableFunction RypeLogRatioTableFunction::GetFunction() {
	TableFunction tf("rype_log_ratio", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute,
	                 Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["id_column"] = LogicalType::VARCHAR;
	tf.named_parameters["skip_threshold"] = LogicalType::DOUBLE;

	return tf;
}

// ============================================================================
// Register
// ============================================================================
void RypeLogRatioTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
