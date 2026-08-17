#include "rype_log_ratio.hpp"
#include "rype_input_stream.hpp"
#include "catalog_utils.hpp"
#include "rype_common.hpp"
#include "duckdb/common/helper.hpp"
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

	// Input stream after the output stream, connection last — see rype_classify.cpp
	// destructor for rationale.
	input_stream.reset();
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

	ParseRypeSharedParams(input, data->max_memory, data->debug);

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

	// Step 5: Build the Arrow input stream — one streaming scan of the caller's
	// relation, carrying the identifier and the sequence in the same row. See
	// rype_input_stream.hpp.
	gstate->id_map = make_uniq<RypeIdMap>(bind_data.id_type);
	RypeInputStreamOptions stream_options;
	stream_options.relation_quoted = table_quoted;
	stream_options.id_column_quoted = id_col_quoted;
	stream_options.source_name = bind_data.sequence_table;
	stream_options.id_column_name = bind_data.id_column;
	stream_options.id_type = bind_data.id_type;
	stream_options.include_pair_column = true;
	stream_options.has_sequence2 = bind_data.has_sequence2;
	// Cap the Arrow batch by bytes; batch_size below governs only how often RYpe
	// runs a classification pass. See the same block in rype_classify.cpp InitGlobal.
	stream_options.batch_bytes = GetRypeArrowBatchBytes(context);
	gstate->input_stream = BuildRypeInputStream(conn, *gstate->id_map, std::move(stream_options));

	// Step 6: Estimate batch size.
	// Log-ratio loads shards from BOTH indices per batch, so use whichever index has larger
	// shards for a conservative memory estimate. rype_calculate_batch_config accounts for shard
	// size in its memory budget, so the index with larger shards yields a smaller batch size.
	//
	// avg_read_length and is_paired come from the head of the stream's own scan
	// rather than separate queries: the caller's relation must be read exactly once
	// (#229). See RypeInputStream::SampleSizing for what that costs on is_paired.
	size_t num_shard_bytes = rype_index_largest_shard_bytes(gstate->numerator_index);
	size_t denom_shard_bytes = rype_index_largest_shard_bytes(gstate->denominator_index);
	const RypeIndex *sizing_index =
	    (denom_shard_bytes > num_shard_bytes) ? gstate->denominator_index : gstate->numerator_index;

	const auto sizing = gstate->input_stream->SampleSizing();
	const size_t memory_budget = ResolveRypeMemoryBudget(context, bind_data.max_memory);
	const size_t batch_size = ResolveRypeBatchSize(context, "rype_log_ratio", sizing_index, sizing.avg_read_length,
	                                               sizing.is_paired ? 1 : 0, memory_budget, bind_data.debug);
	gstate->input_stream->SetBatchRows(batch_size);

	// Step 7: Call RYpe log-ratio classification. The input stream stays owned by
	// this GlobalState — see rype_classify.cpp InitGlobal for why.
	int result = rype_classify_arrow_log_ratio_ex(gstate->numerator_index, gstate->denominator_index,
	                                              &gstate->input_stream->stream, bind_data.skip_threshold, batch_size,
	                                              &gstate->output_stream);

	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe log-ratio classification failed: %s", err ? err : "unknown error");
	}

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

		// Column 0: read_id (lookup from query_id, the surrogate miint assigned)
		int64_t query_id = query_ids[array_idx + query_id_array.offset];
		if (query_id < 0 || static_cast<size_t>(query_id) >= gstate.id_map->size()) {
			// size() - 1 would wrap on an empty map and report a range that contains
			// every possible value, i.e. deny the failure it is reporting.
			throw IOException("RYpe returned invalid query_id %lld (the input carried %zu rows)",
			                  static_cast<long long>(query_id), gstate.id_map->size());
		}
		gstate.id_map->Emit(output.data[0], i, static_cast<idx_t>(query_id));
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
	AddRypeSharedNamedParameters(tf);

	return tf;
}

// ============================================================================
// Register
// ============================================================================
void RypeLogRatioTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
