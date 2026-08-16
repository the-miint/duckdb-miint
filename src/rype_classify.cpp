#include "rype_classify.hpp"
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
RypeClassifyTableFunction::GlobalState::~GlobalState() {
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

	// Free RYpe negative set
	if (negative_set) {
		rype_negative_set_free(negative_set);
	}

	// Free RYpe index
	if (index) {
		rype_index_free(index);
	}

	// Drop the materialized temp table BEFORE releasing the connection that owns it.
	// This is REQUIRED, not just an early release of memory: input_connection inherits
	// the caller's TEMP catalog (#193) so the table lives in the user's session and is
	// NOT reaped when the connection is torn down. DropRypeTempTable is a no-op if
	// either is unset (e.g., bind failure) and warns rather than swallowing a failed
	// drop, since a failure is now a visible leak into the user's catalog.
	if (input_connection) {
		DropRypeTempTable(*input_connection, tmp_table_name);
	}

	// Release sub-connection LAST — RYpe's input stream (released above via
	// output_stream.release) holds a ResultArrowArrayStreamWrapper whose
	// QueryResult has a non-owning optional_ptr<ClientContext> into this connection.
	input_connection.reset();
}

// ============================================================================
// Bind
// ============================================================================
unique_ptr<FunctionData> RypeClassifyTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<Data>();

	// Required: index_path (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("rype_classify requires index_path parameter");
	}
	data->index_path = input.inputs[0].ToString();

	// Required: sequence_table (second positional parameter)
	if (input.inputs.size() < 2) {
		throw BinderException("rype_classify requires sequence_table parameter");
	}
	data->sequence_table = input.inputs[1].ToString();

	// Optional: id_column (defaults to 'read_id')
	auto id_col_param = input.named_parameters.find("id_column");
	if (id_col_param != input.named_parameters.end() && !id_col_param->second.IsNull()) {
		data->id_column = id_col_param->second.ToString();
	}

	// Optional: threshold (defaults to 0.1)
	auto threshold_param = input.named_parameters.find("threshold");
	if (threshold_param != input.named_parameters.end() && !threshold_param->second.IsNull()) {
		data->threshold = threshold_param->second.GetValue<double>();
		if (data->threshold < 0.0 || data->threshold > 1.0) {
			throw BinderException("threshold must be between 0.0 and 1.0");
		}
	}

	// Optional: negative_index
	auto neg_idx_param = input.named_parameters.find("negative_index");
	if (neg_idx_param != input.named_parameters.end() && !neg_idx_param->second.IsNull()) {
		data->negative_index_path = neg_idx_param->second.ToString();
	}

	// Note: We do NOT validate index paths here. Path validation is deferred to
	// rype_index_load() in InitGlobal() which provides proper error messages for
	// missing paths, invalid formats, and corrupted indices. Early validation here
	// would create TOCTOU race conditions and duplicate RYpe's validation logic.

	// Validate sequence table exists and has required columns. Cache whether
	// sequence2 exists, and capture the id column's storage type so read_id
	// mirrors it on output (BIGINT/UUID/VARCHAR) instead of always VARCHAR.
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
unique_ptr<GlobalTableFunctionState> RypeClassifyTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Step 1: Load RYpe index
	gstate->index = rype_index_load(bind_data.index_path.c_str());
	if (!gstate->index) {
		const char *err = rype_get_last_error();
		throw IOException("Failed to load RYpe index '%s': %s", bind_data.index_path, err ? err : "unknown error");
	}

	// Step 2: Load negative set if specified
	if (!bind_data.negative_index_path.empty()) {
		RypeIndex *neg_index = rype_index_load(bind_data.negative_index_path.c_str());
		if (!neg_index) {
			const char *err = rype_get_last_error();
			throw IOException("Failed to load negative index '%s': %s", bind_data.negative_index_path,
			                  err ? err : "unknown error");
		}
		gstate->negative_set = rype_negative_set_create(neg_index);
		rype_index_free(neg_index);
		if (!gstate->negative_set) {
			const char *err = rype_get_last_error();
			throw IOException("Failed to create negative set: %s", err ? err : "unknown error");
		}
	}

	// Step 3: Build read_id mapping and query sequence data for RYpe
	// Store connection in GlobalState so it outlives the RYpe input stream.
	// RYpe lazily consumes the input Arrow stream (backed by ResultArrowArrayStreamWrapper
	// whose QueryResult holds a non-owning optional_ptr<ClientContext>). If the connection
	// is destroyed before RYpe finishes consuming, the dangling pointer causes use-after-free.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->input_connection = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->input_connection);
	auto &conn = *gstate->input_connection;

	// Export BLOB with 64-bit offsets. See ConfigureRypeArrowExport in rype_common.hpp
	// for why this must not be BinaryView (#222).
	ConfigureRypeArrowExport(conn);

	std::string id_col_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.id_column);
	std::string table_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.sequence_table);

	// Materialize the user's sequence_table into a per-call TEMP table with an
	// explicit id column, then read read_ids and the sequence stream from it.
	// See rype_common.hpp for the design rationale (row-index ↔ read_id fix).
	size_t avg_read_length = 0;
	gstate->tmp_table_name =
	    MaterializeRypeInputTempTable(conn, table_quoted, id_col_quoted, bind_data.sequence_table,
	                                  bind_data.has_sequence2, "_rype_classify_", gstate->read_ids, avg_read_length);

	// Estimate batch size. RYpe processes one batch at a time; using
	// STANDARD_VECTOR_SIZE (2048) causes shard I/O to dominate.
	//
	// is_paired follows sequence2 CONTENT, not the column's presence (#199) — see
	// TableHasPairedContent in rype_common.hpp for why, and why it probes the temp table.
	int is_paired = TableHasPairedContent(conn, KeywordHelper::WriteOptionallyQuoted(gstate->tmp_table_name)) ? 1 : 0;
	// is_large_binary=1 tells RYpe to skip its 2 GiB batch cap. That is only sound
	// because ConfigureRypeArrowExport pinned this connection to Arrow LargeBinary
	// (i64 offsets) above — the two must be changed together (#222).
	size_t batch_size = rype_recommend_batch_size(gstate->index, avg_read_length, is_paired, 0, 1);
	if (batch_size == 0) {
		const char *err = rype_get_last_error();
		Printer::Print(StringUtil::Format("Warning: rype_recommend_batch_size failed (%s), using default",
		                                  err ? err : "unknown error"));
		batch_size = STANDARD_VECTOR_SIZE;
	}

	// Build the Arrow input stream from the same temp table, ordered by id.
	auto input_wrapper = BuildRypeArrowInput(conn, gstate->tmp_table_name, /*include_pair_column=*/true, batch_size);
	ArrowArrayStream *input_stream = &input_wrapper->stream;

	// Step 5: Call RYpe classify
	// IMPORTANT: RYpe takes ownership of the stream via ArrowArrayStreamReader::from_raw().
	// When RYpe is done with the stream, it calls the release callback, which deletes
	// the ResultArrowArrayStreamWrapper. We must release() from unique_ptr to avoid double-free,
	// but ONLY after confirming success - if RYpe fails early it may not take ownership.
	int result = rype_classify_arrow(gstate->index, gstate->negative_set, input_stream, bind_data.threshold,
	                                 &gstate->output_stream);

	if (result != 0) {
		// RYpe failed - it may not have taken ownership of the stream, so let
		// unique_ptr destructor clean up the wrapper normally
		const char *err = rype_get_last_error();
		throw IOException("RYpe classification failed: %s", err ? err : "unknown error");
	}

	// Success - RYpe now owns the stream. Transfer ownership out of unique_ptr
	// to prevent double-free when the wrapper is deleted by MyStreamRelease.
	(void)input_wrapper.release();

	// Step 6: Get output schema
	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	// Step 7: Parse Arrow schema for conversion
	ArrowTableFunction::PopulateArrowTableSchema(context, gstate->arrow_table, gstate->output_schema);

	// Verify RYpe's output schema matches expected columns (query_id, bucket_id, score)
	if (gstate->arrow_table.GetColumns().size() != 3) {
		throw IOException("RYpe classify returned %zu columns, expected 3", gstate->arrow_table.GetColumns().size());
	}

	return gstate;
}

// ============================================================================
// InitLocal
// ============================================================================
unique_ptr<LocalTableFunctionState> RypeClassifyTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>(context.client);
}

// ============================================================================
// Execute
// ============================================================================
void RypeClassifyTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();
	auto &bind_data = data_p.bind_data->Cast<Data>();

	// No mutex needed - MaxThreads() returns 1, enforcing single-threaded execution

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

	// RYpe output schema: query_id (Int64), bucket_id (UInt32), score (Float64)
	// Our output schema: read_id (mirrors the id_column's type -- VARCHAR, BIGINT or UUID;
	// see Bind, which sets types[0] = id_type), bucket_id (UINTEGER), bucket_name (VARCHAR),
	// score (DOUBLE)

	// --- Column 0 (read_id) and Column 2 (bucket_name): manual transformation ---
	// These require per-row lookups so they cannot be zero-copy.
	auto &query_id_array = *batch.children[0];
	auto &bucket_id_array = *batch.children[1];

	if (!query_id_array.buffers[1] || !bucket_id_array.buffers[1]) {
		throw IOException("Arrow array has null data buffer in RYpe classify output");
	}

	auto query_ids = reinterpret_cast<const int64_t *>(query_id_array.buffers[1]);
	auto bucket_ids = reinterpret_cast<const uint32_t *>(bucket_id_array.buffers[1]);

	for (idx_t i = 0; i < to_output; i++) {
		idx_t array_idx = gstate.batch_offset + i + batch.offset;

		// Column 0: read_id (lookup from query_id which is a row index)
		int64_t query_id = query_ids[array_idx + query_id_array.offset];
		if (query_id < 0 || static_cast<size_t>(query_id) >= gstate.read_ids.size()) {
			throw IOException("RYpe returned invalid query_id %lld (expected 0-%zu)", static_cast<long long>(query_id),
			                  gstate.read_ids.size() - 1);
		}
		EmitIdCell(output.data[0], i, gstate.read_ids[query_id], bind_data.id_type);

		// Column 2: bucket_name (lookup from index using bucket_id)
		uint32_t bucket_id = bucket_ids[array_idx + bucket_id_array.offset];
		const char *name = rype_bucket_name(gstate.index, bucket_id);
		if (!name) {
			throw IOException("RYpe returned unknown bucket_id %u - index may be corrupted", bucket_id);
		}
		FlatVector::GetData<string_t>(output.data[2])[i] = StringVector::AddString(output.data[2], name);
	}

	// --- Column 1 (bucket_id) and Column 3 (score): zero-copy via Arrow conversion ---
	// Arrow col 1 (UInt32) → DuckDB col 1, Arrow col 2 (Float64) → DuckDB col 3
	auto &arrow_columns = gstate.arrow_table.GetColumns();

	// bucket_id: Arrow col 1 → DuckDB col 1
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

	// score: Arrow col 2 → DuckDB col 3
	{
		auto &array = *batch.children[2];
		auto &arrow_type = *arrow_columns.at(2);
		auto &array_state = lstate.GetState(2);
		array_state.owned_data = gstate.current_chunk;
		ArrowToDuckDBConversion::SetValidityMask(output.data[3], array, gstate.batch_offset, to_output, batch.offset,
		                                         -1);
		ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[3], array, gstate.batch_offset, array_state, to_output,
		                                             arrow_type);
	}

	gstate.batch_offset += to_output;
}

// ============================================================================
// GetFunction
// ============================================================================
TableFunction RypeClassifyTableFunction::GetFunction() {
	TableFunction tf("rype_classify", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                 InitLocal);

	// Named parameters
	tf.named_parameters["id_column"] = LogicalType::VARCHAR;
	tf.named_parameters["threshold"] = LogicalType::DOUBLE;
	tf.named_parameters["negative_index"] = LogicalType::VARCHAR;

	return tf;
}

// ============================================================================
// Register
// ============================================================================
void RypeClassifyTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
