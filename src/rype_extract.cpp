#include "rype_extract.hpp"
#include "rype_input_stream.hpp"
#include "catalog_utils.hpp"
#include "rype_common.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/main/config.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

// ============================================================================
// GlobalState destructor
// ============================================================================
RypeExtractGlobalState::~RypeExtractGlobalState() {
	// Release shared_ptr to current batch. Any DuckDB Vectors that still
	// reference this batch via ArrowAuxiliaryData will keep it alive.
	current_chunk.reset();

	// Clear arrow_table BEFORE releasing output_schema — it holds pointers into schema data.
	arrow_table = ArrowTableSchema();

	if (output_schema.release) {
		output_schema.release(&output_schema);
	}
	if (output_stream.release) {
		output_stream.release(&output_stream);
	}

	// Input stream after the output stream, connection last — see rype_classify.cpp
	// destructor for rationale.
	input_stream.reset();
	input_connection.reset();
}

// ============================================================================
// Shared helpers
// ============================================================================

// Shared Bind logic for both extraction functions.
static unique_ptr<RypeExtractData> BindExtraction(ClientContext &context, TableFunctionBindInput &input,
                                                  const std::string &function_name) {
	auto data = make_uniq<RypeExtractData>();

	if (input.inputs.size() < 3) {
		throw BinderException("%s requires sequence_table, k, and w parameters", function_name);
	}

	data->sequence_table = input.inputs[0].ToString();
	data->k = static_cast<size_t>(input.inputs[1].GetValue<int64_t>());
	data->w = static_cast<size_t>(input.inputs[2].GetValue<int64_t>());

	// RY-space uses 1 bit per base (purine/pyrimidine), so k-mer fits in u64 only
	// for k ∈ {16, 32, 64}. These match hardware word sizes for efficient hashing.
	if (data->k != 16 && data->k != 32 && data->k != 64) {
		throw BinderException("k must be 16, 32, or 64 (got %zu)", data->k);
	}
	if (data->w == 0) {
		throw BinderException("w must be > 0");
	}

	auto salt_param = input.named_parameters.find("salt");
	if (salt_param != input.named_parameters.end() && !salt_param->second.IsNull()) {
		data->salt = salt_param->second.GetValue<uint64_t>();
	}

	auto id_col_param = input.named_parameters.find("id_column");
	if (id_col_param != input.named_parameters.end() && !id_col_param->second.IsNull()) {
		data->id_column = id_col_param->second.ToString();
	}

	ParseRypeSharedParams(input, data->max_memory, data->debug);

	// Capture the id column's storage type so read_id mirrors it on output
	// (BIGINT/UUID/VARCHAR) instead of always VARCHAR. Extraction is single-end,
	// so has_sequence2 is irrelevant here.
	data->id_type = ValidateSequenceTable(context, data->sequence_table, data->id_column).id_type;

	return data;
}

// Build the id map and Arrow input stream for extraction. The stream is owned by
// the returned GlobalState and stays owned by it — see the ownership note in
// rype_input_stream.hpp.
static unique_ptr<RypeExtractGlobalState>
BuildExtractionInputStream(ClientContext &context, const RypeExtractData &bind_data, const char *label) {
	auto gstate = make_uniq<RypeExtractGlobalState>();

	// Store connection in GlobalState — see rype_classify.cpp InitGlobal for rationale.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->input_connection = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->input_connection);
	auto &conn = *gstate->input_connection;

	// Export BLOB with 64-bit offsets — see ConfigureRypeArrowExport in rype_common.hpp (#222).
	ConfigureRypeArrowExport(conn);

	std::string id_col_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.id_column);
	std::string table_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.sequence_table);

	// One streaming scan of the caller's relation, carrying the identifier and the
	// sequence in the same row. Extraction is single-sequence, so no pair column.
	// See rype_input_stream.hpp.
	gstate->id_map = make_uniq<RypeIdMap>(bind_data.id_type);
	RypeInputStreamOptions stream_options;
	stream_options.relation_quoted = table_quoted;
	stream_options.id_column_quoted = id_col_quoted;
	stream_options.source_name = bind_data.sequence_table;
	stream_options.id_column_name = bind_data.id_column;
	stream_options.id_type = bind_data.id_type;
	stream_options.include_pair_column = false;
	// Extraction runs per batch with no index pass behind it, so a smaller Arrow
	// batch costs nothing here — it just bounds the sequence bytes resident at
	// once, which the row-derived batch_size below cannot.
	stream_options.batch_bytes = GetRypeArrowBatchBytes(context);
	gstate->input_stream = BuildRypeInputStream(conn, *gstate->id_map, std::move(stream_options));

	// Estimate batch size — extraction has no index overhead, just sequence data
	// + minimizer lists. avg_read_length comes from the head of the stream's own
	// scan rather than a separate query: the caller's relation must be read
	// exactly once (#229). Extraction ignores sequence2 entirely, so is_paired
	// does not enter the estimate.
	const size_t avg_read_length = gstate->input_stream->SampleSizing().avg_read_length;
	size_t minimizers_per_read = (bind_data.w > 0 && avg_read_length > bind_data.k)
	                                 ? ((avg_read_length - bind_data.k + 1) / bind_data.w + 1)
	                                 : 1;
	size_t record_cost = avg_read_length + minimizers_per_read * sizeof(uint64_t);

	// Budget comes from the same resolver the classify paths use, so `max_memory`
	// means the same thing on every RYpe function and the DuckDB-aware default
	// applies here too (#204). A resolver result of 0 means detection failed;
	// fall back to the row floor below rather than dividing by it.
	const size_t budget = ResolveRypeMemoryBudget(context, bind_data.max_memory);
	size_t batch_size = budget / record_cost;
	if (batch_size < 1000) {
		batch_size = STANDARD_VECTOR_SIZE;
	}
	gstate->input_stream->SetBatchRows(batch_size);

	// Extraction sizes its own batch rather than going through
	// ResolveRypeBatchSize: there is no index, so RYpe's shard-aware estimate
	// does not apply and there is no RypeIndex to pass it. Report the same shape
	// of numbers anyway so `debug` means one thing across the RYpe functions (#204).
	if (bind_data.debug) {
		miint::EmitWarning(context,
		                   "%s debug: extraction batch = %llu reads; memory budget %.2f GB, "
		                   "estimated %llu bytes per read (%llu sequence + %llu minimizers); "
		                   "avg read length %llu",
		                   label, (unsigned long long)batch_size, double(budget) / 1e9, (unsigned long long)record_cost,
		                   (unsigned long long)avg_read_length,
		                   (unsigned long long)(minimizers_per_read * sizeof(uint64_t)),
		                   (unsigned long long)avg_read_length);
	}

	return gstate;
}

// Shared Execute logic: fetches batches from the Arrow output stream and
// uses DuckDB's built-in Arrow conversion for zero-copy on list columns.
// Column 0 (surrogate id → read_id) is resolved through gstate.id_map, which
// already knows the caller's id type.
// num_list_cols is the number of list columns starting at index 1.
static void ExecuteExtraction(RypeExtractGlobalState &gstate, RypeExtractLocalState &lstate, DataChunk &output,
                              idx_t num_list_cols) {
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

		auto wrapper = make_shared_ptr<ArrowArrayWrapper>();
		if (gstate.output_stream.get_next(&gstate.output_stream, &wrapper->arrow_array) != 0) {
			const char *err = gstate.output_stream.get_last_error(&gstate.output_stream);
			throw IOException("Error getting next batch from RYpe: %s", err ? err : "unknown error");
		}

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

	// Column 0: id (Int64) → read_id (mirrors id_column type) — manual transformation.
	// Offset calculation follows rype_classify pattern: parent batch offset + child array offset.
	auto &id_array = *batch.children[0];
	if (!id_array.buffers[1]) {
		throw IOException("Arrow array for query_id has null data buffer");
	}
	auto query_ids = reinterpret_cast<const int64_t *>(id_array.buffers[1]);

	for (idx_t i = 0; i < to_output; i++) {
		idx_t array_idx = gstate.batch_offset + i + batch.offset;
		int64_t query_id = query_ids[array_idx + id_array.offset];

		if (query_id < 0 || static_cast<size_t>(query_id) >= gstate.id_map->size()) {
			// size() - 1 would wrap on an empty map and report a range that contains
			// every possible value, i.e. deny the failure it is reporting.
			throw IOException("RYpe returned invalid query_id %lld (the input carried %zu rows)",
			                  static_cast<long long>(query_id), gstate.id_map->size());
		}
		gstate.id_map->Emit(output.data[0], i, static_cast<idx_t>(query_id));
	}

	// Columns 1..N: List<UInt64> — use DuckDB's built-in Arrow-to-DuckDB conversion.
	// This achieves zero-copy for the uint64 child data via DirectConversion/FlatVector::SetData.
	auto &arrow_columns = gstate.arrow_table.GetColumns();
	for (idx_t col = 1; col <= num_list_cols; col++) {
		auto &array = *batch.children[col];
		auto &arrow_type = *arrow_columns.at(col);
		auto &array_state = lstate.GetState(col);

		// Keep batch alive for zero-copy — Vectors will reference this via ArrowAuxiliaryData
		array_state.owned_data = gstate.current_chunk;

		ArrowToDuckDBConversion::SetValidityMask(output.data[col], array, gstate.batch_offset, to_output, batch.offset,
		                                         -1);
		ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[col], array, gstate.batch_offset, array_state,
		                                             to_output, arrow_type);
	}

	gstate.batch_offset += to_output;
}

// ============================================================================
// rype_extract_minimizer_set
// ============================================================================
unique_ptr<FunctionData> RypeExtractMinimizerSetTableFunction::Bind(ClientContext &context,
                                                                    TableFunctionBindInput &input,
                                                                    vector<LogicalType> &return_types,
                                                                    vector<string> &names) {
	auto data = BindExtraction(context, input, "rype_extract_minimizer_set");

	data->names = {"read_id", "fwd_set", "rc_set"};
	data->types = {data->id_type, LogicalType::LIST(LogicalType::UBIGINT), LogicalType::LIST(LogicalType::UBIGINT)};

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> RypeExtractMinimizerSetTableFunction::InitGlobal(ClientContext &context,
                                                                                      TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<RypeExtractData>();

	auto gstate = BuildExtractionInputStream(context, bind_data, "rype_extract_minimizer_set");

	int result = rype_extract_minimizer_set_arrow(&gstate->input_stream->stream, bind_data.k, bind_data.w,
	                                              bind_data.salt, &gstate->output_stream);
	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe minimizer set extraction failed: %s", err ? err : "unknown error");
	}

	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	ArrowTableFunction::PopulateArrowTableSchema(context, gstate->arrow_table, gstate->output_schema);

	// Verify RYpe's output schema matches our declared columns (id + 2 list columns)
	if (gstate->arrow_table.GetColumns().size() != bind_data.names.size()) {
		throw IOException("RYpe returned %zu columns, expected %zu", gstate->arrow_table.GetColumns().size(),
		                  bind_data.names.size());
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState>
RypeExtractMinimizerSetTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                GlobalTableFunctionState *global_state) {
	return make_uniq<RypeExtractLocalState>(context.client);
}

void RypeExtractMinimizerSetTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                   DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<RypeExtractGlobalState>();
	auto &lstate = data_p.local_state->Cast<RypeExtractLocalState>();
	// Output: read_id, fwd_set, rc_set → 2 list columns
	ExecuteExtraction(gstate, lstate, output, 2);
}

TableFunction RypeExtractMinimizerSetTableFunction::GetFunction() {
	TableFunction tf("rype_extract_minimizer_set", {LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::BIGINT},
	                 Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["salt"] = LogicalType::UBIGINT;
	tf.named_parameters["id_column"] = LogicalType::VARCHAR;
	AddRypeSharedNamedParameters(tf);
	return tf;
}

void RypeExtractMinimizerSetTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

// ============================================================================
// rype_extract_strand_minimizers
// ============================================================================
unique_ptr<FunctionData> RypeExtractStrandMinimizersTableFunction::Bind(ClientContext &context,
                                                                        TableFunctionBindInput &input,
                                                                        vector<LogicalType> &return_types,
                                                                        vector<string> &names) {
	auto data = BindExtraction(context, input, "rype_extract_strand_minimizers");

	data->names = {"read_id", "fwd_hashes", "fwd_positions", "rc_hashes", "rc_positions"};
	data->types = {data->id_type, LogicalType::LIST(LogicalType::UBIGINT), LogicalType::LIST(LogicalType::UBIGINT),
	               LogicalType::LIST(LogicalType::UBIGINT), LogicalType::LIST(LogicalType::UBIGINT)};

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState>
RypeExtractStrandMinimizersTableFunction::InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind_data = input.bind_data->Cast<RypeExtractData>();

	auto gstate = BuildExtractionInputStream(context, bind_data, "rype_extract_strand_minimizers");

	int result = rype_extract_strand_minimizers_arrow(&gstate->input_stream->stream, bind_data.k, bind_data.w,
	                                                  bind_data.salt, &gstate->output_stream);
	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe strand minimizer extraction failed: %s", err ? err : "unknown error");
	}

	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	ArrowTableFunction::PopulateArrowTableSchema(context, gstate->arrow_table, gstate->output_schema);

	// Verify RYpe's output schema matches our declared columns (id + 4 list columns)
	if (gstate->arrow_table.GetColumns().size() != bind_data.names.size()) {
		throw IOException("RYpe returned %zu columns, expected %zu", gstate->arrow_table.GetColumns().size(),
		                  bind_data.names.size());
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState>
RypeExtractStrandMinimizersTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<RypeExtractLocalState>(context.client);
}

void RypeExtractStrandMinimizersTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                       DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<RypeExtractGlobalState>();
	auto &lstate = data_p.local_state->Cast<RypeExtractLocalState>();
	// Output: read_id, fwd_hashes, fwd_positions, rc_hashes, rc_positions → 4 list columns
	ExecuteExtraction(gstate, lstate, output, 4);
}

TableFunction RypeExtractStrandMinimizersTableFunction::GetFunction() {
	TableFunction tf("rype_extract_strand_minimizers", {LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::BIGINT},
	                 Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["salt"] = LogicalType::UBIGINT;
	tf.named_parameters["id_column"] = LogicalType::VARCHAR;
	AddRypeSharedNamedParameters(tf);
	return tf;
}

void RypeExtractStrandMinimizersTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
