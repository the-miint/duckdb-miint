#include "rype_extract.hpp"
#include "rype_common.hpp"
#include "duckdb/common/arrow/result_arrow_wrapper.hpp"
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
}

// ============================================================================
// LocalState helpers
// ============================================================================
ArrowArrayScanState &RypeExtractLocalState::GetState(idx_t col_idx) {
	auto it = array_states.find(col_idx);
	if (it == array_states.end()) {
		auto state = make_uniq<ArrowArrayScanState>(context);
		auto &ref = *state;
		array_states.emplace(col_idx, std::move(state));
		return ref;
	}
	return *it->second;
}

void RypeExtractLocalState::ResetStates() {
	for (auto &state : array_states) {
		state.second->Reset();
	}
}

RypeExtractLocalState::~RypeExtractLocalState() {
	array_states.clear();
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

	ValidateSequenceTable(context, data->sequence_table, data->id_column);

	return data;
}

// Build read_ids vector and create Arrow input stream for extraction.
static unique_ptr<RypeExtractGlobalState>
BuildExtractionInputStream(ClientContext &context, const RypeExtractData &bind_data,
                           ArrowArrayStream **out_input_stream,
                           unique_ptr<ResultArrowArrayStreamWrapper> &out_wrapper) {
	auto gstate = make_uniq<RypeExtractGlobalState>();

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	std::string id_col_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.id_column);
	std::string table_quoted = KeywordHelper::WriteOptionallyQuoted(bind_data.sequence_table);

	// Collect read_ids
	std::string id_query = "SELECT " + id_col_quoted + " FROM " + table_quoted;
	auto id_result = conn.Query(id_query);
	if (id_result->HasError()) {
		throw InvalidInputException("Failed to read from sequence table '%s': %s", bind_data.sequence_table,
		                            id_result->GetError());
	}

	gstate->read_ids.reserve(id_result->RowCount());
	auto &id_materialized = id_result->Cast<MaterializedQueryResult>();
	while (auto chunk = id_materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			gstate->read_ids.push_back(chunk->data[0].GetValue(i).ToString());
		}
	}

	// Query sequence data — extraction only uses single sequence (no pair_sequence).
	// RYpe extraction expects: id (Int64), sequence (Binary)
	std::string query =
	    "SELECT (row_number() OVER () - 1)::BIGINT as id, sequence1::BLOB as sequence FROM " + table_quoted;

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from sequence table '%s': %s", bind_data.sequence_table,
		                            query_result->GetError());
	}

	out_wrapper = make_uniq<ResultArrowArrayStreamWrapper>(std::move(query_result), STANDARD_VECTOR_SIZE);
	*out_input_stream = &out_wrapper->stream;

	return gstate;
}

// Shared Execute logic: fetches batches from the Arrow output stream and
// uses DuckDB's built-in Arrow conversion for zero-copy on list columns.
// Column 0 (id → read_id) requires manual transformation.
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

	// Column 0: id (Int64) → read_id (VARCHAR) — manual transformation.
	// Offset calculation follows rype_classify pattern: parent batch offset + child array offset.
	auto &id_array = *batch.children[0];
	if (!id_array.buffers[1]) {
		throw IOException("Arrow array for query_id has null data buffer");
	}
	auto query_ids = reinterpret_cast<const int64_t *>(id_array.buffers[1]);

	for (idx_t i = 0; i < to_output; i++) {
		idx_t array_idx = gstate.batch_offset + i + batch.offset;
		int64_t query_id = query_ids[array_idx + id_array.offset];

		if (query_id < 0 || static_cast<size_t>(query_id) >= gstate.read_ids.size()) {
			throw IOException("RYpe returned invalid query_id %lld (expected 0-%zu)", static_cast<long long>(query_id),
			                  gstate.read_ids.size() - 1);
		}
		FlatVector::GetData<string_t>(output.data[0])[i] =
		    StringVector::AddString(output.data[0], gstate.read_ids[query_id]);
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
	data->types = {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UBIGINT),
	               LogicalType::LIST(LogicalType::UBIGINT)};

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

	ArrowArrayStream *input_stream = nullptr;
	unique_ptr<ResultArrowArrayStreamWrapper> input_wrapper;
	auto gstate = BuildExtractionInputStream(context, bind_data, &input_stream, input_wrapper);

	int result = rype_extract_minimizer_set_arrow(input_stream, bind_data.k, bind_data.w, bind_data.salt,
	                                              &gstate->output_stream);
	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe minimizer set extraction failed: %s", err ? err : "unknown error");
	}

	// RYpe took ownership of input_stream — release unique_ptr to avoid double-free
	(void)input_wrapper.release();

	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	ArrowTableFunction::PopulateArrowTableSchema(DBConfig::GetConfig(context), gstate->arrow_table,
	                                             gstate->output_schema);

	// Verify RYpe's output schema matches our declared columns (id + 2 list columns)
	if (gstate->arrow_table.GetColumns().size() != bind_data.names.size()) {
		throw IOException("RYpe returned %zu columns, expected %zu", gstate->arrow_table.GetColumns().size(),
		                  bind_data.names.size());
	}

	gstate->schema_initialized = true;
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
	data->types = {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UBIGINT),
	               LogicalType::LIST(LogicalType::UBIGINT), LogicalType::LIST(LogicalType::UBIGINT),
	               LogicalType::LIST(LogicalType::UBIGINT)};

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

	ArrowArrayStream *input_stream = nullptr;
	unique_ptr<ResultArrowArrayStreamWrapper> input_wrapper;
	auto gstate = BuildExtractionInputStream(context, bind_data, &input_stream, input_wrapper);

	int result = rype_extract_strand_minimizers_arrow(input_stream, bind_data.k, bind_data.w, bind_data.salt,
	                                                  &gstate->output_stream);
	if (result != 0) {
		const char *err = rype_get_last_error();
		throw IOException("RYpe strand minimizer extraction failed: %s", err ? err : "unknown error");
	}

	(void)input_wrapper.release();

	if (gstate->output_stream.get_schema(&gstate->output_stream, &gstate->output_schema) != 0) {
		const char *err = gstate->output_stream.get_last_error(&gstate->output_stream);
		throw IOException("Failed to get RYpe output schema: %s", err ? err : "unknown error");
	}

	ArrowTableFunction::PopulateArrowTableSchema(DBConfig::GetConfig(context), gstate->arrow_table,
	                                             gstate->output_schema);

	// Verify RYpe's output schema matches our declared columns (id + 4 list columns)
	if (gstate->arrow_table.GetColumns().size() != bind_data.names.size()) {
		throw IOException("RYpe returned %zu columns, expected %zu", gstate->arrow_table.GetColumns().size(),
		                  bind_data.names.size());
	}

	gstate->schema_initialized = true;
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
	return tf;
}

void RypeExtractStrandMinimizersTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
