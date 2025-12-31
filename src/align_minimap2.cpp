#include "align_minimap2.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

// Batch size for reading queries
static constexpr idx_t QUERY_BATCH_SIZE = 1024;

unique_ptr<FunctionData> AlignMinimap2TableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                          vector<LogicalType> &return_types,
                                                          vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required positional parameters: query_table, subject_table
	if (input.inputs.size() < 2) {
		throw BinderException("align_minimap2 requires query_table and subject_table parameters");
	}

	data->query_table = input.inputs[0].ToString();
	data->subject_table = input.inputs[1].ToString();

	// Validate query and subject tables/views exist
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, true /* allow_paired */);
	ValidateSequenceTableSchema(context, data->subject_table, false /* allow_paired */);

	// Parse optional named parameters
	auto per_subject_param = input.named_parameters.find("per_subject_database");
	if (per_subject_param != input.named_parameters.end() && !per_subject_param->second.IsNull()) {
		data->per_subject_database = per_subject_param->second.GetValue<bool>();
	}

	auto preset_param = input.named_parameters.find("preset");
	if (preset_param != input.named_parameters.end() && !preset_param->second.IsNull()) {
		data->config.preset = preset_param->second.ToString();
	}

	auto max_secondary_param = input.named_parameters.find("max_secondary");
	if (max_secondary_param != input.named_parameters.end() && !max_secondary_param->second.IsNull()) {
		data->config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
	}

	auto k_param = input.named_parameters.find("k");
	if (k_param != input.named_parameters.end() && !k_param->second.IsNull()) {
		data->config.k = k_param->second.GetValue<int32_t>();
	}

	auto w_param = input.named_parameters.find("w");
	if (w_param != input.named_parameters.end() && !w_param->second.IsNull()) {
		data->config.w = w_param->second.GetValue<int32_t>();
	}

	auto eqx_param = input.named_parameters.find("eqx");
	if (eqx_param != input.named_parameters.end() && !eqx_param->second.IsNull()) {
		data->config.eqx = eqx_param->second.GetValue<bool>();
	}

	// Pre-load all subjects at bind time (required for indexing)
	data->subjects = ReadSubjectTable(context, data->subject_table);

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMinimap2TableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Create aligner with config
	gstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);

	// Build index if not per-subject mode
	if (!data.per_subject_database) {
		gstate->aligner->build_index(data.subjects);
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMinimap2TableFunction::InitLocal(ExecutionContext &context,
                                                                           TableFunctionInitInput &input,
                                                                           GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// Helper functions to set result vectors (reused from read_alignments pattern)
static void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values, idx_t offset,
                                  idx_t count) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = StringVector::AddString(result_vector, values[offset + j]);
	}
}

static void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values, idx_t offset,
                                          idx_t count) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(count);

	for (idx_t j = 0; j < count; j++) {
		result_data[j] = StringVector::AddString(result_vector, values[offset + j]);
		if (!values[offset + j].empty()) {
			validity.SetValid(j);
		}
	}
}

static void SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values, idx_t offset, idx_t count) {
	auto result_data = FlatVector::GetData<uint8_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

static void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values, idx_t offset,
                                  idx_t count) {
	auto result_data = FlatVector::GetData<uint16_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

static void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values, idx_t offset, idx_t count) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

static void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values, idx_t offset,
                                         idx_t count) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(count);

	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
		// Tags return -1 when not present
		if (values[offset + j] != -1) {
			validity.SetValid(j);
		}
	}
}

void AlignMinimap2TableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	std::lock_guard<std::mutex> lock(global_state.lock);

	// Check if we're done
	if (global_state.done) {
		output.SetCardinality(0);
		return;
	}

	// Determine how many results we can output this call
	idx_t available = global_state.result_buffer.size() - global_state.buffer_offset;

	// If buffer is empty or exhausted, fill it with more alignments
	while (available == 0) {
		// Clear buffer for new batch
		global_state.result_buffer.clear();
		global_state.buffer_offset = 0;

		if (bind_data.per_subject_database) {
			// Per-subject mode: load all queries once, then align against each subject

			// Load all queries into memory on first use (avoids re-reading table for each subject)
			if (!global_state.queries_loaded) {
				idx_t query_offset = 0;
				miint::SequenceRecordBatch batch;
				bool has_more = true;
				while (has_more) {
					batch.clear();
					has_more = ReadQueryBatch(context, bind_data.query_table, bind_data.query_schema,
					                          QUERY_BATCH_SIZE, query_offset, batch);
					// Append batch to all_queries
					for (size_t i = 0; i < batch.size(); i++) {
						global_state.all_queries.read_ids.push_back(std::move(batch.read_ids[i]));
						global_state.all_queries.comments.push_back(std::move(batch.comments[i]));
						global_state.all_queries.sequences1.push_back(std::move(batch.sequences1[i]));
						global_state.all_queries.quals1.push_back(std::move(batch.quals1[i]));
						if (batch.is_paired) {
							global_state.all_queries.sequences2.push_back(std::move(batch.sequences2[i]));
							global_state.all_queries.quals2.push_back(std::move(batch.quals2[i]));
						}
					}
					global_state.all_queries.is_paired = batch.is_paired;
				}
				global_state.queries_loaded = true;
			}

			// Check if we've processed all subjects
			if (global_state.current_subject_idx >= bind_data.subjects.size()) {
				global_state.done = true;
				output.SetCardinality(0);
				return;
			}

			// Build index for current subject
			global_state.aligner->build_single_index(bind_data.subjects[global_state.current_subject_idx]);

			// Align all queries against this subject (queries already in memory)
			if (!global_state.all_queries.empty()) {
				global_state.aligner->align(global_state.all_queries, global_state.result_buffer);
			}

			// Move to next subject
			global_state.current_subject_idx++;

		} else {
			// Standard mode: stream queries against single index
			miint::SequenceRecordBatch query_batch;
			bool has_more = ReadQueryBatch(context, bind_data.query_table, bind_data.query_schema, QUERY_BATCH_SIZE,
			                               global_state.current_query_offset, query_batch);

			if (query_batch.empty() && !has_more) {
				global_state.done = true;
				output.SetCardinality(0);
				return;
			}

			if (!query_batch.empty()) {
				global_state.aligner->align(query_batch, global_state.result_buffer);
			}
		}

		available = global_state.result_buffer.size() - global_state.buffer_offset;
	}

	// Output up to STANDARD_VECTOR_SIZE results
	idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	idx_t offset = global_state.buffer_offset;
	auto &batch = global_state.result_buffer;

	// Set result vectors
	idx_t field_idx = 0;
	SetResultVectorString(output.data[field_idx++], batch.read_ids, offset, output_count);
	SetResultVectorUInt16(output.data[field_idx++], batch.flags, offset, output_count);
	SetResultVectorString(output.data[field_idx++], batch.references, offset, output_count);
	SetResultVectorInt64(output.data[field_idx++], batch.positions, offset, output_count);
	SetResultVectorInt64(output.data[field_idx++], batch.stop_positions, offset, output_count);
	SetResultVectorUInt8(output.data[field_idx++], batch.mapqs, offset, output_count);
	SetResultVectorString(output.data[field_idx++], batch.cigars, offset, output_count);
	SetResultVectorString(output.data[field_idx++], batch.mate_references, offset, output_count);
	SetResultVectorInt64(output.data[field_idx++], batch.mate_positions, offset, output_count);
	SetResultVectorInt64(output.data[field_idx++], batch.template_lengths, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_as_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xs_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_ys_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xn_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xm_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xo_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xg_values, offset, output_count);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_nm_values, offset, output_count);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_yt_values, offset, output_count);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_md_values, offset, output_count);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_sa_values, offset, output_count);

	output.SetCardinality(output_count);
	global_state.buffer_offset += output_count;
}

TableFunction AlignMinimap2TableFunction::GetFunction() {
	auto tf = TableFunction("align_minimap2", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                        InitLocal);

	// Named parameters
	tf.named_parameters["per_subject_database"] = LogicalType::BOOLEAN;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["w"] = LogicalType::INTEGER;
	tf.named_parameters["eqx"] = LogicalType::BOOLEAN;

	return tf;
}

void AlignMinimap2TableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
