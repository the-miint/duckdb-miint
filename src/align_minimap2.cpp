#include "align_minimap2.hpp"
#include "align_common.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

unique_ptr<FunctionData> AlignMinimap2TableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                          vector<LogicalType> &return_types,
                                                          vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required: query_table (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("align_minimap2 requires query_table parameter");
	}

	data->query_table = input.inputs[0].ToString();

	// Parse subject_table named parameter
	auto subject_param = input.named_parameters.find("subject_table");
	if (subject_param != input.named_parameters.end() && !subject_param->second.IsNull()) {
		data->subject_table = subject_param->second.ToString();
	}

	// Parse index_path named parameter
	auto index_path_param = input.named_parameters.find("index_path");
	if (index_path_param != input.named_parameters.end() && !index_path_param->second.IsNull()) {
		data->index_path = index_path_param->second.ToString();
	}

	// VALIDATION: Exactly one of subject_table or index_path must be provided
	bool has_subject = !data->subject_table.empty();
	bool has_index = data->using_prebuilt_index();

	if (!has_subject && !has_index) {
		throw BinderException("align_minimap2 requires either subject_table or index_path parameter");
	}
	if (has_subject && has_index) {
		throw BinderException(
		    "align_minimap2: Cannot specify both subject_table and index_path. "
		    "Use subject_table to build index from sequences, or index_path to load pre-built index.");
	}

	// Validate query table/view exists
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, true /* allow_paired */);

	// Parse optional named parameters
	auto per_subject_param = input.named_parameters.find("per_subject_database");
	if (per_subject_param != input.named_parameters.end() && !per_subject_param->second.IsNull()) {
		data->per_subject_database = per_subject_param->second.GetValue<bool>();

		// Validate per_subject_database incompatible with index_path
		if (data->per_subject_database && data->using_prebuilt_index()) {
			throw BinderException("per_subject_database mode is incompatible with index_path. "
			                      "Pre-built indexes contain all subjects.");
		}
	}

	// Parse minimap2 config parameters (preset, max_secondary, k, w, eqx)
	ParseMinimap2ConfigParams(input.named_parameters, data->config, data->using_prebuilt_index());

	// Handle subject_table vs index_path modes
	if (data->using_prebuilt_index()) {
		// Validate index file exists
		auto &fs = FileSystem::GetFileSystem(context);
		if (!fs.FileExists(data->index_path)) {
			throw BinderException("Index file does not exist: %s", data->index_path);
		}

		// Validate it's a valid minimap2 index
		// NOTE: This is advisory only (TOCTOU) - actual load may still fail
		if (!miint::Minimap2Aligner::is_index_file(data->index_path)) {
			throw BinderException("File is not a valid minimap2 index: %s", data->index_path);
		}

		// Note: subjects vector remains empty in this mode
	} else {
		// Traditional mode: validate subject table and load subjects
		ValidateSequenceTableSchema(context, data->subject_table, false /* allow_paired */);
		data->subjects = ReadSubjectTable(context, data->subject_table);
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

unique_ptr<GlobalTableFunctionState> AlignMinimap2TableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Create aligner with config
	gstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);

	// Load or build index based on mode
	if (data.using_prebuilt_index()) {
		// Load pre-built index from file
		try {
			gstate->aligner->load_index(data.index_path);
		} catch (const std::exception &e) {
			throw IOException("Failed to load minimap2 index from '%s': %s", data.index_path, e.what());
		}
	} else {
		// Traditional mode: build index from subjects
		if (!data.per_subject_database) {
			gstate->aligner->build_index(data.subjects);
		}
		// Note: per_subject mode builds index per-subject in Execute()
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMinimap2TableFunction::InitLocal(ExecutionContext &context,
                                                                           TableFunctionInitInput &input,
                                                                           GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
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
					                          MINIMAP2_QUERY_BATCH_SIZE, query_offset, batch);
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
			bool has_more = ReadQueryBatch(context, bind_data.query_table, bind_data.query_schema, MINIMAP2_QUERY_BATCH_SIZE,
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
	OutputSAMRecordBatch(output, global_state.result_buffer, global_state.buffer_offset, output_count);
	global_state.buffer_offset += output_count;
}

TableFunction AlignMinimap2TableFunction::GetFunction() {
	// Only query_table is a positional parameter
	auto tf = TableFunction("align_minimap2", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["subject_table"] = LogicalType::VARCHAR;
	tf.named_parameters["index_path"] = LogicalType::VARCHAR;
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
