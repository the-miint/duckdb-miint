#include "align_minimap2.hpp"
#include "align_common.hpp"
#include "shard_debug.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

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
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);

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

	// Parse debug parameter
	auto debug_param = input.named_parameters.find("debug");
	if (debug_param != input.named_parameters.end() && !debug_param->second.IsNull()) {
		data->debug = debug_param->second.GetValue<bool>();
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
		// Traditional mode: load subjects (ReadSubjectTable validates internally)
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
	gstate->debug = data.debug;
	gstate->start_time = std::chrono::steady_clock::now();

	if (data.per_subject_database) {
		// Per-subject mode: single-threaded, builds index per-subject in Execute()
		gstate->per_subject_mode = true;
		gstate->num_threads = 1;
		auto ps = std::make_unique<PerSubjectModeState>();
		ps->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
		gstate->per_subject = std::move(ps);
		SHARD_DBG(*gstate, "InitGlobal: per_subject_mode, num_threads=1");
	} else if (data.using_prebuilt_index()) {
		// Prebuilt index: load into SharedMinimap2Index for multi-threaded access
		auto st = std::make_unique<StandardModeState>();
		try {
			st->shared_index = std::make_shared<miint::SharedMinimap2Index>(data.index_path, data.config);
		} catch (const std::exception &e) {
			throw IOException("Failed to load minimap2 index from '%s': %s", data.index_path, e.what());
		}
		gstate->standard = std::move(st);
		gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		SHARD_DBG_MEM(*gstate, "InitGlobal: prebuilt index loaded, MaxThreads()=%zu",
		              static_cast<size_t>(gstate->num_threads));
	} else {
		// Standard mode with subject table: build shared index
		auto st = std::make_unique<StandardModeState>();
		st->shared_index = miint::Minimap2Aligner::BuildSharedIndex(data.subjects, data.config);
		gstate->standard = std::move(st);
		gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		SHARD_DBG_MEM(*gstate, "InitGlobal: shared index built from %zu subjects, MaxThreads()=%zu",
		              data.subjects.size(), static_cast<size_t>(gstate->num_threads));
	}

	// Materialize all query read_ids for standard mode (avoids slow OFFSET pagination)
	if (!gstate->per_subject_mode) {
		auto id_start = std::chrono::steady_clock::now();
		auto &db = DatabaseInstance::GetDatabase(context);
		Connection conn(db);
		std::string query = "SELECT read_id FROM " + KeywordHelper::WriteOptionallyQuoted(data.query_table);
		auto query_result = conn.Query(query);
		if (query_result->HasError()) {
			throw InvalidInputException("Failed to read query IDs from '%s': %s", data.query_table,
			                            query_result->GetError());
		}
		auto &materialized = query_result->Cast<MaterializedQueryResult>();
		auto &ids = gstate->standard->all_query_ids;
		while (true) {
			auto chunk = materialized.Fetch();
			if (!chunk || chunk->size() == 0) {
				break;
			}
			UnifiedVectorFormat id_data;
			chunk->data[0].ToUnifiedFormat(chunk->size(), id_data);
			auto id_strings = UnifiedVectorFormat::GetData<string_t>(id_data);
			for (idx_t i = 0; i < chunk->size(); i++) {
				auto idx = id_data.sel->get_index(i);
				if (id_data.validity.RowIsValid(idx)) {
					ids.push_back(id_strings[idx].GetString());
				}
			}
		}
		auto id_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - id_start).count();
		SHARD_DBG(*gstate, "InitGlobal: materialized %zu query IDs in %ldms", ids.size(), static_cast<long>(id_ms));
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMinimap2TableFunction::InitLocal(ExecutionContext &context,
                                                                          TableFunctionInitInput &input,
                                                                          GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();

	if (!gstate.per_subject_mode) {
		// Standard mode: create per-thread aligner and attach shared index
		auto &data = input.bind_data->Cast<Data>();
		lstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
		lstate->aligner->attach_shared_index(gstate.standard->shared_index);
		auto thread_num = gstate.init_local_count.fetch_add(1) + 1;
		SHARD_DBG(gstate, "InitLocal: thread %zu of %zu initialized", static_cast<size_t>(thread_num),
		          static_cast<size_t>(gstate.num_threads));
	}

	return lstate;
}

// Per-subject mode: single-threaded with global mutex
static void ExecutePerSubject(ClientContext &context, const AlignMinimap2TableFunction::Data &bind_data,
                              AlignMinimap2TableFunction::PerSubjectModeState &ps, DataChunk &output) {
	std::lock_guard<std::mutex> lock(ps.lock);

	if (ps.done) {
		output.SetCardinality(0);
		return;
	}

	idx_t available = ps.result_buffer.size() - ps.buffer_offset;

	while (available == 0) {
		ps.result_buffer.clear();
		ps.buffer_offset = 0;

		// Load all queries into memory on first use
		if (!ps.queries_loaded) {
			ps.all_queries.is_paired = bind_data.query_schema.has_sequence2;
			idx_t query_offset = 0;
			miint::SequenceRecordBatch batch;
			bool has_more = true;
			while (has_more) {
				batch.clear();
				has_more = ReadQueryBatch(context, bind_data.query_table, bind_data.query_schema,
				                          MINIMAP2_QUERY_BATCH_SIZE, query_offset, batch);
				for (size_t i = 0; i < batch.size(); i++) {
					ps.all_queries.read_ids.push_back(std::move(batch.read_ids[i]));
					ps.all_queries.comments.push_back(std::move(batch.comments[i]));
					ps.all_queries.sequences1.push_back(std::move(batch.sequences1[i]));
					ps.all_queries.quals1.push_back(std::move(batch.quals1[i]));
					if (ps.all_queries.is_paired) {
						ps.all_queries.sequences2.push_back(std::move(batch.sequences2[i]));
						ps.all_queries.quals2.push_back(std::move(batch.quals2[i]));
					}
				}
			}
			ps.queries_loaded = true;
		}

		if (ps.current_subject_idx >= bind_data.subjects.size()) {
			ps.done = true;
			output.SetCardinality(0);
			return;
		}

		ps.aligner->build_single_index(bind_data.subjects[ps.current_subject_idx]);

		if (!ps.all_queries.empty()) {
			ps.aligner->align(ps.all_queries, ps.result_buffer);
		}

		ps.current_subject_idx++;
		available = ps.result_buffer.size() - ps.buffer_offset;
	}

	idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	OutputSAMRecordBatch(output, ps.result_buffer, ps.buffer_offset, output_count);
	ps.buffer_offset += output_count;
}

// Standard mode: multi-threaded with per-thread aligner and atomic batch claiming
static void ExecuteStandard(ClientContext &context, const AlignMinimap2TableFunction::Data &bind_data,
                            AlignMinimap2TableFunction::GlobalState &gstate,
                            AlignMinimap2TableFunction::LocalState &lstate, DataChunk &output) {
	auto &st = *gstate.standard;
	idx_t id_count = st.all_query_ids.size();

	while (true) {
		// 1. If local result_buffer has data, output a chunk
		idx_t available = lstate.result_buffer.size() - lstate.buffer_offset;
		if (available > 0) {
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, lstate.result_buffer, lstate.buffer_offset, output_count);
			lstate.buffer_offset += output_count;
			return;
		}

		// 2. Buffer exhausted — claim a new batch atomically
		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		idx_t my_offset = st.next_query_offset.fetch_add(st.batch_size);

		if (my_offset >= id_count) {
			SHARD_DBG(gstate, "ExecuteStandard: DONE (offset %zu >= id_count %zu)", static_cast<size_t>(my_offset),
			          static_cast<size_t>(id_count));
			output.SetCardinality(0);
			return;
		}

		idx_t batch_count = std::min(st.batch_size, id_count - my_offset);
		SHARD_DBG(gstate, "ExecuteStandard: claimed batch at offset %zu count %zu", static_cast<size_t>(my_offset),
		          static_cast<size_t>(batch_count));

		// 3. Read sequences by joining a temp table of IDs against the query table
		miint::SequenceRecordBatch query_batch;
		auto read_start = std::chrono::steady_clock::now();
		ReadBatchByIds(context, bind_data.query_table, bind_data.query_schema, st.all_query_ids, my_offset, batch_count,
		               query_batch);
		auto read_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - read_start)
		        .count();
		SHARD_DBG(gstate, "ExecuteStandard: read %zu queries in %ldms (offset=%zu)", query_batch.size(),
		          static_cast<long>(read_ms), static_cast<size_t>(my_offset));

		if (query_batch.empty()) {
			continue;
		}

		// 4. Align using per-thread aligner
		auto align_start = std::chrono::steady_clock::now();
		lstate.aligner->align(query_batch, lstate.result_buffer);
		auto align_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
		        .count();
		SHARD_DBG(gstate, "ExecuteStandard: aligned %zu queries -> %zu results in %ldms", query_batch.size(),
		          lstate.result_buffer.size(), static_cast<long>(align_ms));
	}
}

void AlignMinimap2TableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.per_subject_mode) {
		ExecutePerSubject(context, bind_data, *gstate.per_subject, output);
	} else {
		auto &lstate = data_p.local_state->Cast<LocalState>();
		ExecuteStandard(context, bind_data, gstate, lstate, output);
	}
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
	tf.named_parameters["debug"] = LogicalType::BOOLEAN;

	return tf;
}

void AlignMinimap2TableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
