#include "align_minimap2_sharded.hpp"
#include "align_common.hpp"
#include "duckdb/common/file_system.hpp"

namespace duckdb {

// Build minimap2 ShardInfo from raw shard name/counts
// Validates index files exist and are valid minimap2 indexes
static std::vector<ShardInfo> BuildMinimap2ShardInfos(ClientContext &context, const std::string &table_name,
                                                      const std::string &shard_directory, FileSystem &fs) {
	// Get raw shard names and counts from shared utility
	auto raw_shards = ReadShardNameCounts(context, table_name);

	std::vector<ShardInfo> shards;
	shards.reserve(raw_shards.size());

	for (const auto &raw : raw_shards) {
		ShardInfo info;
		info.name = raw.name;
		info.read_count = raw.count;

		// Build index path: shard_directory/shard_name.mmi
		info.index_path = shard_directory;
		if (!info.index_path.empty() && info.index_path.back() != '/') {
			info.index_path += '/';
		}
		info.index_path += info.name + ".mmi";

		// Fail fast: check if .mmi file exists
		if (!fs.FileExists(info.index_path)) {
			throw BinderException("Shard index file does not exist: %s", info.index_path);
		}

		// Validate it's a valid minimap2 index
		if (!miint::Minimap2Aligner::is_index_file(info.index_path)) {
			throw BinderException("File is not a valid minimap2 index: %s", info.index_path);
		}

		shards.push_back(std::move(info));
	}

	return shards;
}

unique_ptr<FunctionData> AlignMinimap2ShardedTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                 vector<LogicalType> &return_types,
                                                                 vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required: query_table (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("align_minimap2_sharded requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	// Required: shard_directory named parameter
	auto shard_dir_param = input.named_parameters.find("shard_directory");
	if (shard_dir_param == input.named_parameters.end() || shard_dir_param->second.IsNull()) {
		throw BinderException("align_minimap2_sharded requires shard_directory parameter");
	}
	data->shard_directory = shard_dir_param->second.ToString();

	// Required: read_to_shard named parameter
	auto read_to_shard_param = input.named_parameters.find("read_to_shard");
	if (read_to_shard_param == input.named_parameters.end() || read_to_shard_param->second.IsNull()) {
		throw BinderException("align_minimap2_sharded requires read_to_shard parameter");
	}
	data->read_to_shard_table = read_to_shard_param->second.ToString();

	// Validate shard_directory exists
	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.DirectoryExists(data->shard_directory)) {
		throw BinderException("Shard directory does not exist: %s", data->shard_directory);
	}

	// Validate query table/view exists
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, true /* allow_paired */);

	// Validate read_to_shard table schema
	ValidateReadToShardSchema(context, data->read_to_shard_table);

	// Parse minimap2 config parameters (preset, max_secondary, eqx)
	// Always warn about k/w since we use pre-built indexes
	ParseMinimap2ConfigParams(input.named_parameters, data->config, true /* warn_prebuilt_index */);

	// Read shard counts and validate .mmi files exist (fail fast)
	data->shards = BuildMinimap2ShardInfos(context, data->read_to_shard_table, data->shard_directory, fs);

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMinimap2ShardedTableFunction::InitGlobal(ClientContext &context,
                                                                                   TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->shard_count = data.shards.size();
	return gstate;
}

unique_ptr<LocalTableFunctionState>
AlignMinimap2ShardedTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                             GlobalTableFunctionState *global_state) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	// Create per-thread aligner with config
	lstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
	return lstate;
}

void AlignMinimap2ShardedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	while (true) {
		// Check if we have buffered results to output
		idx_t available = local_state.result_buffer.size() - local_state.buffer_offset;

		if (available > 0) {
			// Output up to STANDARD_VECTOR_SIZE results
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, local_state.result_buffer, local_state.buffer_offset, output_count);
			local_state.buffer_offset += output_count;
			return;
		}

		// Buffer is empty, need to get more results

		// Claim a shard if we don't have one
		if (!local_state.has_shard) {
			std::lock_guard<std::mutex> lock(global_state.lock);
			if (global_state.next_shard_idx >= bind_data.shards.size()) {
				// No more shards to process
				output.SetCardinality(0);
				return;
			}
			local_state.current_shard_idx = global_state.next_shard_idx++;
			local_state.has_shard = true;
			local_state.current_read_offset = 0;

			// Load index for new shard
			auto &shard = bind_data.shards[local_state.current_shard_idx];
			try {
				local_state.aligner->load_index(shard.index_path);
			} catch (const std::exception &e) {
				throw IOException("Failed to load minimap2 index from '%s': %s", shard.index_path, e.what());
			}
		}

		// Read next batch of queries for current shard
		auto &shard = bind_data.shards[local_state.current_shard_idx];
		miint::SequenceRecordBatch query_batch;

		// Use large batch size to reduce per-batch overhead (each batch creates a
		// new Connection and re-executes a JOIN query against the shard table)
		bool has_more = ReadShardQueryBatch(context, bind_data.query_table, bind_data.read_to_shard_table, shard.name,
		                                    bind_data.query_schema, SHARDED_QUERY_BATCH_SIZE,
		                                    local_state.current_read_offset, query_batch);

		if (query_batch.empty() && !has_more) {
			// Shard exhausted - release and claim next available shard
			local_state.has_shard = false;
			continue;
		}

		// Align batch
		local_state.result_buffer.clear();
		local_state.buffer_offset = 0;

		if (!query_batch.empty()) {
			local_state.aligner->align(query_batch, local_state.result_buffer);
			// Filter out unmapped reads
			FilterMappedOnly(local_state.result_buffer);
		}

		// Loop back to output results
	}
}

TableFunction AlignMinimap2ShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_minimap2_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["eqx"] = LogicalType::BOOLEAN;

	return tf;
}

void AlignMinimap2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
