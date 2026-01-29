#include "align_bowtie2_sharded.hpp"
#include "align_common.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/printer.hpp"

namespace duckdb {

// Build bowtie2 ShardInfo from raw shard name/counts
// Validates index files exist and are valid bowtie2 indexes
static std::vector<Bowtie2ShardInfo> BuildBowtie2ShardInfos(ClientContext &context, const std::string &table_name,
                                                            const std::string &shard_directory) {
	// Get raw shard names and counts from shared utility
	auto raw_shards = ReadShardNameCounts(context, table_name);

	std::vector<Bowtie2ShardInfo> shards;
	shards.reserve(raw_shards.size());

	for (const auto &raw : raw_shards) {
		Bowtie2ShardInfo info;
		info.name = raw.name;
		info.read_count = raw.count;

		// Build index prefix: shard_directory/shard_name/index
		info.index_prefix = shard_directory;
		if (!info.index_prefix.empty() && info.index_prefix.back() != '/') {
			info.index_prefix += '/';
		}
		info.index_prefix += info.name + "/index";

		// Fail fast: check if bowtie2 index files exist
		if (!miint::Bowtie2Aligner::is_index_prefix(info.index_prefix)) {
			throw BinderException("No valid bowtie2 index found at prefix: %s. "
			                      "Expected files like %s.1.bt2, %s.rev.1.bt2, etc.",
			                      info.index_prefix, info.index_prefix, info.index_prefix);
		}

		shards.push_back(std::move(info));
	}

	return shards;
}

unique_ptr<FunctionData> AlignBowtie2ShardedTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                vector<LogicalType> &return_types,
                                                                vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required: query_table (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("align_bowtie2_sharded requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	// Required: shard_directory named parameter
	auto shard_dir_param = input.named_parameters.find("shard_directory");
	if (shard_dir_param == input.named_parameters.end() || shard_dir_param->second.IsNull()) {
		throw BinderException("align_bowtie2_sharded requires shard_directory parameter");
	}
	data->shard_directory = shard_dir_param->second.ToString();

	// Required: read_to_shard named parameter
	auto read_to_shard_param = input.named_parameters.find("read_to_shard");
	if (read_to_shard_param == input.named_parameters.end() || read_to_shard_param->second.IsNull()) {
		throw BinderException("align_bowtie2_sharded requires read_to_shard parameter");
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

	// Parse bowtie2 config parameters
	ParseBowtie2ConfigParams(input.named_parameters, data->config);

	// In sharded mode, parallelism comes from multiple bowtie2 processes (one per shard),
	// so each bowtie2 process uses a single thread to avoid CPU oversubscription.
	// Warn user if they specified a threads value other than 1.
	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		int32_t threads_val = threads_param->second.GetValue<int32_t>();
		if (threads_val != 1) {
			Printer::Print("WARNING: Parameter 'threads' is ignored in sharded mode. "
			               "Parallelism comes from multiple single-threaded bowtie2 processes "
			               "(one per shard). Use DuckDB's SET threads=N to control shard parallelism.\n");
		}
	}
	data->config.threads = 1;

	// Read shard counts and validate index files exist (fail fast)
	data->shards = BuildBowtie2ShardInfos(context, data->read_to_shard_table, data->shard_directory);

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignBowtie2ShardedTableFunction::InitGlobal(ClientContext &context,
                                                                                  TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->shard_count = data.shards.size();
	return gstate;
}

unique_ptr<LocalTableFunctionState>
AlignBowtie2ShardedTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                            GlobalTableFunctionState *global_state) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	// Create per-thread aligner with config
	lstate->aligner = std::make_unique<miint::Bowtie2Aligner>(data.config);
	return lstate;
}

void AlignBowtie2ShardedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
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
			local_state.finished_aligning = false;

			// Load index for new shard
			auto &shard = bind_data.shards[local_state.current_shard_idx];
			try {
				local_state.aligner->load_index(shard.index_prefix);
			} catch (const std::exception &e) {
				throw IOException("Failed to load bowtie2 index from '%s': %s", shard.index_prefix, e.what());
			}
		}

		// If we've finished aligning for this shard, move to next
		if (local_state.finished_aligning) {
			// Reset state for claiming a new shard
			local_state.has_shard = false;
			// Reset aligner for reuse with next shard
			local_state.aligner->reset();
			continue;
		}

		// Read next batch of queries for current shard
		auto &shard = bind_data.shards[local_state.current_shard_idx];
		miint::SequenceRecordBatch query_batch;

		// Use STANDARD_VECTOR_SIZE as query batch size to bound result buffer memory
		bool has_more = ReadShardQueryBatch(context, bind_data.query_table, bind_data.read_to_shard_table, shard.name,
		                                    bind_data.query_schema, STANDARD_VECTOR_SIZE,
		                                    local_state.current_read_offset, query_batch);

		// Clear buffer for new results
		local_state.result_buffer.clear();
		local_state.buffer_offset = 0;

		if (query_batch.empty() && !has_more) {
			// Shard queries exhausted - call finish() to get remaining results from bowtie2
			local_state.aligner->finish(local_state.result_buffer);
			local_state.finished_aligning = true;

			// Filter out unmapped reads
			FilterMappedOnly(local_state.result_buffer);

			// If still no results after finish(), move to next shard
			if (local_state.result_buffer.empty()) {
				local_state.has_shard = false;
				local_state.aligner->reset();
				continue;
			}
		} else if (!query_batch.empty()) {
			// Align batch
			local_state.aligner->align(query_batch, local_state.result_buffer);
			// Filter out unmapped reads
			FilterMappedOnly(local_state.result_buffer);
		}

		// Loop back to output results (or get more if none)
	}
}

TableFunction AlignBowtie2ShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["local"] = LogicalType::BOOLEAN;
	tf.named_parameters["threads"] = LogicalType::INTEGER;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["extra_args"] = LogicalType::VARCHAR;
	tf.named_parameters["quiet"] = LogicalType::BOOLEAN;

	return tf;
}

void AlignBowtie2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
