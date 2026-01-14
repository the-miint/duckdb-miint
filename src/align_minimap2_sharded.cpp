#include "align_minimap2_sharded.hpp"
#include "align_minimap2_common.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"

namespace duckdb {

// Helper to validate read_to_shard table schema
static void ValidateReadToShardSchema(ClientContext &context, const std::string &table_name) {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("Table or view '%s' does not exist", table_name);
	}

	vector<string> col_names;
	vector<LogicalType> col_types;

	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			col_names.push_back(col.Name());
			col_types.push_back(col.Type());
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		col_names = view.names;
		col_types = view.types;
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}

	// Build case-insensitive name map
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Check required columns
	auto it_read_id = name_to_idx.find("read_id");
	if (it_read_id == name_to_idx.end()) {
		throw BinderException("read_to_shard table '%s' missing required column 'read_id'", table_name);
	}
	if (col_types[it_read_id->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'read_id' in read_to_shard table '%s' must be VARCHAR", table_name);
	}

	auto it_shard = name_to_idx.find("shard_name");
	if (it_shard == name_to_idx.end()) {
		throw BinderException("read_to_shard table '%s' missing required column 'shard_name'", table_name);
	}
	if (col_types[it_shard->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'shard_name' in read_to_shard table '%s' must be VARCHAR", table_name);
	}
}

// Helper to read shard counts from read_to_shard table
// Errors on NULL shard_name values
static std::vector<ShardInfo> ReadShardCounts(ClientContext &context, const std::string &table_name,
                                               const std::string &shard_directory, FileSystem &fs) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Query shard counts ordered by count descending (largest first)
	std::string query = "SELECT shard_name, COUNT(*) as cnt FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(table_name) +
	                    " GROUP BY shard_name ORDER BY cnt DESC";

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read shard counts from '%s': %s", table_name,
		                            query_result->GetError());
	}

	std::vector<ShardInfo> shards;
	auto &materialized = query_result->Cast<MaterializedQueryResult>();

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &shard_name_vec = chunk->data[0];
		auto &count_vec = chunk->data[1];

		UnifiedVectorFormat shard_data, count_data;
		shard_name_vec.ToUnifiedFormat(chunk->size(), shard_data);
		count_vec.ToUnifiedFormat(chunk->size(), count_data);

		auto shard_names = UnifiedVectorFormat::GetData<string_t>(shard_data);
		auto counts = UnifiedVectorFormat::GetData<int64_t>(count_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			auto shard_idx = shard_data.sel->get_index(i);
			auto count_idx = count_data.sel->get_index(i);

			// Error on NULL shard_name
			if (!shard_data.validity.RowIsValid(shard_idx)) {
				throw BinderException("read_to_shard table '%s' contains NULL shard_name values", table_name);
			}

			ShardInfo info;
			info.name = shard_names[shard_idx].GetString();
			info.read_count = static_cast<idx_t>(counts[count_idx]);

			// Build index path
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
	}

	if (shards.empty()) {
		throw BinderException("read_to_shard table '%s' is empty or has no valid shard names", table_name);
	}

	return shards;
}

// Helper to filter out unmapped reads from result batch (in-place)
static void FilterMappedOnly(miint::SAMRecordBatch &batch) {
	if (batch.empty()) {
		return;
	}

	// In-place compaction: use write_idx to track where to write mapped reads
	idx_t write_idx = 0;
	for (idx_t read_idx = 0; read_idx < batch.size(); read_idx++) {
		// Filter out unmapped reads (flag & 0x4)
		if ((batch.flags[read_idx] & 0x4) != 0) {
			continue;
		}

		// Move element to write position if needed
		if (write_idx != read_idx) {
			batch.read_ids[write_idx] = std::move(batch.read_ids[read_idx]);
			batch.flags[write_idx] = batch.flags[read_idx];
			batch.references[write_idx] = std::move(batch.references[read_idx]);
			batch.positions[write_idx] = batch.positions[read_idx];
			batch.stop_positions[write_idx] = batch.stop_positions[read_idx];
			batch.mapqs[write_idx] = batch.mapqs[read_idx];
			batch.cigars[write_idx] = std::move(batch.cigars[read_idx]);
			batch.mate_references[write_idx] = std::move(batch.mate_references[read_idx]);
			batch.mate_positions[write_idx] = batch.mate_positions[read_idx];
			batch.template_lengths[write_idx] = batch.template_lengths[read_idx];
			batch.tag_as_values[write_idx] = batch.tag_as_values[read_idx];
			batch.tag_xs_values[write_idx] = batch.tag_xs_values[read_idx];
			batch.tag_ys_values[write_idx] = batch.tag_ys_values[read_idx];
			batch.tag_xn_values[write_idx] = batch.tag_xn_values[read_idx];
			batch.tag_xm_values[write_idx] = batch.tag_xm_values[read_idx];
			batch.tag_xo_values[write_idx] = batch.tag_xo_values[read_idx];
			batch.tag_xg_values[write_idx] = batch.tag_xg_values[read_idx];
			batch.tag_nm_values[write_idx] = batch.tag_nm_values[read_idx];
			batch.tag_yt_values[write_idx] = std::move(batch.tag_yt_values[read_idx]);
			batch.tag_md_values[write_idx] = std::move(batch.tag_md_values[read_idx]);
			batch.tag_sa_values[write_idx] = std::move(batch.tag_sa_values[read_idx]);
		}
		write_idx++;
	}

	// Truncate vectors to final size
	batch.read_ids.resize(write_idx);
	batch.flags.resize(write_idx);
	batch.references.resize(write_idx);
	batch.positions.resize(write_idx);
	batch.stop_positions.resize(write_idx);
	batch.mapqs.resize(write_idx);
	batch.cigars.resize(write_idx);
	batch.mate_references.resize(write_idx);
	batch.mate_positions.resize(write_idx);
	batch.template_lengths.resize(write_idx);
	batch.tag_as_values.resize(write_idx);
	batch.tag_xs_values.resize(write_idx);
	batch.tag_ys_values.resize(write_idx);
	batch.tag_xn_values.resize(write_idx);
	batch.tag_xm_values.resize(write_idx);
	batch.tag_xo_values.resize(write_idx);
	batch.tag_xg_values.resize(write_idx);
	batch.tag_nm_values.resize(write_idx);
	batch.tag_yt_values.resize(write_idx);
	batch.tag_md_values.resize(write_idx);
	batch.tag_sa_values.resize(write_idx);
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
	data->shards = ReadShardCounts(context, data->read_to_shard_table, data->shard_directory, fs);

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

unique_ptr<LocalTableFunctionState> AlignMinimap2ShardedTableFunction::InitLocal(ExecutionContext &context,
                                                                                   TableFunctionInitInput &input,
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

		// Use STANDARD_VECTOR_SIZE as query batch size to bound result buffer memory
		// (each query can produce multiple alignments)
		bool has_more = ReadShardQueryBatch(context, bind_data.query_table, bind_data.read_to_shard_table,
		                                     shard.name, bind_data.query_schema, STANDARD_VECTOR_SIZE,
		                                     local_state.current_read_offset, query_batch);

		if (query_batch.empty() && !has_more) {
			// Shard exhausted - this thread is done (no work stealing)
			// Each thread owns its shard for the lifetime of the query
			// This ensures each index is loaded exactly once
			output.SetCardinality(0);
			return;
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
