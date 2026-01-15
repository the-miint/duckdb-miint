#pragma once
/*
 * Shared utilities for alignment table functions:
 *   - align_minimap2, align_minimap2_sharded
 *   - align_bowtie2, align_bowtie2_sharded
 *
 * Provides common config parsing, output schema, and result output logic.
 */

#include "Minimap2Aligner.hpp"
#include "Bowtie2Aligner.hpp"
#include "SAMRecord.hpp"
#include "align_result_utils.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"

namespace duckdb {

// Batch size for reading queries (shared across alignment functions)
static constexpr idx_t ALIGNMENT_QUERY_BATCH_SIZE = 1024;
// Backward compatibility alias
static constexpr idx_t MINIMAP2_QUERY_BATCH_SIZE = ALIGNMENT_QUERY_BATCH_SIZE;

// Get the standard alignment output column names
inline std::vector<std::string> GetAlignmentOutputNames() {
	return {"read_id", "flags",          "reference",     "position",        "stop_position", "mapq",
	        "cigar",   "mate_reference", "mate_position", "template_length", "tag_as",        "tag_xs",
	        "tag_ys",  "tag_xn",         "tag_xm",        "tag_xo",          "tag_xg",        "tag_nm",
	        "tag_yt",  "tag_md",         "tag_sa"};
}

// Get the standard alignment output column types
inline std::vector<LogicalType> GetAlignmentOutputTypes() {
	return {LogicalType::VARCHAR,   // read_id
	        LogicalType::USMALLINT, // flags
	        LogicalType::VARCHAR,   // reference
	        LogicalType::BIGINT,    // position
	        LogicalType::BIGINT,    // stop_position
	        LogicalType::UTINYINT,  // mapq
	        LogicalType::VARCHAR,   // cigar
	        LogicalType::VARCHAR,   // mate_reference
	        LogicalType::BIGINT,    // mate_position
	        LogicalType::BIGINT,    // template_length
	        LogicalType::BIGINT,    // tag_as
	        LogicalType::BIGINT,    // tag_xs
	        LogicalType::BIGINT,    // tag_ys
	        LogicalType::BIGINT,    // tag_xn
	        LogicalType::BIGINT,    // tag_xm
	        LogicalType::BIGINT,    // tag_xo
	        LogicalType::BIGINT,    // tag_xg
	        LogicalType::BIGINT,    // tag_nm
	        LogicalType::VARCHAR,   // tag_yt
	        LogicalType::VARCHAR,   // tag_md
	        LogicalType::VARCHAR};  // tag_sa
}

// Parse minimap2 config parameters from named_parameters map
// Set warn_prebuilt_index=true to warn when k/w are specified but will be ignored
inline void ParseMinimap2ConfigParams(const named_parameter_map_t &params, miint::Minimap2Config &config,
                                       bool warn_prebuilt_index = false) {
	auto preset_param = params.find("preset");
	if (preset_param != params.end() && !preset_param->second.IsNull()) {
		config.preset = preset_param->second.ToString();
	}

	auto max_secondary_param = params.find("max_secondary");
	if (max_secondary_param != params.end() && !max_secondary_param->second.IsNull()) {
		config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
	}

	auto k_param = params.find("k");
	if (k_param != params.end() && !k_param->second.IsNull()) {
		config.k = k_param->second.GetValue<int32_t>();
		if (warn_prebuilt_index) {
			Printer::Print("WARNING: Parameter 'k' is ignored when using pre-built index. "
			               "The k-mer size is baked into the index.\n");
		}
	}

	auto w_param = params.find("w");
	if (w_param != params.end() && !w_param->second.IsNull()) {
		config.w = w_param->second.GetValue<int32_t>();
		if (warn_prebuilt_index) {
			Printer::Print("WARNING: Parameter 'w' is ignored when using pre-built index. "
			               "The window size is baked into the index.\n");
		}
	}

	auto eqx_param = params.find("eqx");
	if (eqx_param != params.end() && !eqx_param->second.IsNull()) {
		config.eqx = eqx_param->second.GetValue<bool>();
	}
}

// Parse bowtie2 config parameters from named_parameters map
inline void ParseBowtie2ConfigParams(const named_parameter_map_t &params, miint::Bowtie2Config &config) {
	auto preset_param = params.find("preset");
	if (preset_param != params.end() && !preset_param->second.IsNull()) {
		config.preset = preset_param->second.ToString();
	}

	auto local_param = params.find("local");
	if (local_param != params.end() && !local_param->second.IsNull()) {
		config.local = local_param->second.GetValue<bool>();
	}

	auto threads_param = params.find("threads");
	if (threads_param != params.end() && !threads_param->second.IsNull()) {
		config.threads = threads_param->second.GetValue<int32_t>();
	}

	auto max_secondary_param = params.find("max_secondary");
	if (max_secondary_param != params.end() && !max_secondary_param->second.IsNull()) {
		config.max_secondary = max_secondary_param->second.GetValue<int32_t>();
	}

	auto extra_args_param = params.find("extra_args");
	if (extra_args_param != params.end() && !extra_args_param->second.IsNull()) {
		config.extra_args = extra_args_param->second.ToString();
	}

	auto quiet_param = params.find("quiet");
	if (quiet_param != params.end() && !quiet_param->second.IsNull()) {
		config.quiet = quiet_param->second.GetValue<bool>();
	}
}

// Output SAMRecordBatch to DataChunk using standard alignment schema
// Returns the number of records output
inline idx_t OutputSAMRecordBatch(DataChunk &output, const miint::SAMRecordBatch &batch, idx_t offset, idx_t count) {
	idx_t field_idx = 0;
	SetAlignResultString(output.data[field_idx++], batch.read_ids, offset, count);
	SetAlignResultUInt16(output.data[field_idx++], batch.flags, offset, count);
	SetAlignResultString(output.data[field_idx++], batch.references, offset, count);
	SetAlignResultInt64(output.data[field_idx++], batch.positions, offset, count);
	SetAlignResultInt64(output.data[field_idx++], batch.stop_positions, offset, count);
	SetAlignResultUInt8(output.data[field_idx++], batch.mapqs, offset, count);
	SetAlignResultString(output.data[field_idx++], batch.cigars, offset, count);
	SetAlignResultString(output.data[field_idx++], batch.mate_references, offset, count);
	SetAlignResultInt64(output.data[field_idx++], batch.mate_positions, offset, count);
	SetAlignResultInt64(output.data[field_idx++], batch.template_lengths, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_as_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xs_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_ys_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xn_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xm_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xo_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_xg_values, offset, count);
	SetAlignResultInt64Nullable(output.data[field_idx++], batch.tag_nm_values, offset, count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_yt_values, offset, count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_md_values, offset, count);
	SetAlignResultStringNullable(output.data[field_idx++], batch.tag_sa_values, offset, count);

	output.SetCardinality(count);
	return count;
}

// Filter out unmapped reads from result batch (in-place)
// Removes any records with the unmapped flag (0x4) set
inline void FilterMappedOnly(miint::SAMRecordBatch &batch) {
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

// =============================================================================
// Sharded alignment utilities
// =============================================================================

// Validate that read_to_shard table has required columns (read_id VARCHAR, shard_name VARCHAR)
inline void ValidateReadToShardSchema(ClientContext &context, const std::string &table_name) {
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

// Raw shard count info returned by ReadShardNameCounts
struct ShardNameCount {
	std::string name;
	idx_t count;
};

// Read shard names and counts from read_to_shard table
// Returns pairs of (shard_name, count) ordered by count descending (largest first)
// Throws if any shard_name is NULL or if the table is empty
inline std::vector<ShardNameCount> ReadShardNameCounts(ClientContext &context, const std::string &table_name) {
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

	std::vector<ShardNameCount> shards;
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

			ShardNameCount info;
			info.name = shard_names[shard_idx].GetString();
			info.count = static_cast<idx_t>(counts[count_idx]);
			shards.push_back(std::move(info));
		}
	}

	if (shards.empty()) {
		throw BinderException("read_to_shard table '%s' is empty or has no valid shard names", table_name);
	}

	return shards;
}

} // namespace duckdb
