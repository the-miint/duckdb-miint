#pragma once
/*
 * Shared utilities for minimap2 alignment table functions (align_minimap2, align_minimap2_sharded)
 * Provides common config parsing, output schema, and result output logic.
 */

#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "align_result_utils.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"

namespace duckdb {

// Batch size for reading queries (shared across minimap2 functions)
static constexpr idx_t MINIMAP2_QUERY_BATCH_SIZE = 1024;

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

} // namespace duckdb
