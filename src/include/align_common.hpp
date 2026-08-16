#pragma once
/*
 * Shared utilities for the minimap2 alignment table functions
 * (align_minimap2, align_minimap2_sharded). The bowtie2 table functions
 * route through the gpl-boundary daemon and live in
 * `align_bowtie2_daemon_common.{hpp,cpp}`.
 */

#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "align_result_utils.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
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
#include <cstdlib>
#include <string>

namespace duckdb {

// Batch size for reading queries (shared across alignment functions)
static constexpr idx_t ALIGNMENT_QUERY_BATCH_SIZE = 1024;
// Backward compatibility alias
static constexpr idx_t MINIMAP2_QUERY_BATCH_SIZE = ALIGNMENT_QUERY_BATCH_SIZE;

// Get the standard alignment output column names
inline std::vector<std::string> GetAlignmentOutputNames() {
	return {"read_id",        "flags",         "reference",       "position", "stop_position", "mapq",   "cigar",
	        "mate_reference", "mate_position", "template_length", "tag_as",   "tag_xs",        "tag_ys", "tag_xn",
	        "tag_xm",         "tag_xo",        "tag_xg",          "tag_nm",   "tag_yt",        "tag_md", "tag_sa"};
}

// Get the standard alignment output column types.
// `query_id_type` drives `read_id`; `subject_id_type` drives `reference` and
// `mate_reference`. Both must be VARCHAR, BIGINT, or UUID. No default arguments —
// every caller must make the choice explicit so the audit can't slip.
inline std::vector<LogicalType> GetAlignmentOutputTypes(const LogicalType &query_id_type,
                                                        const LogicalType &subject_id_type) {
	return {query_id_type,          // read_id (mirrors query)
	        LogicalType::USMALLINT, // flags
	        subject_id_type,        // reference (mirrors subject)
	        LogicalType::BIGINT,    // position
	        LogicalType::BIGINT,    // stop_position
	        LogicalType::UTINYINT,  // mapq
	        LogicalType::VARCHAR,   // cigar
	        subject_id_type,        // mate_reference (mirrors subject)
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

// Parse minimap2's -f spec (high-occurrence minimizer filter) into config.
//
// Deliberately a transcription of minimap2's own CLI parsing (ext/minimap2/main.c, case 'f'):
//
//     x = strtod(arg, &p);
//     if (x < 1.0) opt.mid_occ_frac = x, opt.mid_occ = 0;
//     else         opt.mid_occ = (int)(x + .499);
//     if (*p == ',') opt.max_occ = (int)(strtod(p+1, &p) + .499);
//
// The `x < 1.0` split is load-bearing and easy to get wrong. A value below 1 is a FRACTION: it
// sets mid_occ_frac and zeroes mid_occ so that mm_mapopt_update derives the threshold from the
// index's own minimizer distribution. So `occ_filter := 0` does not mean "threshold of zero" and
// must not be written to mid_occ directly -- it means "filter the top 0 fraction", and
// mm_idx_cal_max_occ returns INT32_MAX for f <= 0, which mm_mapopt_update then clamps to
// max_mid_occ (1000000 under the default `sr` preset). That is minimap2's way of disabling the
// filter. Note that presets which narrow max_mid_occ to 500 (lr:hq, map-hifi, map-ccs, lr:hqae,
// map-iclr) therefore clamp `occ_filter := 0` to 500 rather than disabling it -- inherited from
// minimap2, not introduced here.
inline void ParseOccFilterSpec(const std::string &spec, miint::Minimap2Config &config) {
	// Split on ',' BEFORE parsing, rather than letting strtod stop at it. strtod honours
	// LC_NUMERIC, so in a locale where ',' is the decimal separator it would read "1000,5000" as
	// the single value 1000.5 and consume the comma -- silently discarding the second value instead
	// of failing. Splitting first makes the two-value form locale-independent.
	const auto comma = spec.find(',');
	const std::string first_spec = spec.substr(0, comma);
	const bool has_second = (comma != std::string::npos);

	// Rejects NaN and infinity as well as negatives and anything past int32: the comparison is
	// written as !(in range) so that NaN, for which every ordered comparison is false, fails here.
	// This matters more than it looks -- `occ_filter := 1e12` is a plausible way to ask for "off",
	// and casting it to int32_t is undefined (x86 yields INT32_MIN), which mm_mapopt_update would
	// then read as mid_occ <= 0 and quietly re-derive the DEFAULT filter: the exact silent masking
	// #187 exists to remove.
	auto parse_value = [&spec](const std::string &text, const char *what) {
		const char *begin = text.c_str();
		char *end = nullptr;
		double v = std::strtod(begin, &end);
		if (end == begin || *end != '\0') {
			throw InvalidInputException("%s must be a number or 'INT1,INT2' (minimap2's -f), got '%s'", what, spec);
		}
		if (!(v >= 0.0 && v <= 2147483647.0)) {
			throw InvalidInputException("%s must be between 0 and 2147483647, got '%s'", what, spec);
		}
		return v;
	};

	const double x = parse_value(first_spec, "occ_filter");
	if (x < 1.0) {
		config.occ_mid_frac = static_cast<float>(x);
		config.occ_mid = 0;
	} else {
		config.occ_mid = static_cast<int32_t>(x + .499);
	}

	if (has_second) {
		// Assigned unconditionally, including 0, because main.c:331 does: `-f 100,0` means max_occ=0
		// (no re-chain pass), which is distinguishable from "no second value given" only by the
		// presence of the comma. Hence the -1 "unset" sentinel on occ_max rather than 0.
		const double y = parse_value(spec.substr(comma + 1), "occ_filter second value");
		config.occ_max = static_cast<int32_t>(y + .499);
	}

	config.occ_filter_set = true;
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

	auto min_qcov_param = params.find("min_chain_coverage");
	if (min_qcov_param != params.end() && !min_qcov_param->second.IsNull()) {
		config.min_chain_coverage = min_qcov_param->second.GetValue<float>();
		if (config.min_chain_coverage < 0.0f || config.min_chain_coverage > 1.0f) {
			throw InvalidInputException("min_chain_coverage must be between 0.0 and 1.0");
		}
	}

	auto occ_param = params.find("occ_filter");
	if (occ_param != params.end() && !occ_param->second.IsNull()) {
		ParseOccFilterSpec(occ_param->second.ToString(), config);
	}

	auto include_unmapped_param = params.find("include_unmapped");
	if (include_unmapped_param != params.end() && !include_unmapped_param->second.IsNull()) {
		config.include_unmapped = include_unmapped_param->second.GetValue<bool>();
	}
}

// Output SAMRecordBatch to DataChunk using the standard alignment schema.
// `query_id_type` and `subject_id_type` drive the three id columns
// (read_id / reference / mate_reference). Both must be VARCHAR, BIGINT, or UUID.
// For non-VARCHAR subject output, the SAM "=" sentinel in mate_reference is
// resolved to the row's `reference` value before emit (the literal "=" has no
// BIGINT/UUID encoding). VARCHAR output preserves "=" as-is, matching
// pre-existing user-observable behavior.
// Returns the number of records output.
inline idx_t OutputSAMRecordBatch(DataChunk &output, const miint::SAMRecordBatch &batch, idx_t offset, idx_t count,
                                  const LogicalType &query_id_type, const LogicalType &subject_id_type) {
	idx_t field_idx = 0;
	EmitIdColumnFromStrings(output.data[field_idx++], batch.read_ids, offset, count, query_id_type);
	SetAlignResultUInt16(output.data[field_idx++], batch.flags, offset, count);
	EmitIdColumnFromStrings(output.data[field_idx++], batch.references, offset, count, subject_id_type);
	SetAlignResultInt64(output.data[field_idx++], batch.positions, offset, count);
	SetAlignResultInt64(output.data[field_idx++], batch.stop_positions, offset, count);
	SetAlignResultUInt8(output.data[field_idx++], batch.mapqs, offset, count);
	SetAlignResultString(output.data[field_idx++], batch.cigars, offset, count);

	// Emit mate_reference. For non-VARCHAR subjects, resolve any "=" sentinel
	// (mate maps to same reference as primary) into the row's reference value
	// via a local copy — the codec rejects "=" for BIGINT/UUID and there is no
	// in-band way to encode it.
	if (subject_id_type.id() != LogicalTypeId::VARCHAR) {
		std::vector<std::string> resolved_mate_refs;
		resolved_mate_refs.reserve(count);
		for (idx_t j = 0; j < count; j++) {
			const auto &mr = batch.mate_references[offset + j];
			resolved_mate_refs.push_back(mr == "=" ? batch.references[offset + j] : mr);
		}
		EmitIdColumnFromStrings(output.data[field_idx++], resolved_mate_refs, 0, count, subject_id_type);
	} else {
		EmitIdColumnFromStrings(output.data[field_idx++], batch.mate_references, offset, count, subject_id_type);
	}

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

	// An unmapped row (flag 0x4) has no reference, coordinates, MAPQ or CIGAR to report, so emit
	// SQL NULL rather than the SAM text sentinels: `WHERE reference IS NULL` is then the test for
	// "measured, did not align", which is the whole point of include_unmapped (#185). The tag
	// columns already null themselves via their -1 / "" sentinels.
	//
	// This loop is inert unless include_unmapped is on. Only align_minimap2 and
	// align_minimap2_sharded call this function, and neither ever emitted a row with 0x4 set --
	// both skip every reg with rid < 0 -- so no pre-existing output changes.
	//
	// Must run after the emitters: EmitIdCell and the nullable tag setters write validity per row
	// and would overwrite these. The plain int64/uint8/string setters do NOT touch validity at all,
	// so for those columns what actually keeps stale NULLs from leaking across chunks is
	// DataChunk::Reset() clearing validity before each GetData -- a dependency that was irrelevant
	// before this loop existed, since nothing ever wrote a NULL into them.
	static constexpr miint::SAMRecordField UNMAPPED_NULL_FIELDS[] = {
	    miint::SAMRecordField::REFERENCE, miint::SAMRecordField::POSITION, miint::SAMRecordField::STOP_POSITION,
	    miint::SAMRecordField::MAPQ, miint::SAMRecordField::CIGAR};
	for (idx_t j = 0; j < count; j++) {
		if ((batch.flags[offset + j] & 0x4) == 0) {
			continue;
		}
		for (auto field : UNMAPPED_NULL_FIELDS) {
			FlatVector::SetNull(output.data[static_cast<idx_t>(field)], j, true);
		}
	}

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

// Validate that the read_to_shard table has the required columns.
//   - `shard_name` is always VARCHAR.
//   - `read_id` defaults to VARCHAR (back-compat). When `expected_read_id_type`
//     is non-INVALID, the column must match it exactly — supports BIGINT once
//     the caller has captured the query table's id type. The strict equality
//     check keeps the downstream shard join (BuildShardedQueryReadsSelect for
//     minimap2, OpenCurrentShardStream for bowtie2) well-typed without relying on
//     implicit casts.
inline void ValidateReadToShardSchema(ClientContext &context, const std::string &table_name,
                                      const LogicalType &expected_read_id_type = LogicalType(LogicalTypeId::INVALID)) {
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
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		col_names = col_info->names;
		col_types = col_info->types;
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
	const auto &actual_id_type = col_types[it_read_id->second];
	if (expected_read_id_type.id() == LogicalTypeId::INVALID) {
		if (actual_id_type.id() != LogicalTypeId::VARCHAR) {
			throw BinderException("Column 'read_id' in read_to_shard table '%s' must be VARCHAR", table_name);
		}
	} else {
		if (!IsAllowedIdType(actual_id_type)) {
			throw BinderException("Column 'read_id' in read_to_shard table '%s' must be %s", table_name,
			                      AllowedIdTypeList());
		}
		if (actual_id_type.id() != expected_read_id_type.id()) {
			throw BinderException("Column 'read_id' in read_to_shard table '%s' is %s but must match the query "
			                      "table's read_id type (%s)",
			                      table_name, actual_id_type.ToString(), expected_read_id_type.ToString());
		}
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
//
// `join_query_table`, when non-empty, restricts the count to reads that actually
// exist in that query relation (`read_to_shard JOIN <query_table> USING read_id`)
// instead of counting mapping rows. Callers that use the result to verify how many
// reads a shard *should* deliver need this form: a mapping legitimately listing
// reads absent from the query relation would otherwise look like data loss (#229).
inline std::vector<ShardNameCount> ReadShardNameCounts(ClientContext &context, const std::string &table_name,
                                                       const std::string &join_query_table = "") {
	auto conn = MakeReadOnlyHelperConnection(context);

	// Query shard counts ordered by count descending (largest first)
	std::string from = KeywordHelper::WriteOptionallyQuoted(table_name) + " rts";
	if (!join_query_table.empty()) {
		from += " JOIN " + KeywordHelper::WriteOptionallyQuoted(join_query_table) + " q ON q.read_id = rts.read_id";
	}
	std::string query = "SELECT rts.shard_name AS shard_name, COUNT(*) as cnt FROM " + from +
	                    " GROUP BY rts.shard_name ORDER BY cnt DESC";

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read shard counts from '%s': %s", table_name, query_result->GetError());
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
