#include "alignment_slice.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "AlignmentSlicer.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/helper.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

namespace duckdb {

using ColRole = AlignmentSliceTableFunction::ColRole;

const vector<AlignmentSliceTableFunction::ColumnInfo> &AlignmentSliceTableFunction::GetRecognizedColumns() {
	// name, type, required, is_tag, invalidate_on_trim, is_id
	static const vector<ColumnInfo> columns = {
	    {"read_id", LogicalType::VARCHAR, false, false, false, true},
	    {"flags", LogicalType::USMALLINT, false, false, false, false},
	    {"reference", LogicalType::VARCHAR, false, false, false, true},
	    {"position", LogicalType::BIGINT, true, false, false, false},
	    {"stop_position", LogicalType::BIGINT, true, false, false, false},
	    {"mapq", LogicalType::UTINYINT, false, false, false, false},
	    {"cigar", LogicalType::VARCHAR, true, false, false, false},
	    {"mate_reference", LogicalType::VARCHAR, false, false, false, true},
	    {"mate_position", LogicalType::BIGINT, false, false, false, false},
	    {"template_length", LogicalType::BIGINT, false, false, true, false},
	    {"tag_as", LogicalType::BIGINT, false, true, false, false},
	    {"tag_xs", LogicalType::BIGINT, false, true, false, false},
	    {"tag_ys", LogicalType::BIGINT, false, true, false, false},
	    {"tag_xn", LogicalType::BIGINT, false, true, false, false},
	    {"tag_xm", LogicalType::BIGINT, false, true, false, false},
	    {"tag_xo", LogicalType::BIGINT, false, true, false, false},
	    {"tag_xg", LogicalType::BIGINT, false, true, false, false},
	    {"tag_nm", LogicalType::BIGINT, false, true, false, false},
	    {"tag_yt", LogicalType::VARCHAR, false, true, false, false},
	    {"tag_md", LogicalType::VARCHAR, false, true, false, false},
	    {"tag_sa", LogicalType::VARCHAR, false, true, false, false},
	    {"sequence", LogicalType::VARCHAR, false, false, false, false},
	    {"qual", LogicalType::LIST(LogicalType::UTINYINT), false, false, false, false},
	};
	return columns;
}

unique_ptr<FunctionData> AlignmentSliceTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<Data>();

	// Extract positional parameters
	data->table_name = input.inputs[0].GetValue<string>();
	data->region_start = input.inputs[1].GetValue<int64_t>();
	data->region_stop = input.inputs[2].GetValue<int64_t>();

	// Extract named parameters
	data->include_deletions = false;
	if (input.named_parameters.count("include_deletions")) {
		data->include_deletions = input.named_parameters.at("include_deletions").GetValue<bool>();
	}

	// Validate region
	if (data->region_start < 1) {
		throw BinderException("alignment_slice: region start must be >= 1, got %lld", data->region_start);
	}
	if (data->region_start > data->region_stop) {
		throw BinderException("alignment_slice: region start (%lld) must be <= region stop (%lld)", data->region_start,
		                      data->region_stop);
	}

	// Discover input table schema
	auto table_info = GetTableOrViewColumns(context, data->table_name, "Alignment table");
	auto &col_names = table_info.names;

	// Build case-insensitive name-to-index map for input columns
	std::unordered_map<string, idx_t> input_name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		input_name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Match recognized columns against input schema.
	// Track which input columns are present so we can build a targeted SELECT.
	auto &recognized = GetRecognizedColumns();
	// input_col_present[col.name] = index in input table, or -1 if absent
	std::unordered_map<string, int> input_col_present;
	for (const auto &col : recognized) {
		auto found = input_name_to_idx.find(col.name);
		if (found != input_name_to_idx.end()) {
			input_col_present[col.name] = static_cast<int>(found->second);
		} else if (col.required) {
			throw BinderException("alignment_slice: input table '%s' missing required column '%s'", data->table_name,
			                      col.name);
		} else {
			input_col_present[col.name] = -1;
		}
	}

	// Build SELECT query with only the columns we need (avoids pulling unnecessary columns).
	// Track the mapping from recognized column name to index in the SELECT result.
	vector<string> select_cols;
	std::unordered_map<string, int> select_col_idx; // col.name -> index in SELECT result
	for (const auto &col : recognized) {
		if (input_col_present[col.name] >= 0) {
			select_col_idx[col.name] = static_cast<int>(select_cols.size());
			select_cols.push_back(KeywordHelper::WriteOptionallyQuoted(col.name));
		}
	}
	data->select_query = "SELECT " + StringUtil::Join(select_cols, ", ") + " FROM " +
	                     KeywordHelper::WriteOptionallyQuoted(data->table_name);

	// Store slicer-relevant indices into the SELECT result
	data->select_cigar_idx = select_col_idx.at("cigar");
	data->select_pos_idx = select_col_idx.at("position");
	data->select_stop_pos_idx = select_col_idx.at("stop_position");
	data->select_seq_idx = select_col_idx.count("sequence") ? select_col_idx.at("sequence") : -1;
	data->select_qual_idx = select_col_idx.count("qual") ? select_col_idx.at("qual") : -1;
	data->select_ref_idx = select_col_idx.count("reference") ? select_col_idx.at("reference") : -1;

	// Build output schema: recognized columns present in input, in recognized order.
	// Precompute per-column roles for fast dispatch in Execute (no string comparisons).
	for (const auto &col : recognized) {
		if (input_col_present[col.name] >= 0) {
			data->output_names.push_back(col.name);

			// Identifier columns (read_id, reference, mate_reference) preserve the
			// input's storage type (VARCHAR/BIGINT/UUID) instead of coercing to
			// VARCHAR, matching align_minimap2 (query id_type drives read_id;
			// subject id_type drives reference + mate_reference). Without this a
			// BIGINT/UUID id is silently stringified on the pass-through path.
			// mate_reference is pass-through here (no "="/"*" resolution), so
			// mirroring the column's own type is always correct.
			LogicalType out_type = col.type;
			if (col.is_id) {
				const auto &actual_type = table_info.types[static_cast<idx_t>(input_col_present[col.name])];
				if (!IsAllowedIdType(actual_type)) {
					throw BinderException("alignment_slice: '%s' column must be %s, got %s", col.name,
					                      AllowedIdTypeList(), actual_type.ToString());
				}
				out_type = actual_type;
			}
			data->output_types.push_back(out_type);
			data->output_input_idx.push_back(select_col_idx.at(col.name));

			ColRole role;
			if (col.is_tag || col.invalidate_on_trim) {
				role = ColRole::TAG_OR_INVALIDATE;
			} else if (col.name == "cigar") {
				role = ColRole::CIGAR;
			} else if (col.name == "position") {
				role = ColRole::POSITION;
			} else if (col.name == "stop_position") {
				role = ColRole::STOP_POSITION;
			} else if (col.name == "sequence") {
				role = ColRole::SEQUENCE;
			} else if (col.name == "qual") {
				role = ColRole::QUAL;
			} else {
				role = ColRole::PASS_THROUGH;
			}
			data->output_roles.push_back(role);
		}
	}

	names = data->output_names;
	return_types = data->output_types;

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignmentSliceTableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Create the slicer
	gstate->slicer = make_uniq<miint::AlignmentSlicer>(data.region_start, data.region_stop, data.include_deletions);

	// Separate connection avoids deadlocking the current context.
	// SendQuery streams chunks lazily — only one chunk in memory at a time.
	auto &db = DatabaseInstance::GetDatabase(context);
	gstate->conn = make_uniq<Connection>(db);
	InheritTempObjects(context, *gstate->conn);
	gstate->query_result = gstate->conn->SendQuery(data.select_query);

	return gstate;
}

void AlignmentSliceTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	idx_t out_row = 0;
	idx_t num_out_cols = data.output_names.size();

	while (out_row < STANDARD_VECTOR_SIZE) {
		// Fetch next chunk from stream if needed
		if (!gstate.current_chunk || gstate.chunk_offset >= gstate.current_chunk->size()) {
			if (gstate.stream_exhausted) {
				break;
			}
			gstate.current_chunk = gstate.query_result->Fetch();
			gstate.chunk_offset = 0;
			if (!gstate.current_chunk || gstate.current_chunk->size() == 0) {
				gstate.stream_exhausted = true;
				break;
			}
		}

		auto &chunk = *gstate.current_chunk;
		idx_t i = gstate.chunk_offset;
		gstate.chunk_offset++;

		// Single-reference check (before NULL skip — all rows must be same reference)
		if (data.select_ref_idx >= 0) {
			auto ref_val = chunk.data[data.select_ref_idx].GetValue(i);
			if (!ref_val.IsNull()) {
				string ref_name = ref_val.GetValue<string>();
				if (!gstate.has_seen_reference) {
					gstate.seen_reference = ref_name;
					gstate.has_seen_reference = true;
				} else if (ref_name != gstate.seen_reference) {
					throw InvalidInputException("alignment_slice: multiple references found in input ('%s' and '%s'). "
					                            "Filter to a single reference before slicing.",
					                            gstate.seen_reference, ref_name);
				}
			}
		}

		// Extract required fields, skip rows with NULLs
		auto cigar_val = chunk.data[data.select_cigar_idx].GetValue(i);
		auto pos_val = chunk.data[data.select_pos_idx].GetValue(i);
		auto stop_pos_val = chunk.data[data.select_stop_pos_idx].GetValue(i);

		if (cigar_val.IsNull() || pos_val.IsNull() || stop_pos_val.IsNull()) {
			continue;
		}

		string cigar = cigar_val.GetValue<string>();
		int64_t position = pos_val.GetValue<int64_t>();
		int64_t stop_position = stop_pos_val.GetValue<int64_t>();

		// Build SliceInput
		miint::SliceInput slice_input;
		slice_input.cigar = cigar;
		slice_input.position = position;
		slice_input.stop_position = stop_position;

		if (data.select_seq_idx >= 0) {
			auto seq_val = chunk.data[data.select_seq_idx].GetValue(i);
			if (!seq_val.IsNull()) {
				slice_input.sequence = seq_val.GetValue<string>();
			}
		}

		if (data.select_qual_idx >= 0) {
			auto qual_val = chunk.data[data.select_qual_idx].GetValue(i);
			if (!qual_val.IsNull()) {
				auto qual_list = ListValue::GetChildren(qual_val);
				for (const auto &q : qual_list) {
					slice_input.quality.push_back(q.GetValue<uint8_t>());
				}
			}
		}

		// Slice
		auto result = gstate.slicer->Slice(slice_input);
		if (!result.overlaps) {
			continue;
		}

		// Determine if trimming occurred by comparing CIGAR (the canonical check —
		// position comparison can miss cases where only CIGAR changes, e.g. soft clips)
		bool was_trimmed = (result.cigar != cigar);

		// Build output row using precomputed roles (no string comparisons)
		for (idx_t col = 0; col < num_out_cols; col++) {
			int src_idx = data.output_input_idx[col];

			switch (data.output_roles[col]) {
			case ColRole::TAG_OR_INVALIDATE:
				if (was_trimmed) {
					FlatVector::SetNull(output.data[col], out_row, true);
				} else {
					output.data[col].SetValue(out_row, chunk.data[src_idx].GetValue(i));
				}
				break;
			case ColRole::CIGAR:
				output.data[col].SetValue(out_row, Value(result.cigar));
				break;
			case ColRole::POSITION:
				output.data[col].SetValue(out_row, Value::BIGINT(result.position));
				break;
			case ColRole::STOP_POSITION:
				output.data[col].SetValue(out_row, Value::BIGINT(result.stop_position));
				break;
			case ColRole::SEQUENCE:
				if (!result.sequence.empty()) {
					output.data[col].SetValue(out_row, Value(result.sequence));
				} else if (!slice_input.sequence.empty()) {
					// Slicer returned empty from non-empty input (length mismatch) — NULL
					FlatVector::SetNull(output.data[col], out_row, true);
				} else {
					output.data[col].SetValue(out_row, chunk.data[src_idx].GetValue(i));
				}
				break;
			case ColRole::QUAL:
				if (!result.quality.empty()) {
					vector<Value> qual_vals;
					for (auto q : result.quality) {
						qual_vals.push_back(Value::UTINYINT(q));
					}
					output.data[col].SetValue(out_row, Value::LIST(LogicalType::UTINYINT, std::move(qual_vals)));
				} else if (!slice_input.quality.empty()) {
					FlatVector::SetNull(output.data[col], out_row, true);
				} else {
					output.data[col].SetValue(out_row, chunk.data[src_idx].GetValue(i));
				}
				break;
			case ColRole::PASS_THROUGH:
				output.data[col].SetValue(out_row, chunk.data[src_idx].GetValue(i));
				break;
			}
		}

		out_row++;
	}

	output.SetCardinality(out_row);
}

void AlignmentSliceTableFunction::Register(ExtensionLoader &loader) {
	auto tf = TableFunction("alignment_slice", {LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::BIGINT},
	                        Execute, Bind, InitGlobal);
	tf.named_parameters["include_deletions"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(tf);
}

} // namespace duckdb
