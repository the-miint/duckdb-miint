// SPDX-License-Identifier: MIT
//
// Sample-specific data mapping for ENASamplesInsert. The PhysicalOperator
// scaffolding (Sink, Source, Finalize, GlobalSinkState, GlobalSourceState)
// lives in ENAObjectInsertOperator (ena_object_insert_op.hpp). Keep this file
// focused on the sample-row → SampleSpec parse and the RETURNING projection.

#include "ena_samples_insert_op.hpp"

#include "ena_insert_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"

namespace duckdb {

namespace {

constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_TITLE = 1;
constexpr idx_t COL_DESCRIPTION = 2;
constexpr idx_t COL_TAXON_ID = 3;
constexpr idx_t COL_SCIENTIFIC_NAME = 4;
constexpr idx_t COL_CHECKLIST = 5;
constexpr idx_t COL_ATTRIBUTES = 6;
constexpr idx_t COL_ATTRIBUTE_UNITS = 7;
constexpr idx_t COL_ERS_ACCESSION = 8;
constexpr idx_t COL_SAMEA_ACCESSION = 9;
constexpr idx_t TABLE_COLUMN_COUNT = 10;

// Pull a MAP(VARCHAR,VARCHAR) value into a (tag,value) vector. NULL maps and
// empty maps both yield no entries. Iteration order follows DuckDB's
// map-storage order, which matches the caller's MAP-literal order in practice
// — the same guarantee the envelope builder already documents for project
// attributes. `column_label` shows up in error messages so users can tell
// which MAP column the bad entry came from (`attributes` vs `attribute_units`).
std::vector<std::pair<std::string, std::string>> ExtractKeyValueMap(const Value &v, const char *column_label) {
	std::vector<std::pair<std::string, std::string>> out;
	if (v.IsNull()) {
		return out;
	}
	const auto &entries = ListValue::GetChildren(v);
	out.reserve(entries.size());
	for (const auto &entry : entries) {
		const auto &kv_children = StructValue::GetChildren(entry);
		if (kv_children.size() != 2) {
			throw InvalidInputException("INSERT INTO ena.samples: '%s' MAP entries must have key+value pairs",
			                            column_label);
		}
		const auto key = ValueToVarchar(kv_children[0]);
		const auto value = ValueToVarchar(kv_children[1]);
		if (key.empty()) {
			throw InvalidInputException("INSERT INTO ena.samples: '%s' map keys must be non-empty", column_label);
		}
		out.emplace_back(key, value);
	}
	return out;
}

} // namespace

ENABuiltSpecs<miint::SampleSpec>
ENASamplesInsert::BuildFromBuffer(ColumnDataCollection &buffer,
                                  const physical_index_vector_t<idx_t> &column_index_map) {
	ENABuiltSpecs<miint::SampleSpec> out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.samples requires the 'alias' column");
	}
	const idx_t taxon_idx = ResolveInputColumn(column_index_map, input_columns, COL_TAXON_ID);
	if (taxon_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.samples requires the 'taxon_id' column");
	}
	// `column_index_map` is per-statement, not per-row; resolve every optional
	// column index up front so the row loop doesn't repeat the work.
	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
	const idx_t description_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESCRIPTION);
	const idx_t scientific_idx = ResolveInputColumn(column_index_map, input_columns, COL_SCIENTIFIC_NAME);
	const idx_t checklist_idx = ResolveInputColumn(column_index_map, input_columns, COL_CHECKLIST);
	const idx_t attrs_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTES);
	const idx_t units_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTE_UNITS);

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::SampleSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.samples: 'alias' must be non-empty");
			}
			const auto taxon_val = chunk.data[taxon_idx].GetValue(row);
			if (taxon_val.IsNull()) {
				throw InvalidInputException("INSERT INTO ena.samples: 'taxon_id' must be non-null for alias '%s'",
				                            spec.alias);
			}
			spec.taxon_id = taxon_val.GetValue<int64_t>();
			if (spec.taxon_id <= 0) {
				throw InvalidInputException("INSERT INTO ena.samples: 'taxon_id' must be > 0 for alias '%s' (got %lld)",
				                            spec.alias, static_cast<long long>(spec.taxon_id));
			}
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			if (description_idx != DConstants::INVALID_INDEX) {
				spec.description = ValueToVarchar(chunk.data[description_idx].GetValue(row));
			}
			if (scientific_idx != DConstants::INVALID_INDEX) {
				spec.scientific_name = ValueToVarchar(chunk.data[scientific_idx].GetValue(row));
			}
			if (checklist_idx != DConstants::INVALID_INDEX) {
				spec.checklist = ValueToVarchar(chunk.data[checklist_idx].GetValue(row));
			}
			if (attrs_idx != DConstants::INVALID_INDEX) {
				spec.attributes = ExtractKeyValueMap(chunk.data[attrs_idx].GetValue(row), "attributes");
			}
			if (units_idx != DConstants::INVALID_INDEX) {
				spec.attribute_units = ExtractKeyValueMap(chunk.data[units_idx].GetValue(row), "attribute_units");
			}
			// HOLD on samples is not exposed in Phase 5 (no hold_until_date
			// column on ena.samples).
			out.specs.push_back(std::move(spec));
		}
	}
	return out;
}

void ENASamplesInsert::AppendReturningRows(ColumnDataCollection &return_collection,
                                           const vector<LogicalType> &return_types,
                                           const std::vector<miint::SampleSpec> &specs,
                                           const miint::ENASamplesSubmissionOutcome &outcome) {
	D_ASSERT(return_types.size() == TABLE_COLUMN_COUNT);
	D_ASSERT(outcome.rows.size() == specs.size());

	DataChunk chunk;
	chunk.Initialize(Allocator::DefaultAllocator(), return_types);
	for (idx_t i = 0; i < outcome.rows.size(); i++) {
		const auto &row = outcome.rows[i];
		const auto &spec = specs[i];
		const auto idx = chunk.size();
		chunk.data[COL_ALIAS].SetValue(idx, Value(row.alias));
		chunk.data[COL_TITLE].SetValue(idx, Value(spec.title));
		chunk.data[COL_DESCRIPTION].SetValue(idx, Value(spec.description));
		chunk.data[COL_TAXON_ID].SetValue(idx, Value::INTEGER(NumericCast<int32_t>(spec.taxon_id)));
		chunk.data[COL_SCIENTIFIC_NAME].SetValue(idx, Value(spec.scientific_name));
		chunk.data[COL_CHECKLIST].SetValue(idx, Value(spec.checklist));
		// attributes column emits NULL on RETURNING in Phase 5 — the user-
		// supplied attribute list is preserved verbatim in
		// `ena.submission_log.request_payload`, and the server's echo (which
		// may include system-injected attributes the user did not write) is
		// in `submission_log.receipt`. Same trade-off as ena.projects in
		// Phase 4. Real MAP-value emission is a future task.
		chunk.data[COL_ATTRIBUTES].SetValue(idx, Value(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)));
		chunk.data[COL_ATTRIBUTE_UNITS].SetValue(idx,
		                                         Value(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)));
		chunk.data[COL_ERS_ACCESSION].SetValue(idx, Value(row.ers_accession));
		chunk.data[COL_SAMEA_ACCESSION].SetValue(idx, row.samea_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                          : Value(row.samea_accession));
		chunk.SetCardinality(idx + 1);
		if (chunk.size() == STANDARD_VECTOR_SIZE) {
			return_collection.Append(chunk);
			chunk.Reset();
		}
	}
	if (chunk.size() > 0) {
		return_collection.Append(chunk);
	}
}

} // namespace duckdb
