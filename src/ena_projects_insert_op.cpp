// SPDX-License-Identifier: MIT
//
// Project-specific data mapping for ENAProjectsInsert. The PhysicalOperator
// scaffolding (Sink, Source, Finalize, GlobalSinkState, GlobalSourceState)
// lives in ENAObjectInsertOperator (ena_object_insert_op.hpp). Keep this file
// focused on the project-row → ProjectSpec parse and the RETURNING projection.

#include "ena_projects_insert_op.hpp"

#include "ena_insert_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/date.hpp"
#include "duckdb/common/types/value.hpp"

namespace duckdb {

namespace {

constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_TITLE = 1;
constexpr idx_t COL_DESCRIPTION = 2;
constexpr idx_t COL_PROJECT_TYPE = 3;
constexpr idx_t COL_UMBRELLA_CHILDREN = 4;
constexpr idx_t COL_ATTRIBUTES = 5;
constexpr idx_t COL_PRJEB_ACCESSION = 6;
constexpr idx_t COL_ERP_ACCESSION = 7;
constexpr idx_t COL_HOLD_UNTIL_DATE = 8;
constexpr idx_t TABLE_COLUMN_COUNT = 9;

} // namespace

ENABuiltSpecs<miint::ProjectSpec>
ENAProjectsInsert::BuildFromBuffer(ColumnDataCollection &buffer,
                                   const physical_index_vector_t<idx_t> &column_index_map) {
	ENABuiltSpecs<miint::ProjectSpec> out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.projects requires the 'alias' column");
	}
	// `column_index_map` is per-statement, not per-row; resolve every optional
	// column index up front so the row loop doesn't repeat the work.
	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
	const idx_t description_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESCRIPTION);
	const idx_t project_type_idx = ResolveInputColumn(column_index_map, input_columns, COL_PROJECT_TYPE);
	const idx_t hold_idx = ResolveInputColumn(column_index_map, input_columns, COL_HOLD_UNTIL_DATE);
	const idx_t umbrella_idx = ResolveInputColumn(column_index_map, input_columns, COL_UMBRELLA_CHILDREN);
	const idx_t attrs_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTES);

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::ProjectSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.projects: 'alias' must be non-empty");
			}
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			if (description_idx != DConstants::INVALID_INDEX) {
				spec.description = ValueToVarchar(chunk.data[description_idx].GetValue(row));
			}
			if (project_type_idx != DConstants::INVALID_INDEX) {
				spec.project_type = ValueToVarchar(chunk.data[project_type_idx].GetValue(row));
			}
			if (hold_idx != DConstants::INVALID_INDEX) {
				const auto row_hold = ValueToDateString(chunk.data[hold_idx].GetValue(row));
				if (!row_hold.empty()) {
					if (out.hold_until_date.empty()) {
						out.hold_until_date = row_hold;
					} else if (out.hold_until_date != row_hold) {
						throw InvalidInputException(
						    "INSERT INTO ena.projects: per-row hold_until_date values must agree across the "
						    "statement (got '%s' and '%s')",
						    out.hold_until_date, row_hold);
					}
				}
			}
			// umbrella_children / attributes / prjeb_accession / erp_accession: not
			// yet wired in Phase 4. Reject explicit non-null umbrella_children so
			// users don't think it's silently honoured.
			if (umbrella_idx != DConstants::INVALID_INDEX) {
				const auto v = chunk.data[umbrella_idx].GetValue(row);
				if (!v.IsNull() && !ListValue::GetChildren(v).empty()) {
					throw NotImplementedException(
					    "INSERT INTO ena.projects: 'umbrella_children' is not supported in this build");
				}
			}
			if (attrs_idx != DConstants::INVALID_INDEX) {
				const auto v = chunk.data[attrs_idx].GetValue(row);
				if (!v.IsNull() && !ListValue::GetChildren(v).empty()) {
					throw NotImplementedException(
					    "INSERT INTO ena.projects: 'attributes' is not supported in this build");
				}
			}
			out.specs.push_back(std::move(spec));
		}
	}
	return out;
}

void ENAProjectsInsert::AppendReturningRows(ColumnDataCollection &return_collection,
                                            const vector<LogicalType> &return_types,
                                            const std::vector<miint::ProjectSpec> &specs,
                                            const miint::ENASubmissionOutcome &outcome) {
	// `LogicalInsert::ResolveTypes` sets `op.types = table.GetTypes()` when
	// return_chunk is true, so we always receive the full 9-column projects
	// schema. The DuckDB executor inserts a PhysicalProjection above us to
	// trim down to the user-requested RETURNING columns. Hard-asserting keeps
	// us honest if that contract ever changes upstream.
	D_ASSERT(return_types.size() == TABLE_COLUMN_COUNT);
	// `outcome.rows` is built in alias-order from `specs`, so positional
	// indexing is sufficient and avoids an O(n) per-row lookup.
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
		chunk.data[COL_PROJECT_TYPE].SetValue(idx, Value(spec.project_type));
		chunk.data[COL_UMBRELLA_CHILDREN].SetValue(idx, Value(LogicalType::LIST(LogicalType::VARCHAR)));
		chunk.data[COL_ATTRIBUTES].SetValue(idx, Value(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)));
		chunk.data[COL_PRJEB_ACCESSION].SetValue(idx, Value(row.prjeb_accession));
		chunk.data[COL_ERP_ACCESSION].SetValue(idx, row.erp_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                      : Value(row.erp_accession));
		if (row.hold_until_date.empty()) {
			chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value(LogicalType::DATE));
		} else {
			date_t d;
			idx_t pos = 0;
			bool special = false;
			if (Date::TryConvertDate(row.hold_until_date.c_str(), row.hold_until_date.size(), pos, d, special) ==
			    DateCastResult::SUCCESS) {
				chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value::DATE(d));
			} else {
				chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value(LogicalType::DATE));
			}
		}
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
