// SPDX-License-Identifier: MIT
//
// Experiment-specific data mapping for ENAExperimentsInsert. The
// PhysicalOperator scaffolding (Sink, Source, Finalize, GlobalSinkState,
// GlobalSourceState) lives in ENAObjectInsertOperator
// (ena_object_insert_op.hpp). Keep this file focused on the
// experiment-row → ExperimentSpec parse and the RETURNING projection.

#include "ena_experiments_insert_op.hpp"

#include "ena_insert_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"

#include <cctype>

namespace duckdb {

namespace {

// Column indices match `BuildENATableInfo(EXPERIMENTS)` in src/ena_storage.cpp.
// Keep both in sync — the order of columns is the contract.
constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_TITLE = 1;
constexpr idx_t COL_STUDY_REF = 2;
constexpr idx_t COL_SAMPLE_DESCRIPTOR = 3;
constexpr idx_t COL_DESIGN_DESCRIPTION = 4;
constexpr idx_t COL_LIBRARY_NAME = 5;
constexpr idx_t COL_LIBRARY_STRATEGY = 6;
constexpr idx_t COL_LIBRARY_SOURCE = 7;
constexpr idx_t COL_LIBRARY_SELECTION = 8;
constexpr idx_t COL_LIBRARY_LAYOUT = 9;
constexpr idx_t COL_PLATFORM = 10;
constexpr idx_t COL_INSTRUMENT_MODEL = 11;
constexpr idx_t COL_ERX_ACCESSION = 12;
constexpr idx_t TABLE_COLUMN_COUNT = 13;

// study_ref / sample_descriptor strings whose shape matches an ENA accession
// (one of the canonical prefixes followed by digits only) route to
// RefDescriptor::accession; everything else is treated as refname (the
// parent's alias). The "digits only" suffix check matters: a user alias like
// `ERPmycoolstudy` would otherwise be silently misclassified as an accession,
// which the server would then fail to find. Real accessions are always
// `<PREFIX><NUMERIC>` (see localdocs/ena-research-webin-v2-deep.md §1.1).
miint::RefDescriptor ToRef(const std::string &value, std::initializer_list<const char *> accession_prefixes) {
	miint::RefDescriptor ref;
	for (const auto *p : accession_prefixes) {
		const std::string prefix(p);
		if (value.size() <= prefix.size() || value.compare(0, prefix.size(), prefix) != 0) {
			continue;
		}
		bool all_digits = true;
		for (size_t i = prefix.size(); i < value.size(); ++i) {
			if (!std::isdigit(static_cast<unsigned char>(value[i]))) {
				all_digits = false;
				break;
			}
		}
		if (all_digits) {
			ref.accession = value;
			return ref;
		}
	}
	ref.refname = value;
	return ref;
}

miint::ENALibraryLayout ParseLayout(const std::string &alias, const std::string &v) {
	if (StringUtil::CIEquals(v, "PAIRED")) {
		return miint::ENALibraryLayout::PAIRED;
	}
	if (StringUtil::CIEquals(v, "SINGLE")) {
		return miint::ENALibraryLayout::SINGLE;
	}
	throw InvalidInputException("INSERT INTO ena.experiments: 'library_layout' must be 'SINGLE' or 'PAIRED' for "
	                            "alias '%s' (got '%s')",
	                            alias, v);
}

} // namespace

ENABuiltSpecs<miint::ExperimentSpec>
ENAExperimentsInsert::BuildFromBuffer(ColumnDataCollection &buffer,
                                      const physical_index_vector_t<idx_t> &column_index_map) {
	ENABuiltSpecs<miint::ExperimentSpec> out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.experiments requires the 'alias' column");
	}

	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
	const idx_t study_idx = ResolveInputColumn(column_index_map, input_columns, COL_STUDY_REF);
	const idx_t sample_idx = ResolveInputColumn(column_index_map, input_columns, COL_SAMPLE_DESCRIPTOR);
	const idx_t design_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESIGN_DESCRIPTION);
	const idx_t libname_idx = ResolveInputColumn(column_index_map, input_columns, COL_LIBRARY_NAME);
	const idx_t strategy_idx = ResolveInputColumn(column_index_map, input_columns, COL_LIBRARY_STRATEGY);
	const idx_t source_idx = ResolveInputColumn(column_index_map, input_columns, COL_LIBRARY_SOURCE);
	const idx_t selection_idx = ResolveInputColumn(column_index_map, input_columns, COL_LIBRARY_SELECTION);
	const idx_t layout_idx = ResolveInputColumn(column_index_map, input_columns, COL_LIBRARY_LAYOUT);
	const idx_t platform_idx = ResolveInputColumn(column_index_map, input_columns, COL_PLATFORM);
	const idx_t instrument_idx = ResolveInputColumn(column_index_map, input_columns, COL_INSTRUMENT_MODEL);

	// Per-statement required-column checks — fail fast before any row work.
	// `column_index_map` is per-statement, not per-row, so these never change
	// inside the scan.
	if (study_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.experiments requires the 'study_ref' column");
	}
	if (sample_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.experiments requires the 'sample_descriptor' column");
	}
	if (strategy_idx == DConstants::INVALID_INDEX || source_idx == DConstants::INVALID_INDEX ||
	    selection_idx == DConstants::INVALID_INDEX || layout_idx == DConstants::INVALID_INDEX ||
	    platform_idx == DConstants::INVALID_INDEX || instrument_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException(
		    "INSERT INTO ena.experiments requires library_strategy, library_source, library_selection, "
		    "library_layout, platform, and instrument_model");
	}

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::ExperimentSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.experiments: 'alias' must be non-empty");
			}

			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			const auto study_val = ValueToVarchar(chunk.data[study_idx].GetValue(row));
			if (study_val.empty()) {
				throw InvalidInputException("INSERT INTO ena.experiments: 'study_ref' must be non-empty (alias '%s')",
				                            spec.alias);
			}
			spec.study_ref = ToRef(study_val, {"PRJEB", "PRJNA", "PRJDB", "ERP"});
			const auto sample_val = ValueToVarchar(chunk.data[sample_idx].GetValue(row));
			if (sample_val.empty()) {
				throw InvalidInputException(
				    "INSERT INTO ena.experiments: 'sample_descriptor' must be non-empty (alias '%s')", spec.alias);
			}
			spec.sample_ref = ToRef(sample_val, {"ERS", "SAMEA", "SAMN", "SAMD"});
			if (design_idx != DConstants::INVALID_INDEX) {
				spec.design_description = ValueToVarchar(chunk.data[design_idx].GetValue(row));
			}
			if (libname_idx != DConstants::INVALID_INDEX) {
				spec.library_name = ValueToVarchar(chunk.data[libname_idx].GetValue(row));
			}
			spec.library_strategy = ValueToVarchar(chunk.data[strategy_idx].GetValue(row));
			spec.library_source = ValueToVarchar(chunk.data[source_idx].GetValue(row));
			spec.library_selection = ValueToVarchar(chunk.data[selection_idx].GetValue(row));
			spec.library_layout = ParseLayout(spec.alias, ValueToVarchar(chunk.data[layout_idx].GetValue(row)));
			spec.platform = ValueToVarchar(chunk.data[platform_idx].GetValue(row));
			spec.instrument_model = ValueToVarchar(chunk.data[instrument_idx].GetValue(row));

			out.specs.push_back(std::move(spec));
		}
	}
	return out;
}

void ENAExperimentsInsert::AppendReturningRows(ColumnDataCollection &return_collection,
                                               const vector<LogicalType> &return_types,
                                               const std::vector<miint::ExperimentSpec> &specs,
                                               const miint::ENAExperimentSubmissionOutcome &outcome) {
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
		chunk.data[COL_STUDY_REF].SetValue(
		    idx, Value(spec.study_ref.accession.empty() ? spec.study_ref.refname : spec.study_ref.accession));
		chunk.data[COL_SAMPLE_DESCRIPTOR].SetValue(
		    idx, Value(spec.sample_ref.accession.empty() ? spec.sample_ref.refname : spec.sample_ref.accession));
		chunk.data[COL_DESIGN_DESCRIPTION].SetValue(idx, Value(spec.design_description));
		chunk.data[COL_LIBRARY_NAME].SetValue(idx, Value(spec.library_name));
		chunk.data[COL_LIBRARY_STRATEGY].SetValue(idx, Value(spec.library_strategy));
		chunk.data[COL_LIBRARY_SOURCE].SetValue(idx, Value(spec.library_source));
		chunk.data[COL_LIBRARY_SELECTION].SetValue(idx, Value(spec.library_selection));
		chunk.data[COL_LIBRARY_LAYOUT].SetValue(
		    idx, Value(spec.library_layout == miint::ENALibraryLayout::PAIRED ? "PAIRED" : "SINGLE"));
		chunk.data[COL_PLATFORM].SetValue(idx, Value(spec.platform));
		chunk.data[COL_INSTRUMENT_MODEL].SetValue(idx, Value(spec.instrument_model));
		chunk.data[COL_ERX_ACCESSION].SetValue(idx, row.erx_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                      : Value(row.erx_accession));
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
