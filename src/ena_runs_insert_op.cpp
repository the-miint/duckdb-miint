// SPDX-License-Identifier: MIT
//
// Run-specific data mapping for ENARunsInsert. The PhysicalOperator
// scaffolding (Sink, Source, Finalize, GlobalSinkState, GlobalSourceState)
// lives in ENAObjectInsertOperator (ena_object_insert_op.hpp). Keep this file
// focused on the run-row → RunSpec parse and the RETURNING projection.

#include "ena_runs_insert_op.hpp"

#include "ena_insert_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"

#include <cctype>

namespace duckdb {

namespace {

// Column indices match `BuildENATableInfo(RUNS)` in src/ena_storage.cpp.
// Keep both in sync — the order of columns is the contract.
constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_EXPERIMENT_REF = 1;
constexpr idx_t COL_TITLE = 2;
constexpr idx_t COL_FILES = 3;
constexpr idx_t COL_ERR_ACCESSION = 4;
constexpr idx_t TABLE_COLUMN_COUNT = 5;

// experiment_ref strings shaped like an ENA accession (`ERX<digits>`) route to
// RefDescriptor::accession; everything else is treated as refname (the
// experiment's alias). The "trailing digits only" check matters: a user alias
// like `ERXperiment-1` would otherwise be silently misclassified as an
// accession the server can't find.
miint::RefDescriptor ToExperimentRef(const std::string &value) {
	miint::RefDescriptor ref;
	if (value.size() > 3 && value.compare(0, 3, "ERX") == 0) {
		bool all_digits = true;
		for (size_t i = 3; i < value.size(); ++i) {
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

// Pull a LIST(STRUCT(filename,filetype,md5)) value into the RunFile vector.
// NULL or empty list → empty vector (the envelope builder rejects empty
// `files` with a clear "must have at least one file" error).
std::vector<miint::RunFile> ExtractFilesList(const Value &v, const std::string &alias) {
	std::vector<miint::RunFile> out;
	if (v.IsNull()) {
		return out;
	}
	const auto &entries = ListValue::GetChildren(v);
	out.reserve(entries.size());
	for (idx_t i = 0; i < entries.size(); i++) {
		const auto &entry = entries[i];
		// A NULL entry (e.g. `[NULL]` or `[{...}, NULL]`) would otherwise
		// reach StructValue::GetChildren on a null Value, which is undefined
		// behaviour. Surface a clear error instead.
		if (entry.IsNull()) {
			throw InvalidInputException(
			    "INSERT INTO ena.runs: 'files' list contains NULL at position %llu (alias '%s')",
			    static_cast<unsigned long long>(i), alias);
		}
		const auto &fields = StructValue::GetChildren(entry);
		if (fields.size() != 3) {
			throw InvalidInputException("INSERT INTO ena.runs: 'files' STRUCT must have (filename,filetype,md5) for "
			                            "alias '%s'",
			                            alias);
		}
		miint::RunFile f;
		f.filename = ValueToVarchar(fields[0]);
		f.filetype = ValueToVarchar(fields[1]);
		f.checksum = ValueToVarchar(fields[2]);
		if (f.filename.empty()) {
			throw InvalidInputException("INSERT INTO ena.runs: 'files.filename' must be non-empty (alias '%s')", alias);
		}
		if (f.checksum.empty()) {
			throw InvalidInputException(
			    "INSERT INTO ena.runs: 'files.md5' must be non-empty (alias '%s', filename '%s')", alias, f.filename);
		}
		out.push_back(std::move(f));
	}
	return out;
}

} // namespace

ENABuiltSpecs<miint::RunSpec> ENARunsInsert::BuildFromBuffer(ColumnDataCollection &buffer,
                                                             const physical_index_vector_t<idx_t> &column_index_map) {
	ENABuiltSpecs<miint::RunSpec> out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.runs requires the 'alias' column");
	}
	const idx_t exp_ref_idx = ResolveInputColumn(column_index_map, input_columns, COL_EXPERIMENT_REF);
	if (exp_ref_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.runs requires the 'experiment_ref' column");
	}
	const idx_t files_idx = ResolveInputColumn(column_index_map, input_columns, COL_FILES);
	if (files_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.runs requires the 'files' column");
	}
	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::RunSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.runs: 'alias' must be non-empty");
			}
			const auto exp_val = ValueToVarchar(chunk.data[exp_ref_idx].GetValue(row));
			if (exp_val.empty()) {
				throw InvalidInputException("INSERT INTO ena.runs: 'experiment_ref' must be non-empty (alias '%s')",
				                            spec.alias);
			}
			spec.experiment_ref = ToExperimentRef(exp_val);
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			spec.files = ExtractFilesList(chunk.data[files_idx].GetValue(row), spec.alias);
			if (spec.files.empty()) {
				throw InvalidInputException("INSERT INTO ena.runs: 'files' must contain at least one file (alias '%s')",
				                            spec.alias);
			}
			out.specs.push_back(std::move(spec));
		}
	}
	return out;
}

void ENARunsInsert::AppendReturningRows(ColumnDataCollection &return_collection,
                                        const vector<LogicalType> &return_types,
                                        const std::vector<miint::RunSpec> &specs,
                                        const miint::ENARunSubmissionOutcome &outcome) {
	D_ASSERT(return_types.size() == TABLE_COLUMN_COUNT);
	D_ASSERT(outcome.rows.size() == specs.size());

	DataChunk chunk;
	chunk.Initialize(Allocator::DefaultAllocator(), return_types);
	for (idx_t i = 0; i < outcome.rows.size(); i++) {
		const auto &row = outcome.rows[i];
		const auto &spec = specs[i];
		const auto idx = chunk.size();
		chunk.data[COL_ALIAS].SetValue(idx, Value(row.alias));
		chunk.data[COL_EXPERIMENT_REF].SetValue(
		    idx,
		    Value(spec.experiment_ref.accession.empty() ? spec.experiment_ref.refname : spec.experiment_ref.accession));
		chunk.data[COL_TITLE].SetValue(idx, Value(spec.title));
		// `files` echoed as NULL on RETURNING — the user-supplied LIST is
		// preserved verbatim in `ena.submission_log.request_payload`. Real
		// LIST(STRUCT) emission is a future task (mirrors the projects /
		// samples attribute MAP trade-off).
		chunk.data[COL_FILES].SetValue(idx, Value(return_types[COL_FILES]));
		chunk.data[COL_ERR_ACCESSION].SetValue(idx, row.err_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                      : Value(row.err_accession));
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
