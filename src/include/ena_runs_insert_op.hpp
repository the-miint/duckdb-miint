// SPDX-License-Identifier: MIT
//
// PhysicalOperator for `INSERT INTO ena.runs ... [RETURNING ...]`.
//
// Sink/Source/Finalize machinery lives in ENAObjectInsertOperator (the CRTP
// base in ena_object_insert_op.hpp). This subclass only owns the per-table
// data mapping: chunk row → RunSpec, RETURNING column population, and the
// per-table identifier strings.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_object_insert_op.hpp"
#include "ena_runs_insert.hpp"

namespace duckdb {

class ENARunsInsert : public ENAObjectInsertOperator<ENARunsInsert, miint::RunSpec, miint::ENARunInsertOptions,
                                                     miint::ENARunSubmissionOutcome> {
public:
	using ENAObjectInsertOperator::ENAObjectInsertOperator;

	static const char *ObjectName() {
		return "runs";
	}
	static const char *ThrowPrefix() {
		return "INSERT INTO ena.runs";
	}
	static string OperatorName() {
		return "ENA_RUNS_INSERT";
	}
	static const std::string &PrimaryAccession(const miint::ENARunInsertResult &row) {
		return row.err_accession;
	}

	static ENABuiltSpecs<miint::RunSpec> BuildFromBuffer(ColumnDataCollection &buffer,
	                                                     const physical_index_vector_t<idx_t> &column_index_map);

	static miint::ENARunSubmissionOutcome Submit(const std::vector<miint::RunSpec> &specs,
	                                             const miint::ENARunInsertOptions &opts,
	                                             const miint::ENAPostFn &post_fn) {
		return miint::SubmitRunInsertOutcome(specs, opts, post_fn);
	}

	static void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
	                                const std::vector<miint::RunSpec> &specs,
	                                const miint::ENARunSubmissionOutcome &outcome);
};

} // namespace duckdb
