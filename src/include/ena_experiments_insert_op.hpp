// SPDX-License-Identifier: MIT
//
// PhysicalOperator for `INSERT INTO ena.experiments ... [RETURNING ...]`.
//
// Sink/Source/Finalize machinery lives in ENAObjectInsertOperator (the CRTP
// base in ena_object_insert_op.hpp). This subclass only owns the per-table
// data mapping: chunk row → ExperimentSpec, RETURNING column population, and
// the per-table identifier strings.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_experiments_insert.hpp"
#include "ena_object_insert_op.hpp"

namespace duckdb {

class ENAExperimentsInsert
    : public ENAObjectInsertOperator<ENAExperimentsInsert, miint::ExperimentSpec, miint::ENAExperimentInsertOptions,
                                     miint::ENAExperimentSubmissionOutcome> {
public:
	using ENAObjectInsertOperator::ENAObjectInsertOperator;

	static const char *ObjectName() {
		return "experiments";
	}
	static const char *ThrowPrefix() {
		return "INSERT INTO ena.experiments";
	}
	static string OperatorName() {
		return "ENA_EXPERIMENTS_INSERT";
	}
	static const std::string &PrimaryAccession(const miint::ENAExperimentInsertResult &row) {
		return row.erx_accession;
	}

	static ENABuiltSpecs<miint::ExperimentSpec> BuildFromBuffer(ColumnDataCollection &buffer,
	                                                            const physical_index_vector_t<idx_t> &column_index_map);

	static miint::ENAExperimentSubmissionOutcome Submit(const std::vector<miint::ExperimentSpec> &specs,
	                                                    const miint::ENAExperimentInsertOptions &opts,
	                                                    const miint::ENAPostFn &post_fn) {
		return miint::SubmitExperimentInsertOutcome(specs, opts, post_fn);
	}

	static void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
	                                const std::vector<miint::ExperimentSpec> &specs,
	                                const miint::ENAExperimentSubmissionOutcome &outcome);
};

} // namespace duckdb
