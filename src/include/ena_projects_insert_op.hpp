// SPDX-License-Identifier: MIT
//
// PhysicalOperator for `INSERT INTO ena.projects ... [RETURNING ...]`.
//
// Sink/Source/Finalize machinery lives in ENAObjectInsertOperator (the CRTP
// base in ena_object_insert_op.hpp). This subclass only owns the per-table
// data mapping: chunk row → ProjectSpec, RETURNING column population, and the
// per-table identifier strings.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_object_insert_op.hpp"
#include "ena_projects_insert.hpp"

namespace duckdb {

class ENAProjectsInsert : public ENAObjectInsertOperator<ENAProjectsInsert, miint::ProjectSpec,
                                                         miint::ENAProjectInsertOptions, miint::ENASubmissionOutcome> {
public:
	using ENAObjectInsertOperator::ENAObjectInsertOperator;

	static const char *ObjectName() {
		return "projects";
	}
	static const char *ThrowPrefix() {
		return "INSERT INTO ena.projects";
	}
	static string OperatorName() {
		return "ENA_PROJECTS_INSERT";
	}

	static ENABuiltSpecs<miint::ProjectSpec> BuildFromBuffer(ColumnDataCollection &buffer,
	                                                         const physical_index_vector_t<idx_t> &column_index_map);

	static miint::ENASubmissionOutcome Submit(const std::vector<miint::ProjectSpec> &specs,
	                                          const miint::ENAProjectInsertOptions &opts,
	                                          const miint::ENAPostFn &post_fn) {
		return miint::SubmitProjectInsertOutcome(specs, opts, post_fn);
	}

	static void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
	                                const std::vector<miint::ProjectSpec> &specs,
	                                const miint::ENASubmissionOutcome &outcome);
};

} // namespace duckdb
