// SPDX-License-Identifier: MIT
//
// PhysicalOperator for `INSERT INTO ena.samples ... [RETURNING ...]`.
//
// Sink/Source/Finalize machinery lives in ENAObjectInsertOperator (the CRTP
// base in ena_object_insert_op.hpp). This subclass only owns the per-table
// data mapping: chunk row → SampleSpec, RETURNING column population, and the
// per-table identifier strings.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_object_insert_op.hpp"
#include "ena_samples_insert.hpp"

namespace duckdb {

class ENASamplesInsert
    : public ENAObjectInsertOperator<ENASamplesInsert, miint::SampleSpec, miint::ENASampleInsertOptions,
                                     miint::ENASamplesSubmissionOutcome> {
public:
	using ENAObjectInsertOperator::ENAObjectInsertOperator;

	static const char *ObjectName() {
		return "samples";
	}
	static const char *ThrowPrefix() {
		return "INSERT INTO ena.samples";
	}
	static string OperatorName() {
		return "ENA_SAMPLES_INSERT";
	}

	static ENABuiltSpecs<miint::SampleSpec> BuildFromBuffer(ColumnDataCollection &buffer,
	                                                        const physical_index_vector_t<idx_t> &column_index_map);

	static miint::ENASamplesSubmissionOutcome Submit(const std::vector<miint::SampleSpec> &specs,
	                                                 const miint::ENASampleInsertOptions &opts,
	                                                 const miint::ENAPostFn &post_fn) {
		return miint::SubmitSampleInsertOutcome(specs, opts, post_fn);
	}

	// Phase 8 Step 8b: validate user-supplied attribute MAPs against the
	// declared ENA checklist (mandatory fields, units, controlled
	// vocabularies). Fetches the checklist XML from the ENA browser API
	// (or the MIINT_ENA_CHECKLIST_URL_BASE override) and caches it for the
	// process lifetime. Specs with empty `checklist` are skipped — the
	// user has opted out of client-side validation.
	static void ValidateBuiltSpecs(const std::vector<miint::SampleSpec> &specs, miint::ENAClient &client);

	static void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
	                                const std::vector<miint::SampleSpec> &specs,
	                                const miint::ENASamplesSubmissionOutcome &outcome);
};

} // namespace duckdb
