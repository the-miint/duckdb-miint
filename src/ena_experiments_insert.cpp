// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.experiments. See
// ena_experiments_insert.hpp for the contract.

#include "ena_experiments_insert.hpp"

#include "ena_receipt_parser.hpp"
#include "ena_submit_outcome.hpp"

#include <stdexcept>

namespace miint {

namespace {

struct ExperimentSubmitTraits {
	static void SetEnvelopeArray(SubmissionSpec &env, const std::vector<ExperimentSpec> &specs) {
		env.experiments = specs;
	}
	static const char *ReceiptObjectType() {
		return "EXPERIMENT";
	}
	static ENAExperimentInsertResult BuildRow(const ExperimentSpec &spec, const ENAObjectReceipt &obj) {
		ENAExperimentInsertResult row;
		row.alias = spec.alias;
		row.erx_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		return row;
	}
};

} // namespace

ENAExperimentSubmissionOutcome SubmitExperimentInsertOutcome(const std::vector<ExperimentSpec> &experiments,
                                                             const ENAExperimentInsertOptions &opts,
                                                             const ENAPostFn &post_fn) {
	return SubmitENAObjectOutcome<ExperimentSubmitTraits, ExperimentSpec, ENAExperimentInsertOptions,
	                              ENAExperimentSubmissionOutcome>(experiments, opts, post_fn);
}

std::vector<ENAExperimentInsertResult> SubmitExperimentInsert(const std::vector<ExperimentSpec> &experiments,
                                                              const ENAExperimentInsertOptions &opts,
                                                              const ENAPostFn &post_fn) {
	auto outcome = SubmitExperimentInsertOutcome(experiments, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenENAErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
