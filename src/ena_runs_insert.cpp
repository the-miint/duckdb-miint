// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.runs. See ena_runs_insert.hpp for the
// contract.

#include "ena_runs_insert.hpp"

#include "ena_receipt_parser.hpp"
#include "ena_submit_outcome.hpp"

#include <stdexcept>

namespace miint {

namespace {

struct RunSubmitTraits {
	static void SetEnvelopeArray(SubmissionSpec &env, const std::vector<RunSpec> &specs) {
		env.runs = specs;
	}
	static const char *ReceiptObjectType() {
		return "RUN";
	}
	// V2 server's JSON dispatcher NPEs for SRA-side objects — XML required.
	static std::string BuildEnvelope(const SubmissionSpec &env) {
		return BuildEnvelopeXML(env);
	}
	static const char *ContentType() {
		return "application/xml";
	}
	static ENARunInsertResult BuildRow(const RunSpec &spec, const ENAObjectReceipt &obj) {
		ENARunInsertResult row;
		row.alias = spec.alias;
		row.err_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		return row;
	}
};

} // namespace

ENARunSubmissionOutcome SubmitRunInsertOutcome(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                               const ENAPostFn &post_fn) {
	return SubmitENAObjectOutcome<RunSubmitTraits, RunSpec, ENARunInsertOptions, ENARunSubmissionOutcome>(runs, opts,
	                                                                                                      post_fn);
}

std::vector<ENARunInsertResult> SubmitRunInsert(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                                const ENAPostFn &post_fn) {
	auto outcome = SubmitRunInsertOutcome(runs, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenENAErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
