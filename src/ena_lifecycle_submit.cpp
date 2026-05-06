// SPDX-License-Identifier: MIT

#include "ena_lifecycle_submit.hpp"

#include "ena_receipt_parser.hpp"

#include <chrono>
#include <stdexcept>

namespace miint {

LifecycleOutcome SubmitLifecycle(ENAAction action, const LifecycleSubmitOptions &opts, const ENAPostFn &post_fn) {
	if (action != ENAAction::CANCEL && action != ENAAction::RELEASE && action != ENAAction::HOLD) {
		throw std::runtime_error("SubmitLifecycle: only CANCEL, RELEASE, and HOLD are supported here; "
		                         "ADD/MODIFY/VALIDATE carry a body and go through the per-table insert path");
	}

	LifecycleOutcome outcome;
	outcome.action = action;
	outcome.hold_until_date = opts.hold_until_date;
	// Echo what we actually sent (accession-or-refname after RefDescriptor
	// precedence) so submission_log captures the resolved target even when
	// the caller supplied only a refname.
	outcome.target = !opts.target.accession.empty() ? opts.target.accession : opts.target.refname;

	// BuildEnvelopeXML enforces the full validation matrix (target presence,
	// hold-date requirement / rejection, whitespace-only rejection,
	// body-content rejection). Throws std::runtime_error on violation —
	// propagates to the caller before any network activity.
	SubmissionSpec env;
	env.action = action;
	env.target_accession = opts.target.accession;
	env.target_refname = opts.target.refname;
	env.hold_until_date = opts.hold_until_date;
	outcome.envelope_payload = BuildEnvelopeXML(env);

	const auto t_start = std::chrono::steady_clock::now();
	outcome.raw_receipt =
	    post_fn(opts.endpoint_url, outcome.envelope_payload, opts.user, opts.password, "application/xml");
	const auto t_end = std::chrono::steady_clock::now();
	outcome.duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	ENAReceipt receipt;
	try {
		receipt = ParseReceiptXML(outcome.raw_receipt);
	} catch (const std::exception &e) {
		outcome.success = false;
		outcome.error_messages.push_back(std::string("ENA receipt parse failure: ") + e.what());
		return outcome;
	}
	outcome.success = receipt.success;
	outcome.era_accession = receipt.submission_accession;
	outcome.error_messages = receipt.errors;
	return outcome;
}

} // namespace miint
