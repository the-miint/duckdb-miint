// SPDX-License-Identifier: MIT
// Pure-data helper for INSERT INTO ena.samples. See ena_samples_insert.hpp.
//
// Phase 5 stub — implementation lands in the Green step.

#include "ena_samples_insert.hpp"

#include "ena_receipt_parser.hpp"

#include <chrono>
#include <stdexcept>
#include <unordered_map>

namespace miint {

namespace {

std::string FlattenSampleErrors(const std::vector<std::string> &errs) {
	std::string out;
	for (size_t i = 0; i < errs.size(); ++i) {
		if (i > 0) {
			out += "; ";
		}
		out += errs[i];
	}
	return out;
}

} // namespace

ENASamplesSubmissionOutcome SubmitSampleInsertOutcome(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn) {
	ENASamplesSubmissionOutcome outcome;
	if (samples.empty()) {
		outcome.success = true;
		return outcome;
	}

	SubmissionSpec env;
	env.action = ENAAction::ADD;
	env.hold_until_date = opts.hold_until_date;
	env.samples = samples;

	outcome.envelope_payload = BuildEnvelopeJSON(env);

	const auto t_start = std::chrono::steady_clock::now();
	outcome.raw_receipt =
	    post_fn(opts.endpoint_url, outcome.envelope_payload, opts.user, opts.password, "application/json");
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

	if (!receipt.success) {
		return outcome;
	}

	std::unordered_map<std::string, const ENAObjectReceipt *> by_alias;
	by_alias.reserve(receipt.objects.size());
	for (const auto &obj : receipt.objects) {
		if (obj.object_type == "SAMPLE") {
			by_alias[obj.alias] = &obj;
		}
	}

	outcome.rows.reserve(samples.size());
	for (const auto &spec : samples) {
		auto it = by_alias.find(spec.alias);
		if (it == by_alias.end()) {
			outcome.success = false;
			outcome.error_messages.push_back("ENA submission: receipt missing SAMPLE entry for alias '" + spec.alias +
			                                 "'");
			return outcome;
		}
		const auto &obj = *it->second;
		ENASampleInsertResult row;
		row.alias = spec.alias;
		row.ers_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		for (const auto &ext : obj.ext_ids) {
			// Case-insensitive: the ENA receipt XSD enumerates this attribute
			// as a free-form string, and seen casings include "biosample" (V2),
			// "BioSample" (some legacy paths), and "BIOSAMPLE" (rare).
			if (EqualsIgnoreCase(ext.type, "biosample")) {
				row.samea_accession = ext.accession;
				break;
			}
		}
		outcome.rows.push_back(std::move(row));
	}
	return outcome;
}

std::vector<ENASampleInsertResult> SubmitSampleInsert(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn) {
	auto outcome = SubmitSampleInsertOutcome(samples, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenSampleErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
