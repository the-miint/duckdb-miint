// SPDX-License-Identifier: MIT
//
// Shared template for the pure-data submit-side of every ENA virtual-table
// INSERT (projects, samples, experiments, runs, analyses). Each per-table
// SubmitXInsertOutcome had near-identical structure: build envelope, POST,
// parse receipt, pair receipt rows with input by alias. The only per-table
// differences are:
//   1. Which `SubmissionSpec` array gets populated.
//   2. Which receipt object_type to filter on ("PROJECT" / "SAMPLE" / ...).
//   3. How to translate one ENAObjectReceipt + its source spec into the
//      per-table result row (which accession field gets which value).
//
// All three are captured by a `Traits` struct with three static methods.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_post_fn.hpp"
#include "ena_receipt_parser.hpp"

#include <chrono>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace miint {

// Joins a list of error messages into a single "; "-separated string. Used by
// the convenience throw-on-error wrappers each per-table file exposes.
inline std::string FlattenENAErrors(const std::vector<std::string> &errs) {
	std::string out;
	for (size_t i = 0; i < errs.size(); ++i) {
		if (i > 0) {
			out += "; ";
		}
		out += errs[i];
	}
	return out;
}

// `Traits` must provide:
//   static void SetEnvelopeArray(SubmissionSpec &env, const std::vector<SpecT> &specs);
//   static const char *ReceiptObjectType();   // "PROJECT" / "SAMPLE" / "EXPERIMENT" / "RUN"
//   static RowT BuildRow(const SpecT &spec, const ENAObjectReceipt &obj);
//   static std::string BuildEnvelope(const SubmissionSpec &env);
//   static const char *ContentType();         // "application/json" or "application/xml"
//
// V2 server caveat: project + sample go via JSON; experiment + run + analysis
// must be XML (V2's JSON dispatcher returns NPE for SRA-side objects). The
// per-table Traits chooses; the post functor dispatches by content type.
//
// `OutcomeT` must expose: success (bool), envelope_payload, raw_receipt,
//   era_accession, error_messages (vector<string>), duration_ms (int64),
//   rows (vector<RowT>).
template <class Traits, class SpecT, class OptsT, class OutcomeT>
OutcomeT SubmitENAObjectOutcome(const std::vector<SpecT> &specs, const OptsT &opts, const ENAPostFn &post_fn) {
	OutcomeT outcome;
	if (specs.empty()) {
		outcome.success = true;
		return outcome;
	}

	SubmissionSpec env;
	env.action = ENAAction::ADD;
	env.hold_until_date = opts.hold_until_date;
	Traits::SetEnvelopeArray(env, specs);
	outcome.envelope_payload = Traits::BuildEnvelope(env);

	const auto t_start = std::chrono::steady_clock::now();
	outcome.raw_receipt =
	    post_fn(opts.endpoint_url, outcome.envelope_payload, opts.user, opts.password, Traits::ContentType());
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
	const std::string target_type = Traits::ReceiptObjectType();
	for (const auto &obj : receipt.objects) {
		if (obj.object_type == target_type) {
			by_alias[obj.alias] = &obj;
		}
	}

	outcome.rows.reserve(specs.size());
	for (const auto &spec : specs) {
		auto it = by_alias.find(spec.alias);
		if (it == by_alias.end()) {
			outcome.success = false;
			outcome.error_messages.push_back("ENA submission: receipt missing " + std::string(target_type) +
			                                 " entry for alias '" + spec.alias + "'");
			return outcome;
		}
		outcome.rows.push_back(Traits::BuildRow(spec, *it->second));
	}
	return outcome;
}

} // namespace miint
