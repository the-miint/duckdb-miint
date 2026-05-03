// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.projects. See ena_projects_insert.hpp.

#include "ena_projects_insert.hpp"

#include "ena_receipt_parser.hpp"

#include <chrono>
#include <stdexcept>
#include <unordered_map>

namespace miint {

namespace {

std::string FlattenErrors(const std::vector<std::string> &errs) {
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

ENASubmissionOutcome SubmitProjectInsertOutcome(const std::vector<ProjectSpec> &projects,
                                                const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn) {
	ENASubmissionOutcome outcome;
	if (projects.empty()) {
		outcome.success = true;
		return outcome;
	}

	SubmissionSpec env;
	env.action = ENAAction::ADD;
	env.hold_until_date = opts.hold_until_date;
	env.projects = projects;

	outcome.envelope_payload = BuildEnvelopeJSON(env);

	const auto t_start = std::chrono::steady_clock::now();
	outcome.raw_receipt =
	    post_fn(opts.endpoint_url, outcome.envelope_payload, opts.user, opts.password, "application/json");
	const auto t_end = std::chrono::steady_clock::now();
	outcome.duration_ms = std::chrono::duration_cast<std::chrono::milliseconds>(t_end - t_start).count();

	// Parse the receipt. Surface a parse failure as success=false with the raw
	// body retained so the submission_log row preserves the server response.
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
		// Server reported failure; raw_receipt is already on the outcome.
		return outcome;
	}

	// Index receipt PROJECT entries by alias for joining back to the request.
	std::unordered_map<std::string, const ENAObjectReceipt *> by_alias;
	by_alias.reserve(receipt.objects.size());
	for (const auto &obj : receipt.objects) {
		if (obj.object_type == "PROJECT") {
			by_alias[obj.alias] = &obj;
		}
	}

	outcome.rows.reserve(projects.size());
	for (const auto &spec : projects) {
		auto it = by_alias.find(spec.alias);
		if (it == by_alias.end()) {
			outcome.success = false;
			outcome.error_messages.push_back("ENA submission: receipt missing PROJECT entry for alias '" + spec.alias +
			                                 "'");
			return outcome;
		}
		const auto &obj = *it->second;
		ENAProjectInsertResult row;
		row.alias = spec.alias;
		row.prjeb_accession = obj.accession;
		row.status = obj.status;
		row.hold_until_date = obj.hold_until_date;
		for (const auto &ext : obj.ext_ids) {
			if (ext.type == "study") {
				row.erp_accession = ext.accession;
				break;
			}
		}
		outcome.rows.push_back(std::move(row));
	}
	return outcome;
}

std::vector<ENAProjectInsertResult> SubmitProjectInsert(const std::vector<ProjectSpec> &projects,
                                                        const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn) {
	auto outcome = SubmitProjectInsertOutcome(projects, opts, post_fn);
	if (!outcome.success) {
		const std::string detail =
		    outcome.error_messages.empty() ? "no error detail" : FlattenErrors(outcome.error_messages);
		throw std::runtime_error("ENA submission failed: " + detail);
	}
	return outcome.rows;
}

} // namespace miint
