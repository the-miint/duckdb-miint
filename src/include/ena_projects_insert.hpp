// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.projects.
//
// Wraps envelope build + HTTP POST + receipt parse so the Sink/Finalize side
// of the PhysicalOperator can stay thin and the round-trip is unit-testable
// with an injected post functor (no DuckDB linkage required by the unit
// tests).

#pragma once

#include "ena_envelope_builder.hpp"

#include <cstdint>
#include <functional>
#include <string>
#include <vector>

namespace miint {

struct ENAProjectInsertOptions {
	std::string endpoint_url; // resolved base URL incl. /submit suffix
	std::string user;         // Webin-XXXXX
	std::string password;
	std::string hold_until_date; // optional, "YYYY-MM-DD"
};

struct ENAProjectInsertResult {
	std::string alias;
	std::string prjeb_accession;
	std::string erp_accession;
	std::string status;
	std::string hold_until_date;
};

// Functor signature: (url, body, user, password, content_type) -> response_body.
// Production wires this to ENAClient::PostJSON; tests inject a fake.
using ENAPostFn = std::function<std::string(const std::string &url, const std::string &body, const std::string &user,
                                            const std::string &password, const std::string &content_type)>;

// Outcome-rich result: contains the per-row insert results plus the receipt
// metadata needed to populate ena.submission_log.
struct ENASubmissionOutcome {
	std::vector<ENAProjectInsertResult> rows;
	std::string envelope_payload; // the JSON we POSTed (passwords not embedded)
	std::string raw_receipt;      // the raw response body
	std::string era_accession;    // server-assigned <SUBMISSION accession>
	bool success = false;
	std::vector<std::string> error_messages;
	int64_t duration_ms = 0;
};

// Build envelope, POST it, parse receipt, return per-row results paired with
// the originating ProjectSpecs by alias. Empty `projects` is a no-op (no POST,
// empty result). Throws std::runtime_error if the receipt reports
// success=false (the message text includes any server-side error messages),
// or if the receipt is missing accessions for one or more requested aliases.
ENASubmissionOutcome SubmitProjectInsertOutcome(const std::vector<ProjectSpec> &projects,
                                                const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn);

// Convenience wrapper used by the Catch2 unit tests — discards the outcome
// envelope and returns just the row vector. Throws on receipt failure.
std::vector<ENAProjectInsertResult> SubmitProjectInsert(const std::vector<ProjectSpec> &projects,
                                                        const ENAProjectInsertOptions &opts, const ENAPostFn &post_fn);

} // namespace miint
