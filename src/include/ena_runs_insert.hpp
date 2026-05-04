// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.runs. Mirrors the projects / samples /
// experiments shape: assemble the envelope, POST via injected functor, parse
// the receipt, return per-row results.

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_post_fn.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct ENARunInsertOptions {
	std::string endpoint_url; // resolved base URL incl. /submit suffix
	std::string user;         // Webin-XXXXX
	std::string password;
	std::string hold_until_date; // optional, "YYYY-MM-DD"
};

struct ENARunInsertResult {
	std::string alias;
	std::string err_accession;
	std::string status;
	std::string hold_until_date;
};

struct ENARunSubmissionOutcome {
	std::vector<ENARunInsertResult> rows;
	std::string envelope_payload;
	std::string raw_receipt;
	std::string era_accession;
	bool success = false;
	std::vector<std::string> error_messages;
	int64_t duration_ms = 0;
};

ENARunSubmissionOutcome SubmitRunInsertOutcome(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                               const ENAPostFn &post_fn);

// Convenience wrapper used by the Catch2 unit tests — discards the outcome
// envelope and returns just the row vector. Throws on logical failure.
std::vector<ENARunInsertResult> SubmitRunInsert(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                                const ENAPostFn &post_fn);

} // namespace miint
