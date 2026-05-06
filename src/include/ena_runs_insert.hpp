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

using ENARunInsertOptions = ENAInsertOptions;

struct ENARunInsertResult {
	std::string alias;
	std::string err_accession;
	std::string status;
	std::string hold_until_date;
};

using ENARunSubmissionOutcome = ENABaseSubmissionOutcome<ENARunInsertResult>;

ENARunSubmissionOutcome SubmitRunInsertOutcome(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                               const ENAPostFn &post_fn);

// Convenience wrapper used by the Catch2 unit tests — discards the outcome
// envelope and returns just the row vector. Throws on logical failure.
std::vector<ENARunInsertResult> SubmitRunInsert(const std::vector<RunSpec> &runs, const ENARunInsertOptions &opts,
                                                const ENAPostFn &post_fn);

} // namespace miint
