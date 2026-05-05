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
#include "ena_post_fn.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

using ENAProjectInsertOptions = ENAInsertOptions;

struct ENAProjectInsertResult {
	std::string alias;
	std::string prjeb_accession;
	std::string erp_accession;
	std::string status;
	std::string hold_until_date;
};

// Outcome-rich result: contains the per-row insert results plus the receipt
// metadata needed to populate ena.submission_log.
using ENASubmissionOutcome = ENABaseSubmissionOutcome<ENAProjectInsertResult>;

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
