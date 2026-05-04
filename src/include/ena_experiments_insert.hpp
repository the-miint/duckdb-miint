// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.experiments. Mirrors
// ena_projects_insert.hpp / ena_samples_insert.hpp so the PhysicalOperator
// side can stay thin and tests can drive the round-trip with an injected
// post functor (no DuckDB linkage required by the unit tests).

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_post_fn.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct ENAExperimentInsertOptions {
	std::string endpoint_url; // resolved base URL incl. /submit suffix
	std::string user;         // Webin-XXXXX
	std::string password;
	std::string hold_until_date; // optional, "YYYY-MM-DD"
};

struct ENAExperimentInsertResult {
	std::string alias;
	std::string erx_accession;
	std::string status;
	std::string hold_until_date;
};

// Outcome-rich result mirroring ENASubmissionOutcome but for the experiments
// receipt path. Reuses the same envelope/raw-receipt/error_messages fields
// since submission_log doesn't distinguish object types beyond the row's
// `object_type` column.
struct ENAExperimentSubmissionOutcome {
	std::vector<ENAExperimentInsertResult> rows;
	std::string envelope_payload;
	std::string raw_receipt;
	std::string era_accession;
	bool success = false;
	std::vector<std::string> error_messages;
	int64_t duration_ms = 0;
};

// Build envelope, POST, parse receipt, return per-row results paired with the
// originating ExperimentSpecs by alias. Empty `experiments` is a no-op.
// Logical failures (server success=false, parse errors, missing alias) are
// surfaced as `success=false` with `error_messages` populated and
// `raw_receipt` retained — matching the project/sample insert contract so
// submission_log stays useful on parse failure.
ENAExperimentSubmissionOutcome SubmitExperimentInsertOutcome(const std::vector<ExperimentSpec> &experiments,
                                                             const ENAExperimentInsertOptions &opts,
                                                             const ENAPostFn &post_fn);

// Convenience wrapper used by the Catch2 unit tests — discards the outcome
// envelope and returns just the row vector. Throws on logical failure.
std::vector<ENAExperimentInsertResult> SubmitExperimentInsert(const std::vector<ExperimentSpec> &experiments,
                                                              const ENAExperimentInsertOptions &opts,
                                                              const ENAPostFn &post_fn);

} // namespace miint
