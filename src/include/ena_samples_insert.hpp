// SPDX-License-Identifier: MIT
//
// Pure-data helper for INSERT INTO ena.samples.
//
// Mirrors the ena_projects_insert.hpp shape so the PhysicalOperator side
// can stay thin and tests can drive the round-trip with an injected post
// functor (no DuckDB linkage required by the unit tests).

#pragma once

#include "ena_envelope_builder.hpp"
#include "ena_post_fn.hpp"

#include <string>
#include <vector>

namespace miint {

using ENASampleInsertOptions = ENAInsertOptions;

struct ENASampleInsertResult {
	std::string alias;
	std::string ers_accession;   // ERS… (primary)
	std::string samea_accession; // SAMEA… (BioSample EXT_ID)
	std::string status;
	std::string hold_until_date;
};

// Outcome-rich result mirroring ENASubmissionOutcome but for the samples
// receipt path. Reuses the same envelope/raw-receipt/error_messages fields
// since submission_log doesn't distinguish object types beyond the row's
// `object_type` column.
using ENASamplesSubmissionOutcome = ENABaseSubmissionOutcome<ENASampleInsertResult>;

// Build envelope, POST, parse receipt, return per-row results paired with
// the originating SampleSpecs by alias. Empty `samples` is a no-op (no POST,
// empty result). Logical failures (server success=false, parse errors,
// missing alias) are surfaced as `success=false` with `error_messages`
// populated and `raw_receipt` retained — matching the project insert
// contract so submission_log stays useful on parse failure.
ENASamplesSubmissionOutcome SubmitSampleInsertOutcome(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn);

// Convenience wrapper used by the Catch2 unit tests — discards the outcome
// envelope and returns just the row vector. Throws on logical failure.
std::vector<ENASampleInsertResult> SubmitSampleInsert(const std::vector<SampleSpec> &samples,
                                                      const ENASampleInsertOptions &opts, const ENAPostFn &post_fn);

} // namespace miint
