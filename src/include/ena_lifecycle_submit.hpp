// SPDX-License-Identifier: MIT
//
// Pure-data submit core for ENA Webin V2 *targeted* lifecycle actions:
// CANCEL, RELEASE, and post-hoc HOLD. ADD/MODIFY/VALIDATE go through the
// per-table SubmitXInsertOutcome path (they carry a body); this file is
// only for the actions that operate on an existing accession via
// `target=` on the action element.
//
// Usage: callers (the ena_cancel/ena_hold/ena_release table functions and
// ENACatalog::PlanDelete dispatch) build a LifecycleSubmitOptions, hand it
// to SubmitLifecycle, and consume LifecycleOutcome. Receipt parsing
// failures and server `success=false` responses populate
// `outcome.error_messages` and leave `outcome.success=false`; transport
// failures propagate as exceptions from the post functor.

#pragma once

#include "ena_envelope_builder.hpp" // ENAAction
#include "ena_post_fn.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct LifecycleSubmitOptions {
	std::string endpoint_url; // resolved base URL incl. /submit suffix
	std::string user;         // Webin-XXXXX
	std::string password;
	// Target the existing object via accession (preferred) or refname/alias.
	// Reusing RefDescriptor keeps the precedence rule ("accession wins when
	// both are set") in one place across the codebase.
	RefDescriptor target;
	// Required for action=HOLD. Rejected for CANCEL/RELEASE — those actions
	// don't take a date and silently echoing one through into submission_log
	// would mislead audit consumers.
	std::string hold_until_date;
};

struct LifecycleOutcome {
	ENAAction action;                        // echoed from input for downstream logging
	std::string target;                      // resolved (accession ?? refname) — what we sent
	std::string hold_until_date;             // echoed from input on HOLD; empty otherwise
	std::string envelope_payload;            // request body — useful for submission_log audit
	std::string raw_receipt;                 // server response — preserved even if parse fails
	std::string era_accession;               // server-assigned ERA submission id
	bool success = false;                    // true iff receipt explicitly reports success
	std::vector<std::string> error_messages; // server errors + parse-error wrapper
	int64_t duration_ms = 0;                 // wall-clock POST round-trip
};

// Build envelope, POST via `post_fn`, parse receipt, return outcome.
//
// Throws std::runtime_error before any network activity if:
//   - `action` is not one of CANCEL / RELEASE / HOLD (others have body shapes
//     that go through the per-table insert path; this is a programming error
//     surfaced loudly so a misrouted call fails fast).
//   - `opts` violates the envelope-builder invariants (missing target, missing
//     date for HOLD, whitespace-only target, etc.). Validation is delegated
//     to BuildEnvelopeXML — see ena_envelope_builder.cpp::ValidateActions.
//
// Receipt parse errors and server `success=false` responses are non-throwing:
// they populate `outcome.error_messages` and leave `outcome.success=false`.
// This lets callers (notably PlanDelete) preserve the receipt in
// submission_log even when the action failed.
//
// Transport failures (the post functor itself raising) propagate as
// `std::exception` (or any subclass thereof — the ENAPostFn signature does
// not constrain the exception type). Callers that need a uniform place to
// read the result must catch `std::exception` at their layer.
LifecycleOutcome SubmitLifecycle(ENAAction action, const LifecycleSubmitOptions &opts, const ENAPostFn &post_fn);

} // namespace miint
