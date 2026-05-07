// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 targeted lifecycle actions:
//   - ena_cancel(secret, accession?, refname?, catalog?)
//   - ena_release(secret, accession?, refname?, catalog?)
//   - ena_hold(secret, accession?, refname?, until, catalog?)
//
// All three return a single row with the same shape:
//   action, target, success, era_accession, hold_until_date,
//   error_messages, duration_ms
//
// `hold_until_date` is empty for CANCEL/RELEASE.
//
// Submission_log integration: if a catalog of type ena is attached (default
// name 'ena', overridable via the `catalog` named param), each invocation
// appends a submission_log row to that catalog. If no matching catalog is
// attached, the function still works but does not log — useful for one-shot
// CANCEL operations from a session without a full ena catalog setup.

#pragma once

namespace duckdb {
class ExtensionLoader;
}

namespace miint {

void RegisterENALifecycleTableFunctions(duckdb::ExtensionLoader &loader);

} // namespace miint
