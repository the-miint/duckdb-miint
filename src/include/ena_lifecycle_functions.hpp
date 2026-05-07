// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 targeted lifecycle actions:
//   - ena_cancel(secret, accession, catalog?)
//   - ena_release(secret, accession, catalog?)
//   - ena_hold(secret, accession, until, catalog?)
//
// All three return a single row with the same shape:
//   action, target, success, era_accession, hold_until_date,
//   error_messages, duration_ms
//
// `hold_until_date` is empty for CANCEL/RELEASE.
//
// `accession` must be the server-assigned ID (PRJEB / ERS / ERX / ERR / ERZ)
// from the original registration. Refname/alias targeting is not supported
// by Webin V2 for cross-submission lifecycle ops; verified live 2026-05-07.
// Within a session, an accession can be recovered from an alias via
// `ena.submission_log.object_aliases` / `object_accessions` (see
// docs/ena.md).
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
