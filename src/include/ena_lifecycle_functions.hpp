// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 targeted lifecycle actions:
//   - ena_cancel(secret, accession | (refname, kind), catalog?)
//   - ena_release(secret, accession | (refname, kind), catalog?)
//   - ena_hold(secret, accession | (refname, kind), until, catalog?)
//
// All three return a single row with the same shape:
//   action, target, success, era_accession, hold_until_date,
//   error_messages, duration_ms
//
// `hold_until_date` is empty for CANCEL/RELEASE.
//
// Target identification — exactly one of:
//   - `accession` — server-assigned ID (PRJEB / ERS / ERX / ERR / ERZ).
//     Webin V2 lifecycle uses this directly. Kind-agnostic — the prefix
//     encodes the object kind.
//   - `refname` — user-supplied alias from the original registration.
//     Translated to the server-assigned accession at execute time via the
//     Webin Reports API (L5). Requires the sibling `kind` named parameter
//     (one of 'projects' / 'samples' / 'experiments' / 'runs') because the
//     Reports API URL is kind-tagged and aliases are unique per-account-per-
//     kind.
//
// Within a session, an accession can also be recovered from an alias via
// `ena.submission_log.object_aliases` / `object_accessions` (see
// docs/ena.md) — useful when the user wants to inspect the accession before
// passing it back in.
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
