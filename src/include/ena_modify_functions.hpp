// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 MODIFY actions on submittable objects.
//   - ena_modify_project(secret, accession, alias, title, description?,
//                         project_type?, is_umbrella?, catalog?)
//   - ena_modify_sample(secret, accession, alias, taxon_id, scientific_name?,
//                        title?, description?, checklist?, attributes?,
//                        attribute_units?, catalog?)
//   - ena_modify_experiment / ena_modify_run — future phases (L4c, L4d).
//
// MODIFY semantics: re-submit the full updated body for an already-registered
// object, identified by its server-assigned `accession` (PRJEB / ERS / ERX /
// ERR). Webin V2 expresses updates via `<ACTION><MODIFY/></ACTION>` rather
// than REST PUT/PATCH; the alias still labels the local reference and is
// emitted alongside the accession on the per-object element.
//
// Each function returns a single row:
//   action, target, success, alias, era_accession, error_messages, duration_ms
// where `target` is the user-supplied accession echoed back (so the row
// confirms which object the server acted on).
//
// Submission_log integration: identical to the lifecycle table functions
// (default catalog 'ena' silently no-ops if not attached; explicit `catalog =>`
// throws on missing/wrong-type). Audit row records `action='MODIFY'`,
// `target=<accession>`, with `object_aliases`/`object_accessions` populated
// from the receipt's per-object children.

#pragma once

namespace duckdb {
class ExtensionLoader;
}

namespace miint {

// Register every shipped ena_modify_* table function with the loader.
// Today wires ena_modify_project + ena_modify_sample; experiments and runs
// land in subsequent phases. The registration site is one function so
// miint_extension.cpp doesn't grow a new line per object kind.
void RegisterENAModifyTableFunctions(duckdb::ExtensionLoader &loader);

} // namespace miint
