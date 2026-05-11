// SPDX-License-Identifier: MIT
//
// Webin Reports API client — alias → primary-accession lookup.
//
// Used by L5 to translate user-supplied aliases into server-assigned
// accessions before lifecycle envelope build, restoring the alias-based UX
// removed in L1d (decision #12). The envelope-builder rejection of
// `target_refname` is intentionally preserved as defense-in-depth; only the
// SQL-surface translates refname → accession before constructing the
// LifecycleSubmitOptions.
//
// Wire shape (live-verified 2026-05-08 via localdocs/ena-live-reports-probe.sh):
//   GET <base>/<kind>/{id}?format=json
//     Authorization: Basic <Webin user:password>
//   200  [{"report":{"id":"<primary-accession>","alias":"<alias>", …}}]   (hit)
//   200  []                                                               (miss)
// Both alias and primary accession resolve via the same {id} segment.
//
// AliasObjectKind → URL path:
//   STUDY      → /projects     (id = PRJEB…   — primary; lifecycle target)
//   SAMPLE     → /samples      (id = ERS…)
//   EXPERIMENT → /experiments  (id = ERX…)
//   RUN        → /runs         (id = ERR…)
// (Note: /studies also exists but returns id = ERP — the secondary
// accession; lifecycle targets the primary, so we don't use it.)
//
// All entry points are pure-data (no DuckDB linkage) and parameterised on a
// URLFetcher callback. Tests in test/cpp/test_ena_reports_client.cpp use a
// mock fetcher; the live wiring binds the fetcher to
// `ENAClient::AuthenticatedGet` at the call site (analogous to
// `ena_alias_check`'s integration in `ena_object_insert_op.hpp`).

#pragma once

#include "ena_alias_check.hpp" // AliasObjectKind, URLFetcher

#include <string>

namespace miint {

// Build the Reports lookup URL for a single (kind, id) probe. `id` accepts
// either an alias or a primary accession; the server resolves both. Throws
// std::invalid_argument on empty id or id containing newline characters
// (header-injection guard).
std::string BuildReportsLookupURL(const std::string &reports_base, AliasObjectKind kind, const std::string &id);

// One-call lookup: GET the single-id Reports endpoint, parse the JSON array,
// return `report.id` of the first element on hit, empty string on miss.
//
// Returns:
//   non-empty primary accession  on HIT (PRJEB / ERS / ERX / ERR)
//   ""                           on MISS (server returned `[]`)
// Throws:
//   std::invalid_argument        on empty alias (programmer error — table-fn
//                                bind layer must whitespace-guard before
//                                calling)
//   std::runtime_error           on transport failure (fetcher exception
//                                propagated), or malformed JSON, or hit
//                                without `report.id`
std::string LookupAccessionByAlias(const std::string &reports_base, AliasObjectKind kind, const std::string &alias,
                                   const URLFetcher &fetcher);

// Default Reports API bases. Reports is endpoint-segregated (each server
// only sees its own account's records), unlike the portal API used by
// ena_alias_check which is unified. The L5 client must hit the same
// endpoint the lifecycle POST will hit, otherwise an alias registered on
// wwwdev would surface as "not found" when the lookup goes to production.
constexpr const char *DEFAULT_REPORTS_BASE_PROD = "https://www.ebi.ac.uk/ena/submit/report";
constexpr const char *DEFAULT_REPORTS_BASE_TEST = "https://wwwdev.ebi.ac.uk/ena/submit/report";

// Resolve the Reports base used by callers: env-var override
// (MIINT_ENA_REPORTS_URL_BASE — test fixture only), falling back to the
// per-endpoint default. `endpoint` is the label from the ENA secret /
// ATTACH options ("production" → www, anything else / empty → wwwdev),
// matching ResolveDefaultENAEndpointURL's mapping for the submit URL.
// Trailing slashes stripped.
std::string ResolveReportsBaseForEndpoint(const std::string &endpoint);

} // namespace miint
