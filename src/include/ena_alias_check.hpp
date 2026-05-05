// SPDX-License-Identifier: MIT
//
// Pre-INSERT alias collision check (Phase 8 Step 8a).
//
// Aliases are unique per (submission_account, object_type) on the ENA side.
// Reusing an alias for the same object type is a hard server error. Catching
// this client-side via a portal-API search round-trip surfaces the conflict
// before any envelope POST and gives the user the same alias name in the
// error message instead of the receipt's free-form string.
//
// All entry points are pure-data (no DuckDB linkage) and parameterised on a
// URLFetcher callback, so the unit tests in test/cpp/test_ena_alias_check.cpp
// substitute a mock for live HTTP. The DuckDB-side glue lives in the
// per-table operator (ena_object_insert_op.hpp PreSubmit hook).

#pragma once

#include <functional>
#include <string>
#include <vector>

namespace miint {

// Object kinds the ENA portal API recognizes. Maps to result=read_<kind> and
// fields=<kind>_alias on the search endpoint.
enum class AliasObjectKind {
	STUDY,      // ena.projects
	SAMPLE,     // ena.samples
	EXPERIMENT, // ena.experiments
	RUN,        // ena.runs
};

// Lowercase singular token used in URL build and error messages.
const char *AliasObjectKindName(AliasObjectKind kind);

// Map ena.<table_name> → AliasObjectKind. Throws std::invalid_argument on
// unknown table names (so a future ena.analyses gets a loud failure instead
// of silently bypassing the check).
AliasObjectKind AliasObjectKindFromTableName(const std::string &table_name);

// Free-function URL fetcher for testability. Implementations must perform an
// authenticated GET — anonymous portal search would not see the user's own
// HOLD/private records. See ENAClient::AuthenticatedGet.
using URLFetcher = std::function<std::string(const std::string &url)>;

// Build the portal API URL for a single alias-collision query batch.
// `portal_base` should be the portal API base (no trailing slash), e.g.
// "https://www.ebi.ac.uk/ena/portal/api". `aliases` must be non-empty,
// each non-empty, and free of `"` and newline characters; otherwise throws
// std::invalid_argument.
std::string BuildAliasCollisionURL(const std::string &portal_base, AliasObjectKind kind,
                                   const std::vector<std::string> &aliases);

// Returns the subset of `aliases` that already exist in the user's submission
// account, in the same order as the input. Empty input → empty output, no
// fetch invoked. Calls the fetcher one or more times (chunking is automatic
// when alias count exceeds an internal per-URL limit).
//
// Aliases that appear twice in the input collapse to one fetch entry and
// produce at most one hit. Aliases the server names but we did not query for
// are ignored.
std::vector<std::string> CheckAliasCollisions(const std::string &portal_base, AliasObjectKind kind,
                                              const std::vector<std::string> &aliases, const URLFetcher &fetcher);

// Default portal API base. May be overridden via the MIINT_ENA_PORTAL_URL_BASE
// env var (test fixture only — production users should not need to set it).
constexpr const char *DEFAULT_PORTAL_BASE = "https://www.ebi.ac.uk/ena/portal/api";

// Resolve the portal base used by the operator: env-var override, falling back
// to DEFAULT_PORTAL_BASE. Trailing slashes stripped.
std::string ResolvePortalBaseFromEnv();

} // namespace miint
