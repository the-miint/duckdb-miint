// SPDX-License-Identifier: MIT
//
// ENA checklist registry + sample-attribute validator (Phase 8 Step 8b).
//
// ENA samples target a checklist (ERC000015 = "GSC MIxS human gut", etc.),
// which lists the field definitions with mandatory/recommended/optional
// status, allowed units, and controlled vocabularies. Validating
// user-supplied attributes against the checklist client-side surfaces
// missing mandatory fields, missing units, and out-of-CV values BEFORE
// the envelope POST. Three of the four bugs caught during the Phase 7.5
// live-wwwdev pass would have been blocked here.
//
// The XML schema (https://www.ebi.ac.uk/ena/browser/api/xml/ERC000NN) keys
// fields by `<LABEL>` (human-readable, e.g. "project name") which is what
// users supply in the SAMPLE_ATTRIBUTE/TAG and in the duckdb-miint
// `attributes` MAP. `<NAME>` is the snake_case form, kept here for
// diagnostics only.

#pragma once

#include <functional>
#include <mutex>
#include <string>
#include <unordered_map>
#include <vector>

namespace miint {

enum class ChecklistFieldMandatory {
	MANDATORY,
	RECOMMENDED,
	OPTIONAL,
};

struct ChecklistFieldDef {
	std::string label; // human-readable; matches user attribute MAP keys
	std::string name;  // snake_case; informational
	ChecklistFieldMandatory mandatory = ChecklistFieldMandatory::OPTIONAL;
	std::vector<std::string> allowed_units;  // empty if no <UNITS> block
	std::vector<std::string> allowed_values; // empty if no <TEXT_CHOICE_FIELD> CV
};

struct ChecklistDef {
	std::string accession; // e.g. "ERC000015"
	std::string label;     // descriptor LABEL, e.g. "GSC MIxS human gut"
	std::vector<ChecklistFieldDef> fields;
};

// Parse a full <CHECKLIST_SET><CHECKLIST/></CHECKLIST_SET> XML body. Throws
// std::runtime_error on malformed XML or missing required elements.
ChecklistDef ParseChecklistXML(const std::string &xml);

// One human-readable issue surfaced by ValidateAttributesAgainstChecklist.
struct ChecklistValidationIssue {
	std::string field;   // field LABEL the issue applies to (or unknown user-supplied key)
	std::string message; // user-facing description of the issue
};

// Validate user-supplied attributes (label → value) and units (label → unit)
// against the parsed checklist. Returns a list of issues; empty vector means
// the input is valid by Phase 8b's bare-minimum rules:
//   1. Every MANDATORY field must be present in `attributes` with a non-empty
//      value.
//   2. Any user attribute key not in the checklist is flagged ("not in
//      checklist"). Reject unknown keys per Phase 8 scope decision.
//   3. Any user attribute whose field declares <UNITS> must have a
//      corresponding `units` entry whose value is in `allowed_units`.
//   4. Any user attribute whose field declares <TEXT_CHOICE_FIELD> must
//      have a value present in `allowed_values` (when non-empty).
std::vector<ChecklistValidationIssue>
ValidateAttributesAgainstChecklist(const ChecklistDef &checklist,
                                   const std::vector<std::pair<std::string, std::string>> &attributes,
                                   const std::vector<std::pair<std::string, std::string>> &units);

// URL builder for the ENA browser API. `base` should not have a trailing
// slash but is tolerated either way. Throws std::invalid_argument on empty
// or non-alphanumeric accession (defends against path traversal / URL
// injection if the value ever flows from user input).
std::string BuildChecklistFetchURL(const std::string &base, const std::string &accession);

// Default browser-API base. May be overridden via MIINT_ENA_CHECKLIST_URL_BASE
// (test fixture only).
constexpr const char *DEFAULT_CHECKLIST_URL_BASE = "https://www.ebi.ac.uk/ena/browser/api/xml";

// Resolve the browser-API base. Reads the env var override on every call;
// trailing slash stripped.
std::string ResolveChecklistBaseFromEnv();

// Process-lifetime cache. Lookup invokes the fetcher only on first miss;
// subsequent lookups for the same accession return the cached parsed
// ChecklistDef by reference. Fetcher exceptions propagate and do NOT poison
// the cache (a failed fetch is retried on the next call).
class ChecklistRegistry {
public:
	using Fetcher = std::function<std::string(const std::string &url)>;

	// Get a process-wide singleton. Production callers use this. Tests may
	// also instantiate ChecklistRegistry directly to keep state isolated.
	static ChecklistRegistry &Instance();

	// Returns a const reference to the parsed ChecklistDef for `accession`.
	// References to cached values stay valid as long as the corresponding
	// entry is not erased — std::unordered_map's rehash invalidates
	// iterators but NOT references to mapped values (per cppreference). We
	// never erase from the cache during normal operation; Clear() is
	// test-only and callers must not hold references across a Clear().
	const ChecklistDef &GetOrFetch(const std::string &accession, const Fetcher &fetcher);

	// Drop all cached checklists. **Test-only.** Invalidates every
	// reference previously returned by GetOrFetch. Production code must
	// not call this — references handed out across the operator path
	// would dangle.
	void Clear();

private:
	std::mutex mu;
	std::unordered_map<std::string, ChecklistDef> cache;
};

} // namespace miint
