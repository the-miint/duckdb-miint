// SPDX-License-Identifier: MIT
//
// ENA Webin V2 JSON submission envelope builder.
//
// Builds the wire-format JSON document that gets POSTed to /ena/submit/webin-v2/submit
// (or /submit/queue). Pure data: takes a SubmissionSpec, returns a compact
// JSON string. No yyjson / no DuckDB dependency, so the unit-test target can
// link this file directly without pulling in the duckdb library.

#pragma once

#include <cstdint>
#include <string>
#include <utility>
#include <vector>

namespace miint {

enum class ENAAction { ADD, MODIFY, CANCEL, HOLD, RELEASE, VALIDATE };

struct ProjectSpec {
	std::string alias;
	std::string title;
	std::string description; // optional
	std::string project_type; // METAGENOMIC, WGS, ... — currently informational
	bool is_umbrella = false;
};

struct SampleSpec {
	std::string alias;
	std::string title;     // optional
	int64_t taxon_id = 0;  // required, > 0
	std::string checklist; // optional, ERC000NN — emitted as ENA-CHECKLIST attribute
	// Application-supplied attribute pairs (tag, value); preserves insertion order.
	std::vector<std::pair<std::string, std::string>> attributes;
};

struct SubmissionSpec {
	ENAAction action = ENAAction::ADD;
	std::string hold_until_date; // optional, "YYYY-MM-DD" or ISO-8601 with TZ
	std::vector<ProjectSpec> projects;
	std::vector<SampleSpec> samples;
	// Future phases append: experiments, runs, analyses
};

// Build the V2 JSON envelope. Compact (no whitespace) for byte-stable output.
// Throws std::runtime_error on invariant violations:
//   - empty project or sample alias
//   - sample taxon_id <= 0
//   - action == HOLD without hold_until_date
//   - hold_until_date set with action == HOLD (the date is inferred; do not
//     also set action=HOLD — pair `action=ADD, hold_until_date="..."` instead)
// Intent is to surface programmer errors before any HTTP traffic.
//
// String inputs (alias, title, description, attribute keys/values) are emitted
// as JSON strings with RFC 8259 §7 escaping. Inputs are treated as opaque
// bytes; valid UTF-8 is the caller's responsibility, since the function does
// not validate or substitute U+FFFD.
std::string BuildEnvelopeJSON(const SubmissionSpec &env);

} // namespace miint
