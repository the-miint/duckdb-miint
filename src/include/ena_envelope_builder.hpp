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

// Wire-format action name ("ADD" / "MODIFY" / ...). Lifted out of the
// envelope builder's anonymous namespace so submission_log writers can
// stringify the same way the wire format does, without duplicating the
// switch.
const char *ActionName(ENAAction a);

// Shared across all four per-table insert paths. Each `Submit*InsertOutcome`
// takes this; per-table aliases (ENAProjectInsertOptions, ENASampleInsertOptions,
// ...) are kept as type aliases for call-site readability.
struct ENAInsertOptions {
	std::string endpoint_url; // resolved base URL incl. /submit suffix
	std::string user;         // Webin-XXXXX
	std::string password;
	std::string hold_until_date; // optional, "YYYY-MM-DD"
	// Wire-format action. ADD is the registration path (assigns accessions);
	// VALIDATE is a server-side dry-run that returns a receipt with no
	// per-object children. MODIFY will follow the same envelope shape as ADD
	// once L4 ships. Lifecycle ops (CANCEL/RELEASE/HOLD) go through
	// SubmitLifecycle and never touch this struct.
	ENAAction action = ENAAction::ADD;
};

// Common shape of every Submit*Outcome. The four per-table outcomes
// (ENASubmissionOutcome / ENASamplesSubmissionOutcome / ... ) are aliases of
// this template specialised on the per-table Result row type. Fields:
//   - rows: server-assigned per-row outputs (alias + accession + ext_ids)
//   - envelope_payload: the request body we POSTed (passwords NOT embedded)
//   - raw_receipt: the raw response body
//   - era_accession: server-assigned <SUBMISSION accession>
//   - success: false if logical failure (server success=false OR parse error)
//   - error_messages: free-form server messages on failure
//   - duration_ms: wall-clock time of the POST + parse round-trip
template <class RowT>
struct ENABaseSubmissionOutcome {
	std::vector<RowT> rows;
	std::string envelope_payload;
	std::string raw_receipt;
	std::string era_accession;
	bool success = false;
	std::vector<std::string> error_messages;
	int64_t duration_ms = 0;
};

struct ProjectSpec {
	std::string alias;
	std::string title;
	std::string description;  // optional
	std::string project_type; // METAGENOMIC, WGS, ... — currently informational
	bool is_umbrella = false;
	// Existing PRJEB accession — required on MODIFY (Webin V2 needs both
	// `alias` and `accession` on the project element to identify the
	// already-registered object), ignored on ADD. The envelope builder emits
	// `"accession":"…"` only when this is non-empty; an ADD with this set is
	// silently passed through (server rejects the contradiction).
	std::string accession;
};

struct SampleSpec {
	std::string alias;
	std::string title;           // optional
	std::string description;     // optional, sample-level description
	int64_t taxon_id = 0;        // required, > 0
	std::string scientific_name; // optional, organism.scientificName
	std::string checklist;       // optional, ERC000NN — emitted as ENA-CHECKLIST attribute
	// Application-supplied attribute pairs (tag, value); preserves insertion order.
	std::vector<std::pair<std::string, std::string>> attributes;
	// Optional per-attribute units (sparse — only entries that need them).
	// Some checklist attributes mandate a `<UNITS>` sibling — e.g. ERC000015
	// rejects `geographic location (latitude)` without `unit: DD`. Lookup is
	// by tag (matches against `attributes[].first`); tags not present here
	// are emitted without a `units` JSON field.
	std::vector<std::pair<std::string, std::string>> attribute_units;
};

// Cross-reference to a parent object (study / sample / experiment). V2 accepts
// either an accession (the server-assigned ID, e.g. PRJEB123) or a refname
// (the parent's alias, useful for first-time bulk submissions where parent and
// children are POSTed together). `accession` wins when both are set.
struct RefDescriptor {
	std::string accession;
	std::string refname;
};

enum class ENALibraryLayout { SINGLE, PAIRED };

// One file referenced by a Run's <DATA_BLOCK><FILES>. The server re-computes
// MD5 after upload and compares against `checksum`; mismatch → validation
// error. `filetype` is one of the SRA.run.xsd enum (see
// localdocs/ena-research-webin-v2-deep.md §6.1) — `fastq` is the common case.
struct RunFile {
	std::string filename;
	std::string filetype; // "fastq", "bam", "cram", ...
	std::string checksum; // MD5 hex (32 chars); SHA-256 also accepted by run XSD
};

struct ExperimentSpec {
	std::string alias;
	std::string title; // optional
	RefDescriptor study_ref;
	RefDescriptor sample_ref;
	std::string design_description; // optional
	std::string library_name;       // optional
	std::string library_strategy;   // required, e.g. WGS, AMPLICON, RNA-Seq
	std::string library_source;     // required, e.g. METAGENOMIC, GENOMIC
	std::string library_selection;  // required, e.g. RANDOM, PCR
	ENALibraryLayout library_layout = ENALibraryLayout::SINGLE;
	std::string platform;         // required, e.g. ILLUMINA, OXFORD_NANOPORE
	std::string instrument_model; // required, free-form (server validates)
};

struct RunSpec {
	std::string alias;
	std::string title; // optional
	RefDescriptor experiment_ref;
	std::vector<RunFile> files; // required, ≥ 1
};

struct SubmissionSpec {
	ENAAction action = ENAAction::ADD;
	std::string hold_until_date; // optional, "YYYY-MM-DD" or ISO-8601 with TZ
	// Targeted lifecycle actions (CANCEL, RELEASE, HOLD-with-date) reference an
	// already-registered object via `target=` on the action element. Body sets
	// (PROJECT_SET / SAMPLE_SET / ...) are not emitted when a target is set —
	// the action itself is the entire payload. `target_accession` wins when
	// both are populated, mirroring `RefDescriptor`. Setting either on
	// ADD/MODIFY/VALIDATE is rejected at validate time.
	std::string target_accession;
	std::string target_refname;
	std::vector<ProjectSpec> projects;
	std::vector<SampleSpec> samples;
	std::vector<ExperimentSpec> experiments;
	std::vector<RunSpec> runs;
	// Future phase: analyses
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
//
// V2 server caveat: this JSON envelope is only accepted for `projects` and
// `samples`. Submitting `experiments` / `runs` via JSON returns HTTP 500
// from a generic NPE in the V2 dispatcher (verified live 2026-05-04). Use
// `BuildEnvelopeXML` for those object types.
std::string BuildEnvelopeJSON(const SubmissionSpec &env);

// Build the V2 XML envelope. Currently scoped to experiments + runs (the
// SRA-side objects whose JSON dispatch is broken on V2). Same `SubmissionSpec`
// input shape as `BuildEnvelopeJSON`; populated `experiments` / `runs` arrays
// are emitted under `<EXPERIMENT_SET>` / `<RUN_SET>`. `projects` / `samples`
// in the spec are ignored here — submit them via `BuildEnvelopeJSON` instead.
//
// Compact, no whitespace, with the `<?xml ... ?>` declaration. String values
// are XML-escaped (`< > & " '`). Invariant violations throw
// std::runtime_error in the same shape as the JSON builder.
std::string BuildEnvelopeXML(const SubmissionSpec &env);

} // namespace miint
