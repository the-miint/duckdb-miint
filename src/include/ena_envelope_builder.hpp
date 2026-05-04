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
	std::string description;  // optional
	std::string project_type; // METAGENOMIC, WGS, ... — currently informational
	bool is_umbrella = false;
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
std::string BuildEnvelopeJSON(const SubmissionSpec &env);

} // namespace miint
