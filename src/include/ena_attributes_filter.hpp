#pragma once

#include "ena_parser.hpp"

#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace miint {

// Abstract filter-tree node used by the `read_ena_attributes` pushdown
// extractor. Kept abstract (independent of `duckdb::Expression`) so the
// extractor can be unit-tested without constructing DuckDB bound
// expressions. The DuckDB→ENAFilterNode translator lives alongside the
// `pushdown_complex_filter` hook in `read_ena_attributes.cpp`.
//
// Policy: the extractor is intentionally strict. Any node that cannot be
// cleanly translated to an ENA `/search?query=...` predicate — unsupported
// operator, OR across atoms, non-`tag`/`value` column, tag not in the
// `ENASearchFieldRegistry` allowlist — forces the whole filter set to fall
// back to the XML path. The structured path only runs when we can prove it
// is equivalent to the XML output.
struct ENAFilterNode {
	enum class Kind { EQUAL, IN, AND_CONJUNCTION, OR_CONJUNCTION, UNSUPPORTED };

	Kind kind;
	// For EQUAL / IN: lowercase column name (the translator lowercases).
	std::string column;
	// For EQUAL: single value. For IN: one per element. Exact user-typed case
	// is preserved; canonicalization (e.g., lowercasing `tag` values against
	// the registry) happens inside ExtractPushdownPredicates.
	std::vector<std::string> values;
	// For AND_CONJUNCTION / OR_CONJUNCTION.
	std::vector<std::unique_ptr<ENAFilterNode>> children;

	static std::unique_ptr<ENAFilterNode> MakeEqual(const std::string &column, const std::string &value);
	static std::unique_ptr<ENAFilterNode> MakeIn(const std::string &column, const std::vector<std::string> &values);
	static std::unique_ptr<ENAFilterNode> MakeAnd(std::vector<std::unique_ptr<ENAFilterNode>> children);
	static std::unique_ptr<ENAFilterNode> MakeOr(std::vector<std::unique_ptr<ENAFilterNode>> children);
	static std::unique_ptr<ENAFilterNode> MakeUnsupported();
};

// Result of analyzing a filter conjunction. `tags.empty()` is the signal that
// pushdown is not possible — the caller should use the XML path.
//
// When `tags` is non-empty, every entry is a canonical lowercase ENA sample
// field name, and every entry of `tag_value_pairs` is `(tag, value)` with
// `tag` guaranteed to appear in `tags`.
struct ENAAttributePushdown {
	std::vector<std::string> tags;
	std::vector<std::pair<std::string, std::string>> tag_value_pairs;
};

// Analyze a conjunct list (as would be produced by DuckDB's
// `pushdown_complex_filter` callback, which pre-splits top-level AND).
// Returns populated `tags` iff EVERY conjunct is pushable AND every
// referenced tag passes `ENASearchFieldRegistry::IsSearchableSampleField`.
// Otherwise returns an empty pushdown (XML fallback).
ENAAttributePushdown ExtractPushdownPredicates(const std::vector<std::unique_ptr<ENAFilterNode>> &conjuncts);

// Build the ENA `/search?result=sample` URL for a batched sample-accession
// list plus a pushdown's tag+value constraints. `sample_accessions` must be
// non-empty. Emits `fields=sample_accession,<tag1>,<tag2>,...` (deduped),
// and a `query` clause joining the accession IN-filter with any pinned
// tag=value equality constraints via URL-encoded AND.
//
// Preconditions:
//   - `sample_accessions` non-empty; each entry passes ENAParser::ValidateAccession.
//   - `tags` non-empty; each entry passes ENAParser::ValidateFields.
//   - Every `tag_value_pairs[i].first` appears in `tags` (so the returned
//     TSV has a column to un-pivot against). Violations throw
//     `std::invalid_argument`.
//
// Pure (no HTTP). Kept here so the pushdown path is testable without a
// network dependency.
std::string BuildStructuredSearchURL(const std::vector<std::string> &sample_accessions,
                                     const std::vector<std::string> &tags,
                                     const std::vector<std::pair<std::string, std::string>> &tag_value_pairs);

// Unpivot a structured-search TSV response into `(sample_accession, tag,
// value)` tuples — one row per (sample, non-empty requested tag). Expects
// `parsed` to contain a `sample_accession` column plus one column per entry
// in `requested_tags` (case-insensitive column match). Missing columns are
// treated as all-empty (no rows emitted for that tag). Empty cells are
// skipped (matches the XML path, which omits absent attributes).
std::vector<SampleAttribute> UnpivotStructuredTSV(const ENATSVResult &parsed,
                                                  const std::vector<std::string> &requested_tags);

} // namespace miint
