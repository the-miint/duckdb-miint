#pragma once

#include "duckdb/main/client_context.hpp"
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace duckdb {

// Read a reference table or view and return map of reference name -> length.
// Used by read_alignments and COPY FORMAT SAM/BAM for headerless SAM support.
//
// Both tables and views are supported. The implementation uses a separate Connection
// to execute a SELECT statement, which DuckDB's binder resolves appropriately for
// either tables or views.
//
// Required columns (by position):
//   - Column 0: reference name (VARCHAR) - must not be NULL
//   - Column 1: reference length (BIGINT, INTEGER, UBIGINT, or UINTEGER) - must not be NULL
//
// Throws InvalidInputException if:
//   - Table/view doesn't exist
//   - Table/view has fewer than 2 columns
//   - Column types are incorrect
//   - Table/view is empty
//   - NULL value in name or length column
//   - Negative length value
//   - UBIGINT length exceeds INT64_MAX
std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name);

// Same validation and error contract as ReadReferenceTable, but returns (name, length) pairs
// SORTED BY REFERENCE NAME instead of a hash map. Use this wherever the ORDER of references is
// observable -- specifically SAM/BAM @SQ header emission, where the index of a reference in the
// header IS its TID and therefore defines what "coordinate-sorted" means.
//
// This exists because @SQ was previously emitted by iterating the unordered_map above, i.e. in
// hash-bucket order: deterministic, but equal to neither the table's row order, nor its reverse,
// nor name order. No ORDER BY the caller could write produced a coordinate-sorted BAM, and
// htslib-based consumers rejected the output (issue #173).
//
// Name order (not table row order) is the contract deliberately: DuckDB does not guarantee the
// row order of a parallel table scan, so row order could not be reproducible. It also makes the
// natural `ORDER BY reference, position` actually yield a coordinate-sorted file -- which
// depends on this sort agreeing with DuckDB's VARCHAR comparison. Both are plain byte-wise
// comparisons (DuckDB's default collation is binary), so they agree.
//
// Karyotype order (chr1..chr9, chr10) is NOT expressible this way, since 'chr10' sorts before
// 'chr2'. That needs an explicit ordinal column and is a recorded follow-up on #173.
//
// Duplicate names are collapsed keeping the FIRST occurrence, matching ReadReferenceTable's
// map-insert behaviour (and required for valid SAM, which cannot repeat an @SQ SN).
std::vector<std::pair<std::string, uint64_t>> ReadReferenceTableSortedByName(ClientContext &context,
                                                                             const std::string &table_name);

} // namespace duckdb
