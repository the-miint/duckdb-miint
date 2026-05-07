// SPDX-License-Identifier: MIT
//
// Shared helpers for the ENA virtual-catalog INSERT operators
// (ENAProjectsInsert, ENASamplesInsert, future ENAExperimentsInsert/etc.).
//
// The chunk-buffering Sink and the RETURNING/count Source path are nearly
// identical across object types — the per-table differences are limited to
// the spec-builder, the submit functor, and the RETURNING column population.
// Routines here cover the parts that don't vary.

#pragma once

#include "duckdb/common/index_vector.hpp"
#include "duckdb/common/types/timestamp.hpp"
#include "duckdb/common/types/value.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

class ClientContext;
class ENACatalog;
struct ENASubmissionLogRow;

// User/password/endpoint resolved from the ENA secret named at ATTACH time.
// Consumed by every ENA insert operator.
struct ResolvedENACredentials {
	string user;
	string password;
	string endpoint_url;
	string endpoint;
	string secret_name;
};

ResolvedENACredentials ResolveENACredentials(ClientContext &context, ENACatalog &catalog);

// Look up a Webin secret by name and resolve its user/password/endpoint
// fields. Used by code paths that don't have an ENACatalog handle
// (notably the lifecycle table functions). `caller` is the user-facing
// function name, embedded in error messages so the user sees
// "ena_cancel: ..." rather than a generic "ENA secret: ...".
ResolvedENACredentials ResolveENACredentialsByName(ClientContext &context, const string &caller,
                                                   const string &secret_name);

// Translate a logical (table) column index into a position in the input
// chunk, honouring `LogicalInsert::column_index_map` semantics. Returns
// `DConstants::INVALID_INDEX` when the column was not provided.
idx_t ResolveInputColumn(const physical_index_vector_t<idx_t> &column_index_map, idx_t input_chunk_columns,
                         idx_t table_column_index);

// Stringify a value as VARCHAR; returns "" for NULL.
string ValueToVarchar(const Value &v);

// Stringify a DATE value as YYYY-MM-DD; returns "" for NULL. Throws if the
// value is not a DATE — surfaces type drift as a clean error rather than a
// silent format mismatch.
string ValueToDateString(const Value &v);

// Generate a fresh UUID v4 string suitable for submission_log.submission_id.
string GenerateSubmissionId();

// Per-row data needed to populate one submission_log entry. The object-type
// dispatcher (e.g. "projects" / "samples") is filled in by the caller.
struct SubmissionLogPayload {
	string object_type;
	string action; // "ADD", "MODIFY", "CANCEL", "HOLD", "RELEASE", "VALIDATE"
	int32_t n_objects;
	bool success;
	int64_t duration_ms;
	string envelope_payload;
	string raw_receipt;
	string era_accession;
	std::vector<std::string> error_messages;
	// Lifecycle target accession or refname; empty for body-style actions
	// (ADD / MODIFY / VALIDATE).
	string target;
	// Per-object alias / primary-accession parallel arrays. Populated by
	// the ADD path so users can later recover an accession from an alias
	// (lifecycle ops on already-registered objects need the accession;
	// see docs/ena.md). Empty on lifecycle ops.
	std::vector<std::string> object_aliases;
	std::vector<std::string> object_accessions;
};

void RecordSubmissionLog(ENACatalog &catalog, const ResolvedENACredentials &creds, const SubmissionLogPayload &payload);

// Look up the named ATTACHed database and return it iff its catalog is an
// ENACatalog. Behaviour depends on whether the catalog name was the
// (silent) default 'ena' or explicitly requested by the user via `catalog =>`:
//   - default name + missing/wrong-type → nullptr (silent skip; matches the
//     one-shot CANCEL/MODIFY UX where the user has no full ENA catalog
//     attached).
//   - explicit name + missing/wrong-type → throw InvalidInputException (the
//     user asked for audit logging; failing silently would create an
//     invisible audit gap).
// `caller` is the user-facing function name — embedded in error messages
// (e.g. "ena_modify_project: catalog 'foo' is not attached").
ENACatalog *FindAttachedENACatalog(ClientContext &context, const string &caller, const string &catalog_name,
                                   bool explicit_name);

// True iff `s` contains only ASCII whitespace (or is empty). Used by the
// table functions to bind-time-reject inputs like `accession => '   '`
// which pass non-empty but would emit garbage to the server.
bool IsENAStringWhitespaceOnly(const string &s);

} // namespace duckdb
