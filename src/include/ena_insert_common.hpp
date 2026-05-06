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
	string action; // "ADD", "MODIFY", ...
	int32_t n_objects;
	bool success;
	int64_t duration_ms;
	string envelope_payload;
	string raw_receipt;
	string era_accession;
	std::vector<std::string> error_messages;
};

void RecordSubmissionLog(ENACatalog &catalog, const ResolvedENACredentials &creds, const SubmissionLogPayload &payload);

} // namespace duckdb
