#pragma once

#include "Minimap2Aligner.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"
#include <mutex>
#include <string>
#include <vector>

namespace duckdb {

// Schema info for a sequence table/view
struct SequenceTableSchema {
	bool has_sequence2 = false;    // True if paired-end (sequence2 column exists)
	bool has_qual1 = false;        // True if quality scores present
	bool has_qual2 = false;        // True if quality scores present for second read
	bool is_physical_table = true; // True if physical table (has rowid), false if view
	// Storage type of the read_id column (VARCHAR or BIGINT). Defaults to
	// INVALID so any code path that constructs SequenceTableSchema without
	// going through ValidateSequenceTableSchema fails loud at the helper
	// dispatch in id_column_utils.hpp rather than silently misreading data.
	LogicalType id_type = LogicalType(LogicalTypeId::INVALID);
};

// Validate that a table/view has required columns for sequence data.
// Returns schema information about what optional columns are present.
// Throws BinderException if required columns are missing or have wrong types.
// `allow_bigint`: if true, accepts both VARCHAR and BIGINT for `read_id` and
// records the discovered type on `schema.id_type`. If false (default), only
// VARCHAR is accepted — preserves backward-compatible behavior for callers
// that haven't yet been audited for BIGINT support.
SequenceTableSchema ValidateSequenceTableSchema(ClientContext &context, const std::string &table_name,
                                                bool allow_bigint = false);

// Read all subjects from a table/view into memory.
// Subjects cannot be paired-end (sequence2 must be NULL for all rows).
// Throws InvalidInputException if sequence2 contains non-NULL values.
// The schema (in particular `id_type`) governs how the read_id column is
// extracted — VARCHAR rows pass through as strings; BIGINT rows are
// stringified into the carrier so the aligner sees a uniform string contract.
std::vector<miint::AlignmentSubject> ReadSubjectTable(ClientContext &context, const std::string &table_name,
                                                      const SequenceTableSchema &schema);

// NOTE: there is deliberately no offset/LIMIT batch reader here. Paging a
// relation with a fresh `ORDER BY ... LIMIT n OFFSET k` query per batch silently
// corrupts results for any relation that is not stable across re-evaluation
// (volatile views, views over a changing table, registered Arrow streams): the
// next page re-evaluates the relation and returns different rows, and an empty
// page is indistinguishable from end-of-input (#229). Use QuerySequenceStream
// below for a single streaming pass instead.

// SELECT yielding (shard_name, read_id, sequence1[, sequence2][, qual1][, qual2])
// for every read that `read_to_shard` assigns to a shard — one pass over the
// query relation. Shared so the snapshot below and any direct read of the same
// rows cannot drift apart.
std::string BuildShardedQueryReadsSelect(const std::string &query_table, const std::string &read_to_shard_table,
                                         const SequenceTableSchema &schema);

// Materialize the shard-assigned reads into a per-call TEMP table, reading the
// query relation exactly ONCE (#229 — see
// docs/internals/reading-tables-views.md § "Read the relation ONCE").
//
// Buffering is intrinsic for a sharded consumer rather than a choice: a read
// assigned to shard N cannot be aligned until shard N's index is loaded, and
// reads for all shards arrive interleaved, so a single-pass relation must be
// either buffered or re-read.
//
// Returns the unquoted TEMP table name. Created on `conn`, which must inherit the
// caller's TEMP catalog — that is what makes the table visible to the per-worker
// connections that also inherit it. The name is uniquified per call and the
// caller MUST drop it via DropHelperTempRelation; both halves are mandatory (see
// catalog_utils.hpp).
//
// Worth calling only for more than one shard: with a single shard the relation is
// already read once, so a snapshot would add a write and a scan for nothing.
std::string MaterializeShardedQueryReads(Connection &conn, const std::string &query_table,
                                         const std::string &read_to_shard_table, const SequenceTableSchema &schema);

// Read every read assigned to `shard_name`. `source_sql` is anything that can
// follow FROM and exposes a shard_name column — either a quoted snapshot table
// name from MaterializeShardedQueryReads, or a parenthesized
// BuildShardedQueryReadsSelect for the single-shard case where no snapshot is
// built. One query either way: no id list, no per-shard temp table, no appender.
void ReadShardReadsFrom(ClientContext &context, const std::string &source_sql, const SequenceTableSchema &schema,
                        const std::string &shard_name, miint::SequenceRecordBatch &output);

// Labels and single-end sequences loaded from a table for vsearch operations.
struct LoadedSingleEndSequences {
	std::vector<std::string> labels;
	std::vector<std::string> sequences;
};

// Load all (read_id, sequence1) pairs from a table via a separate connection.
// Used by vsearch-backed table functions (search, chimera, cluster) that operate
// on single-end sequences only.
// If strict=true: throws on NULL read_id, NULL sequence1, or empty sequence1.
// If strict=false: silently skips those rows.
// Always throws if the result set is empty after filtering.
// function_name is used in error messages (e.g. "cluster_sequences").
LoadedSingleEndSequences LoadSingleEndSequences(ClientContext &context, const std::string &table_name,
                                                const std::string &function_name, bool strict = false);

// Overload that runs against a caller-owned connection with an optional WHERE clause.
// Useful for per-sample callers: they already hold a per-thread Connection in LocalState
// and want to filter by a sample predicate. Pass `where_sql` without the leading "WHERE".
// An empty `where_sql` is equivalent to the context-based overload.
LoadedSingleEndSequences LoadSingleEndSequences(Connection &conn, const std::string &table_name,
                                                const std::string &function_name, bool strict,
                                                const std::string &where_sql);

// Streaming query sequence reader for lazy sub-batching.
// Produces sub-batches on demand from a streaming query result.
// Thread-safe: multiple threads can call FetchSubBatch() concurrently.
class QuerySequenceStream {
public:
	// Owns an internal Connection — convenient for single-pipeline readers.
	QuerySequenceStream(ClientContext &context, const std::string &table_name, const SequenceTableSchema &schema,
	                    idx_t sub_batch_size = STANDARD_VECTOR_SIZE);

	// Uses a caller-owned Connection. Required when the stream needs to see TEMP
	// objects (e.g. views) created by the caller on the same connection — per-sample
	// callers that do CREATE OR REPLACE TEMP VIEW … before streaming must use this.
	QuerySequenceStream(Connection &conn, const std::string &table_name, const SequenceTableSchema &schema,
	                    idx_t sub_batch_size = STANDARD_VECTOR_SIZE);

	// Fetch the next sub-batch. Returns an empty batch when the stream is exhausted.
	// Thread-safe — serializes access to the underlying stream via mutex.
	miint::SequenceRecordBatch FetchSubBatch();

private:
	// Only one of these two is populated: owned_conn_ for the ClientContext&
	// constructor, nullptr for the Connection& constructor. conn_ptr_ always
	// points to the live connection.
	unique_ptr<Connection> owned_conn_;
	Connection *conn_ptr_;
	unique_ptr<QueryResult> stream_;
	SequenceTableSchema schema_;
	idx_t sub_batch_size_;
	miint::SequenceRecordBatch partial_; // Partially-filled sub-batch carried across Fetch() calls
	bool exhausted_ = false;
	mutable std::mutex mutex_;
	// Reusable temp vectors for string extraction (avoids per-chunk allocation)
	std::vector<std::string> temp_read_ids_;
	std::vector<std::string> temp_seq1_;
	std::vector<std::string> temp_seq2_;

	void InitStream(const std::string &table_name);
};

} // namespace duckdb
