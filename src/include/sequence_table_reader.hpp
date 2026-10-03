#pragma once

#include "Minimap2Aligner.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"
#include <memory>
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

// SELECT yielding BuildQueryReadsSelect's columns for every row of
// `snapshot_table` that `routing_snapshot` assigns to `shard_name`, meant to be
// streamed (QuerySequenceStream's SelectSql form), once per shard and index part.
//
// BOTH relations are snapshots, and for the same #229 reason: this SELECT is
// re-run once per shard and part, so naming the user's relations here would read
// each of them many times over. `snapshot_table` is a MaterializeQueryReads
// snapshot and `routing_snapshot` a MaterializeReadToShard one.
//
// `snapshot_table` is the query relation itself, with no shard join — so a read
// routed to N shards is held once, not N times. The earlier shard-keyed snapshot
// of `query JOIN read_to_shard` was the other way round, and all-shards routing
// of 1M HiFi reads over 1000 shards made it ~1000 copies of the corpus before
// the first alignment.
//
// Duplicate (read_id, shard_name) rows in the routing relation align the read
// ONCE against that shard: `IN` is set membership, so a read listed twice for a
// shard is indistinguishable from one listed once. v1.0.0-v1.0.1 joined instead
// and emitted the alignment twice.
//
// Membership is an `IN (subquery)` in the SELECT list, deliberately not a JOIN or
// an IN in WHERE. DuckDB plans a projected IN as a MARK join, which it never
// flips, so the hash table is always built on this shard's read_ids and the
// snapshot (read sequences) is only ever probed. A SEMI join is flipped to
// RIGHT_SEMI whenever the routing side's estimate is larger — the all-shards
// case, where `shard_name = x` is estimated at 20% of a table holding reads x
// shards rows — and would then hash the whole corpus, sequences included, once
// per shard.
std::string BuildShardReadsSelect(const std::string &snapshot_table, const std::string &routing_snapshot,
                                  const SequenceTableSchema &schema, const std::string &shard_name);

// The projection every read of a sequence relation uses: exactly the columns
// `schema` says alignment consumes, in a fixed order. Shared by the snapshot
// builder and the streaming reader so the snapshot's columns and the columns
// later selected back out of it cannot drift apart.
std::string BuildQueryReadsSelect(const std::string &query_table, const SequenceTableSchema &schema);

// Materialize the query relation into a per-call TEMP table, reading it exactly
// ONCE (#229 — see docs/internals/reading-tables-views.md § "Read the relation
// ONCE"), for a consumer that needs to replay it more than once: align_minimap2
// streaming a multi-part prebuilt index (one pass per part), and
// align_minimap2_sharded (one pass per shard and part, each filtered to that
// shard by BuildShardReadsSelect).
//
// Returns the unquoted TEMP table name. Created on `conn`, which must inherit
// the caller's TEMP catalog. The caller MUST drop it via DropHelperTempRelation.
//
// Streams query_table and appends each chunk into the snapshot as it arrives
// (one pass, via SendQuery + Appender) rather than a single CREATE TABLE AS
// SELECT that pulls the whole query relation through the pipeline before this
// call returns — bounds the materialization's own working set to O(one chunk)
// instead of O(corpus size), which is what a full-corpus-sized query relation
// (tens of millions of reads) needs to fit in memory alongside a multi-part
// index's parts.
//
// `out_row_count` receives the number of rows materialized —
// summed from the streamed chunks as they're appended, rather than a second
// query. Lets a caller with zero query rows skip replaying every remaining
// index part for nothing (align_minimap2's multi-part path).
std::string MaterializeQueryReads(Connection &conn, const std::string &query_table, const SequenceTableSchema &schema,
                                  idx_t &out_row_count);

// The same one-pass TEMP snapshot, for align_minimap2_sharded's routing relation
// (`read_to_shard`). Holds (read_id, shard_name) only.
//
// Routing is read as many times as the reads are — once per shard and index part
// — so it needs the #229 guarantee just as much: a routing view built on
// `random()`, or one over a table another connection is writing, would otherwise
// send a read to one set of shards on part 1 and a different set on part 2, and
// the reads it stopped naming would simply never be aligned. Note that bind-time
// shard discovery (ReadShardNameCounts) reads the user's relation separately and
// before this, so the shard LIST can still come from a different observation than
// the membership; this snapshot is what makes membership stable for the scan.
//
// Returns the unquoted TEMP table name; the caller MUST drop it via
// DropHelperTempRelation.
std::string MaterializeReadToShard(Connection &conn, const std::string &read_to_shard_table, idx_t &out_row_count);

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

	// A full SELECT to stream instead of a table name. It must yield
	// BuildQueryReadsSelect's columns for `schema`, in that order (e.g.
	// BuildShardReadsSelect). Wrapped in its own type so it cannot be passed where
	// a table name is expected, or the reverse.
	//
	// `label` names the stream in error messages. Without it the message quotes
	// the generated SQL, which for a shard read is a two-UUID-table MARK-join
	// SELECT — accurate and unreadable. Say "reads for shard 'x'" instead.
	struct SelectSql {
		std::string sql;
		std::string label;
	};

	// Owns an internal Connection that inherits the caller's TEMP objects.
	//
	// `source_keepalive`, if given, is held for the life of the stream. A SELECT
	// over a caller-owned TEMP table must not outlive that table, and the stream
	// is the thing that actually reads it — so it pins it, rather than the
	// caller maintaining a destruction order it cannot enforce.
	QuerySequenceStream(ClientContext &context, const SelectSql &select, const SequenceTableSchema &schema,
	                    idx_t sub_batch_size = STANDARD_VECTOR_SIZE, std::shared_ptr<void> source_keepalive = nullptr);

	// Closes the stream before any member (the owned connection, the keepalive)
	// is destroyed, so neither the connection nor the source table can go away
	// underneath a live StreamQueryResult. Doing this in the body rather than by
	// member declaration order means reordering the members cannot break it.
	~QuerySequenceStream();

	// Fetch the next sub-batch. Returns an empty batch when the stream is exhausted.
	// Thread-safe — serializes access to the underlying stream via mutex.
	// Throws if the query fails mid-stream: DuckDB reports that as the same null
	// chunk as end-of-input, so without the check a failure would be silently
	// truncated rows.
	miint::SequenceRecordBatch FetchSubBatch();

	// Rows returned by FetchSubBatch so far. Exact once FetchSubBatch has
	// returned an empty batch: every earlier batch was counted under the same
	// mutex before it was handed out.
	idx_t RowsDelivered() const;

private:
	// Both ClientContext& forms above are this one: create an owned connection,
	// inherit the caller's TEMP objects, stream `select_sql`. They differ only in
	// how they build the SELECT and how they name it in errors, so the setup
	// lives here once. `source` names what is being read, for error messages.
	QuerySequenceStream(ClientContext &context, const std::string &select_sql, const std::string &source,
	                    const SequenceTableSchema &schema, idx_t sub_batch_size,
	                    std::shared_ptr<void> source_keepalive = nullptr);

	// Only one of these two is populated: owned_conn_ for the ClientContext&
	// constructor, nullptr for the Connection& constructor. conn_ptr_ always
	// points to the live connection.
	unique_ptr<Connection> owned_conn_;
	Connection *conn_ptr_;
	// Whatever the SELECT reads from and does not own — see the SelectSql ctor.
	std::shared_ptr<void> source_keepalive_;
	unique_ptr<QueryResult> stream_;
	SequenceTableSchema schema_;
	idx_t sub_batch_size_;
	miint::SequenceRecordBatch partial_; // Partially-filled sub-batch carried across Fetch() calls
	bool exhausted_ = false;
	idx_t rows_delivered_ = 0;
	mutable std::mutex mutex_;
	// Reusable temp vectors for string extraction (avoids per-chunk allocation)
	std::vector<std::string> temp_read_ids_;
	std::vector<std::string> temp_seq1_;
	std::vector<std::string> temp_seq2_;

	// `source` names what is being read, for error messages only.
	void InitStream(const std::string &select_sql, const std::string &source);
};

} // namespace duckdb
