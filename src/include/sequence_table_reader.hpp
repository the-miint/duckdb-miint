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

// Read a batch of query sequences from a table/view.
// Returns true if there are more rows to read, false if done.
// The batch is appended to the output, caller should clear if needed.
// offset is updated to the next position after reading.
bool ReadQueryBatch(ClientContext &context, const std::string &table_name, const SequenceTableSchema &schema,
                    idx_t batch_size, idx_t &offset, miint::SequenceRecordBatch &output);

// Read all read_ids for a shard from the read_to_shard table, ordered.
// Returns the IDs as a vector of strings for use with ReadBatchByIds.
std::vector<std::string> ReadShardIds(ClientContext &context, const std::string &read_to_shard_table,
                                      const std::string &shard_name);

// Read sequences for a known set of IDs by creating a temp table and joining.
// Reads ids[offset..offset+count] from the pre-materialized ID list,
// loads them into a temp table, and joins against query_table to fetch sequences.
void ReadBatchByIds(ClientContext &context, const std::string &query_table, const SequenceTableSchema &schema,
                    const std::vector<std::string> &ids, idx_t offset, idx_t count, miint::SequenceRecordBatch &output);

// Labels and single-end sequences loaded from a table for vsearch operations.
struct LoadedSingleEndSequences {
	std::vector<std::string> labels;
	std::vector<std::string> sequences;
};

// Load all (read_id, sequence1) pairs from a table via a separate connection.
// Used by vsearch-backed table functions (search, chimera, cluster) and
// align_mafft for single-end sequence loading.
// `schema` governs how the read_id column is extracted — VARCHAR rows pass
// through as strings; BIGINT rows are stringified into the carrier. Callers
// that haven't opted in to BIGINT should pass a schema captured via
// ValidateSequenceTableSchema (without allow_bigint), so id_type is VARCHAR.
// If strict=true: throws on NULL read_id, NULL sequence1, or empty sequence1.
// If strict=false: silently skips those rows.
// Always throws if the result set is empty after filtering.
// function_name is used in error messages (e.g. "cluster_sequences").
LoadedSingleEndSequences LoadSingleEndSequences(ClientContext &context, const std::string &table_name,
                                                const std::string &function_name, const SequenceTableSchema &schema,
                                                bool strict = false);

// Overload that runs against a caller-owned connection with an optional WHERE clause.
// Useful for per-sample callers: they already hold a per-thread Connection in LocalState
// and want to filter by a sample predicate. Pass `where_sql` without the leading "WHERE".
// An empty `where_sql` is equivalent to the context-based overload.
LoadedSingleEndSequences LoadSingleEndSequences(Connection &conn, const std::string &table_name,
                                                const std::string &function_name, const SequenceTableSchema &schema,
                                                bool strict, const std::string &where_sql);

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
