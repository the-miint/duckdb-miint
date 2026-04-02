#pragma once

#include "Minimap2Aligner.hpp"
#include "SequenceRecord.hpp"
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
};

// Validate that a table/view has required columns for sequence data.
// Returns schema information about what optional columns are present.
// Throws BinderException if required columns are missing or have wrong types.
SequenceTableSchema ValidateSequenceTableSchema(ClientContext &context, const std::string &table_name);

// Read all subjects from a table/view into memory.
// Subjects cannot be paired-end (sequence2 must be NULL for all rows).
// Throws InvalidInputException if sequence2 contains non-NULL values.
std::vector<miint::AlignmentSubject> ReadSubjectTable(ClientContext &context, const std::string &table_name);

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
// Used by vsearch-backed table functions (search, chimera, cluster) that operate
// on single-end sequences only.
// If strict=true: throws on NULL read_id, NULL sequence1, or empty sequence1.
// If strict=false: silently skips those rows.
// Always throws if the result set is empty after filtering.
// function_name is used in error messages (e.g. "cluster_sequences").
LoadedSingleEndSequences LoadSingleEndSequences(ClientContext &context, const std::string &table_name,
                                                const std::string &function_name, bool strict = false);

// Streaming query sequence reader for lazy sub-batching.
// Produces sub-batches on demand from a streaming query result.
// Thread-safe: multiple threads can call FetchSubBatch() concurrently.
class QuerySequenceStream {
public:
	QuerySequenceStream(ClientContext &context, const std::string &table_name, const SequenceTableSchema &schema,
	                    idx_t sub_batch_size = STANDARD_VECTOR_SIZE);

	// Fetch the next sub-batch. Returns an empty batch when the stream is exhausted.
	// Thread-safe — serializes access to the underlying stream via mutex.
	miint::SequenceRecordBatch FetchSubBatch();

private:
	Connection conn_;
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
};

} // namespace duckdb
