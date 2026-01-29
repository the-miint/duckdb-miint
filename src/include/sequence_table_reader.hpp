#pragma once

#include "Minimap2Aligner.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/main/client_context.hpp"
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
// If allow_paired is false, will throw if sequence2 column exists and has non-NULL data.
SequenceTableSchema ValidateSequenceTableSchema(ClientContext &context, const std::string &table_name,
                                                bool allow_paired = true);

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

// Read a batch of query sequences for a specific shard.
// Joins query_table with read_to_shard_table, filtering by shard_name.
// Returns true if there are more rows to read, false if done.
// offset is updated to the next position after reading.
bool ReadShardQueryBatch(ClientContext &context, const std::string &query_table, const std::string &read_to_shard_table,
                         const std::string &shard_name, const SequenceTableSchema &schema, idx_t batch_size,
                         idx_t &offset, miint::SequenceRecordBatch &output);

} // namespace duckdb
