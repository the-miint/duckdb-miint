#pragma once

// Shared helpers between align_bowtie2 and align_bowtie2_sharded once both
// route through the gpl-boundary daemon's bowtie2-align tool. Lives in its
// own translation unit so neither caller has to drag in the other's
// table-function plumbing.

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/table_function.hpp"

#include <cstdint>
#include <string>
#include <vector>

// nanoarrow C ABI is opaque here; the source file pulls in the full header.
struct ArrowSchema;
struct ArrowArray;

namespace duckdb {
namespace bt2_daemon {

// gpl-boundary bowtie2-align output schema_version we wire against
// (v0.2.0 commit a28cf56). Strict equality check at session bind.
constexpr uint32_t kAlignSchemaVersion = 2;

constexpr int kNumOutputColumns = 21;

// User-facing column order; MUST match the daemon's bowtie2-align v2 wire
// schema so wire column i → output column i.
extern const char *const kOutputColumnNames[kNumOutputColumns];

// Arrow C Data Interface format string per column; matches the daemon's
// bowtie2-align v2 wire schema. Validated structurally so a type drift on the
// daemon side fails loud rather than reinterpreting bytes through
// `read_fixed<T>` with the wrong width.
extern const char *const kOutputColumnFormats[kNumOutputColumns];

// Populate the 21-column output schema (matching the daemon's wire schema).
void PopulateOutputSchema(std::vector<std::string> &names, std::vector<LogicalType> &types);

// Validate that an Arrow IPC schema returned by the daemon matches the
// expected 21 named columns. Throws IOException on drift.
void ValidateOutputSchema(const ArrowSchema &schema);

// Describe the column layout for a query Arrow batch. The query Arrow
// batch always has read_id+sequence1; sequence2/qual1/qual2 are optional
// based on the user's input table.
struct QueryArrowSchema {
	bool has_sequence2 = false;
	bool has_qual1 = false;
	bool has_qual2 = false;
	int num_columns() const {
		return 2 + (has_sequence2 ? 1 : 0) + (has_qual1 ? 1 : 0) + (has_qual2 ? 1 : 0);
	}
};

// Row-oriented query batch. Optional columns use parallel `_valid` vectors
// (1 = present, 0 = NULL). When the column flag is false in QueryArrowSchema,
// the corresponding vectors are empty / unused.
struct QueryBatch {
	std::vector<std::string> read_ids;
	std::vector<std::string> sequence1;
	std::vector<std::string> sequence2;
	std::vector<int8_t> sequence2_valid;
	std::vector<std::string> qual1;
	std::vector<int8_t> qual1_valid;
	std::vector<std::string> qual2;
	std::vector<int8_t> qual2_valid;
};

// Encode a QueryBatch as an Arrow IPC stream consumable by the daemon's
// bowtie2-align input schema. Throws InternalException on encoder failure.
std::vector<uint8_t> BuildQueryIpc(const QueryBatch &qb, const QueryArrowSchema &schema_flags);

// Drain `to_emit` rows starting at `row_start` from a decoded daemon Arrow
// batch into a DuckDB DataChunk. Assumes `output` has 21 columns matching
// PopulateOutputSchema's types. Caller is responsible for SetCardinality
// before calling (we just write into the vectors).
//
// Tag widening (Int32 → BIGINT) and nullable-Utf8 decoding match the
// daemon's wire schema.
void EmitChunkRows(DataChunk &output, idx_t to_emit, idx_t row_start, const ArrowArray &batch);

} // namespace bt2_daemon
} // namespace duckdb
