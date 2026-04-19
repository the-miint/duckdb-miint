#pragma once

#include "SortMeRNAAligner.hpp"
#include "align_result_utils.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include <string>
#include <vector>

namespace duckdb {

// Output schema for align_sortmerna_rrna — carries fields SAM cannot: e_value,
// identity, coverage, edit_distance, segment_idx. Paired-end produces two rows
// per input row with segment_idx 0 (fwd) and 1 (rev).
inline std::vector<std::string> GetSortMeRNARRNAOutputNames() {
	return {"read_id", "aligned", "strand",   "ref_name", "ref_start",     "ref_end",    "cigar",
	        "score",   "e_value", "identity", "coverage", "edit_distance", "segment_idx"};
}

inline std::vector<LogicalType> GetSortMeRNARRNAOutputTypes() {
	return {LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR, LogicalType::INTEGER,
	        LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER,
	        LogicalType::INTEGER};
}

// Emit [offset, offset+count) rows of `batch` into the output DataChunk.
inline void OutputSortMeRNARRNABatch(DataChunk &output, const miint::SortMeRNAResultBatch &batch, idx_t offset,
                                     idx_t count) {
	idx_t col = 0;
	SetAlignResultString(output.data[col++], batch.read_ids, offset, count);
	SetAlignResultInt32(output.data[col++], batch.aligned, offset, count);
	SetAlignResultInt32(output.data[col++], batch.strands, offset, count);
	SetAlignResultString(output.data[col++], batch.ref_names, offset, count);
	SetAlignResultInt32(output.data[col++], batch.ref_starts, offset, count);
	SetAlignResultInt32(output.data[col++], batch.ref_ends, offset, count);
	SetAlignResultString(output.data[col++], batch.cigars, offset, count);
	SetAlignResultInt32(output.data[col++], batch.scores, offset, count);
	SetAlignResultDouble(output.data[col++], batch.e_values, offset, count);
	SetAlignResultDouble(output.data[col++], batch.identities, offset, count);
	SetAlignResultDouble(output.data[col++], batch.coverages, offset, count);
	SetAlignResultInt32(output.data[col++], batch.edit_distances, offset, count);
	SetAlignResultInt32(output.data[col++], batch.segment_indices, offset, count);
	output.SetCardinality(count);
}

} // namespace duckdb
