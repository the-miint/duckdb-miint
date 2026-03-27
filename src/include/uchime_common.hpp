#pragma once

#include "ChimeraDetector.hpp"

#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

#include <string>
#include <vector>

namespace duckdb {

// Output column names and types for the 18-column uchimeout format.
// Shared between uchime_ref and uchime_denovo.
std::vector<std::string> GetUchimeOutputNames();
std::vector<LogicalType> GetUchimeOutputTypes();

// Write a batch of UchimeResults into a DataChunk.
// Returns the number of rows written (min of count and remaining results).
idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset,
                          idx_t count);

} // namespace duckdb
