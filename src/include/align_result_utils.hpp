#pragma once
/*
 * Shared utility functions for alignment table functions (align_minimap2, align_bowtie2)
 * These set DuckDB result vectors from SAMRecordBatch data with offset/count support.
 */

#include "duckdb/common/types/vector.hpp"
#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// Set a VARCHAR vector from a string vector slice
inline void SetAlignResultString(Vector &result_vector, const std::vector<std::string> &values, idx_t offset,
                                 idx_t count) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = StringVector::AddString(result_vector, values[offset + j]);
	}
}

// Set a VARCHAR vector from a string vector slice, with empty strings as NULL
inline void SetAlignResultStringNullable(Vector &result_vector, const std::vector<std::string> &values, idx_t offset,
                                         idx_t count) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(count);

	for (idx_t j = 0; j < count; j++) {
		result_data[j] = StringVector::AddString(result_vector, values[offset + j]);
		if (!values[offset + j].empty()) {
			validity.SetValid(j);
		}
	}
}

// Set a UTINYINT vector from a uint8_t vector slice
inline void SetAlignResultUInt8(Vector &result_vector, const std::vector<uint8_t> &values, idx_t offset, idx_t count) {
	auto result_data = FlatVector::GetData<uint8_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

// Set a USMALLINT vector from a uint16_t vector slice
inline void SetAlignResultUInt16(Vector &result_vector, const std::vector<uint16_t> &values, idx_t offset,
                                 idx_t count) {
	auto result_data = FlatVector::GetData<uint16_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

// Set a BIGINT vector from an int64_t vector slice
inline void SetAlignResultInt64(Vector &result_vector, const std::vector<int64_t> &values, idx_t offset, idx_t count) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
	}
}

// Set a BIGINT vector from an int64_t vector slice, with -1 as NULL (for SAM tags)
inline void SetAlignResultInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values, idx_t offset,
                                        idx_t count) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(count);

	for (idx_t j = 0; j < count; j++) {
		result_data[j] = values[offset + j];
		// Tags return -1 when not present
		if (values[offset + j] != -1) {
			validity.SetValid(j);
		}
	}
}

} // namespace duckdb
