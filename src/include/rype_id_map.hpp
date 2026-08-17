#pragma once
//
// Surrogate-id -> caller-id map for the RYpe table functions.
//
// RYpe treats the `id` column of its Arrow input as an opaque i64 and echoes it
// back as `query_id`: it indexes every internal structure by the read's position
// within the batch and only reads `query_ids[read_idx]` to stamp the output
// (ext/rype/src/classify/sharded.rs, classify/merge_join.rs). miint therefore
// feeds it a dense surrogate and keeps the caller's real identifier here, indexed
// by that surrogate.
//
// Stored in the id column's own storage type rather than as
// std::vector<std::string>. A string carrier costs 32+ bytes per read plus a
// Value::ToString() on ingress and a re-parse on egress; on a 2M-read corpus that
// is ~67 MB of C++ heap DuckDB's memory_limit does not account for. BIGINT costs
// 8 bytes here and round-trips with no parsing at all.
//
// NULL handling reproduces the pre-existing carrier semantics exactly:
//   BIGINT / UUID : a NULL id round-trips to SQL NULL.
//   VARCHAR       : a NULL id round-trips to the empty string, not NULL.
// The asymmetry predates this class — it fell out of the empty-string carrier the
// old codec used — and test/sql/rype_classify_bigint.test and
// rype_extract_bigint.test pin the BIGINT/UUID half ("NULL ids emit SQL NULL
// instead of aborting the query"). Preserved deliberately rather than fixed in
// passing.

#include "id_column_utils.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/common/types/vector.hpp"

#include <string>
#include <vector>

namespace duckdb {

class RypeIdMap {
public:
	explicit RypeIdMap(LogicalType id_type_p) : id_type(std::move(id_type_p)) {
		if (!IsAllowedIdType(id_type)) {
			throw InternalException("RypeIdMap: unsupported id type '%s' (must be %s)", id_type.ToString(),
			                        AllowedIdTypeList());
		}
		if (id_type.id() == LogicalTypeId::VARCHAR) {
			// Offsets are the usual n+1 form, so entry i spans [offsets[i], offsets[i+1]).
			varchar_offsets.push_back(0);
		}
	}

	//! Append rows [from, to) of `fmt`, a UnifiedVectorFormat over the caller's id
	//! column. Appends in logical row order, so the surrogate for row i is
	//! size() + (i - from) as observed before the call.
	void AppendRange(const UnifiedVectorFormat &fmt, idx_t from, idx_t to) {
		switch (id_type.id()) {
		case LogicalTypeId::BIGINT: {
			auto data = UnifiedVectorFormat::GetData<int64_t>(fmt);
			for (idx_t i = from; i < to; i++) {
				auto idx = fmt.sel->get_index(i);
				const bool is_null = !fmt.validity.RowIsValid(idx);
				nulls.push_back(is_null);
				i64_values.push_back(is_null ? 0 : data[idx]);
			}
			break;
		}
		case LogicalTypeId::UUID: {
			auto data = UnifiedVectorFormat::GetData<hugeint_t>(fmt);
			for (idx_t i = from; i < to; i++) {
				auto idx = fmt.sel->get_index(i);
				const bool is_null = !fmt.validity.RowIsValid(idx);
				nulls.push_back(is_null);
				i128_values.push_back(is_null ? hugeint_t(0) : data[idx]);
			}
			break;
		}
		default: { // VARCHAR
			auto data = UnifiedVectorFormat::GetData<string_t>(fmt);
			for (idx_t i = from; i < to; i++) {
				auto idx = fmt.sel->get_index(i);
				const bool is_null = !fmt.validity.RowIsValid(idx);
				nulls.push_back(is_null);
				if (!is_null) {
					varchar_arena.append(data[idx].GetData(), data[idx].GetSize());
				}
				varchar_offsets.push_back(varchar_arena.size());
			}
			break;
		}
		}
		count += to - from;
	}

	idx_t size() const {
		return count;
	}

	//! Write entry `idx` into row `row` of the flat output vector `out`, whose
	//! type must be the id type this map was constructed with.
	void Emit(Vector &out, idx_t row, idx_t idx) const {
		D_ASSERT(idx < count);
		auto &validity = FlatVector::Validity(out);
		switch (id_type.id()) {
		case LogicalTypeId::BIGINT:
			FlatVector::GetData<int64_t>(out)[row] = i64_values[idx];
			validity.Set(row, !nulls[idx]);
			return;
		case LogicalTypeId::UUID:
			FlatVector::GetData<hugeint_t>(out)[row] = i128_values[idx];
			validity.Set(row, !nulls[idx]);
			return;
		default: { // VARCHAR: a NULL id stored zero bytes, so it emits '' and stays valid
			const auto begin = varchar_offsets[idx];
			const auto length = varchar_offsets[idx + 1] - begin;
			FlatVector::GetData<string_t>(out)[row] =
			    StringVector::AddString(out, varchar_arena.data() + begin, length);
			validity.SetValid(row);
			return;
		}
		}
	}

	//! Human-readable form of entry `idx`. For error messages only — the egress
	//! path uses Emit, which never round-trips through a string.
	std::string ToString(idx_t idx) const {
		D_ASSERT(idx < count);
		if (nulls[idx]) {
			return "NULL";
		}
		switch (id_type.id()) {
		case LogicalTypeId::BIGINT:
			return std::to_string(i64_values[idx]);
		case LogicalTypeId::UUID:
			return UUID::ToString(i128_values[idx]);
		default:
			return varchar_arena.substr(varchar_offsets[idx], varchar_offsets[idx + 1] - varchar_offsets[idx]);
		}
	}

private:
	LogicalType id_type;
	idx_t count = 0;

	std::vector<int64_t> i64_values;    // BIGINT
	std::vector<hugeint_t> i128_values; // UUID
	std::string varchar_arena;          // VARCHAR bytes, back to back
	std::vector<idx_t> varchar_offsets; // VARCHAR, size() + 1 entries
	std::vector<bool> nulls;            // bit-packed
};

} // namespace duckdb
