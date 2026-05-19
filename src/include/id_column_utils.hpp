#pragma once
//
// DuckDB-aware wrappers around the id-column codec.
//
// duckdb-miint accepts identifier columns of type VARCHAR or BIGINT and keeps
// them as `std::vector<std::string>` in its in-memory carriers, regardless of
// the SQL-level type. These helpers handle the translation at the ingress
// (DataChunk -> vector<string>) and egress (vector<string> -> Vector)
// boundaries.
//
// Pure parse/format logic lives in id_column_codec.hpp so the unit-test
// target can link it without the duckdb library.

#include "id_column_codec.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// True for the set of column types accepted as an identifier column.
inline bool IsAllowedIdType(const LogicalType &t) {
	return t.id() == LogicalTypeId::VARCHAR || t.id() == LogicalTypeId::BIGINT;
}

// Extract `chunk.data[col_idx]` as strings, dispatching on `id_type`.
// Always resizes `out_values` and `out_nulls` to `chunk.size()`.
// NULL rows get an empty string in `out_values` and `true` in `out_nulls`.
// VARCHAR ingress copies the bytes; BIGINT ingress stringifies via the codec.
// Caller is responsible for ensuring `id_type.id()` matches the column's
// storage type — typically by capturing the type at bind time via
// ValidateSequenceTableSchema and threading it through.
inline void ExtractIdColumnAsStrings(DataChunk &chunk, idx_t col_idx, const LogicalType &id_type,
                                     std::vector<std::string> &out_values, std::vector<bool> &out_nulls) {
	const idx_t n = chunk.size();
	out_values.clear();
	out_values.resize(n);
	out_nulls.clear();
	out_nulls.resize(n, false);

	UnifiedVectorFormat fmt;
	chunk.data[col_idx].ToUnifiedFormat(n, fmt);

	if (id_type.id() == LogicalTypeId::VARCHAR) {
		auto data = UnifiedVectorFormat::GetData<string_t>(fmt);
		for (idx_t i = 0; i < n; i++) {
			auto idx = fmt.sel->get_index(i);
			if (!fmt.validity.RowIsValid(idx)) {
				out_nulls[i] = true;
				continue;
			}
			out_values[i] = data[idx].GetString();
		}
	} else if (id_type.id() == LogicalTypeId::BIGINT) {
		auto data = UnifiedVectorFormat::GetData<int64_t>(fmt);
		for (idx_t i = 0; i < n; i++) {
			auto idx = fmt.sel->get_index(i);
			if (!fmt.validity.RowIsValid(idx)) {
				out_nulls[i] = true;
				continue;
			}
			// Absolute `::miint::` rather than `miint::` — this header is included
			// from translation units that also reference duckdb::miint (e.g.
			// gpl_boundary callers like align_bowtie2_daemon_common.cpp). Without
			// the absolute qualifier, the unqualified lookup resolves to
			// duckdb::miint first and fails to find the codec functions there.
			out_values[i] = ::miint::FormatIdFromInt64(data[idx]);
		}
	} else {
		throw InternalException("ExtractIdColumnAsStrings: unsupported id type '%s' (must be VARCHAR or BIGINT)",
		                        id_type.ToString());
	}
}

// Emit `ids[offset..offset+count)` into `out`, dispatching on `id_type`.
// VARCHAR: emits bytes verbatim (empty string passes through; validity is
//          reset to all-valid since VARCHAR ingress doesn't carry a NULL
//          encoding here).
// BIGINT:  parses each string via ParseIdAsInt64; empty string and "*" become
//          NULL; any other unparseable value throws InvalidInputException.
// The caller must resolve SAM sentinels (e.g. "=" for mate_reference) before
// invoking — the codec deliberately fails loud on "=" to surface that contract.
// Precondition: `out` is a FlatVector sized for at least `count` rows. Both
// branches set the validity mask explicitly so reusing the vector is safe.
inline void EmitIdColumnFromStrings(Vector &out, const std::vector<std::string> &ids, idx_t offset, idx_t count,
                                    const LogicalType &id_type) {
	if (id_type.id() == LogicalTypeId::VARCHAR) {
		auto out_data = FlatVector::GetData<string_t>(out);
		FlatVector::Validity(out).SetAllValid(count);
		for (idx_t j = 0; j < count; j++) {
			out_data[j] = StringVector::AddString(out, ids[offset + j]);
		}
		return;
	}
	if (id_type.id() == LogicalTypeId::BIGINT) {
		auto out_data = FlatVector::GetData<int64_t>(out);
		auto &validity = FlatVector::Validity(out);
		validity.SetAllInvalid(count);
		for (idx_t j = 0; j < count; j++) {
			const auto &s = ids[offset + j];
			try {
				auto parsed = ::miint::ParseIdAsInt64(s);
				if (parsed.has_value()) {
					out_data[j] = *parsed;
					validity.SetValid(j);
				} else {
					out_data[j] = 0;
				}
			} catch (const std::exception &e) {
				throw InvalidInputException(
				    "EmitIdColumnFromStrings: cannot parse value '%s' at row %llu as BIGINT: %s", s,
				    static_cast<unsigned long long>(offset + j), e.what());
			}
		}
		return;
	}
	throw InternalException("EmitIdColumnFromStrings: unsupported id type '%s' (must be VARCHAR or BIGINT)",
	                        id_type.ToString());
}

} // namespace duckdb
