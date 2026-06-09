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
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector_operations/vector_operations.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// True for the set of column types accepted as an identifier column.
inline bool IsAllowedIdType(const LogicalType &t) {
	return t.id() == LogicalTypeId::VARCHAR || t.id() == LogicalTypeId::BIGINT || t.id() == LogicalTypeId::UUID;
}

// Human-readable description of the accepted identifier-column types, kept in
// lockstep with IsAllowedIdType so error messages never drift from the actual
// accepted set. Used as the `%s` / type-desc in every "must be ..." id-type
// rejection across the codebase (single source of truth).
inline const char *AllowedIdTypeList() {
	return "VARCHAR, BIGINT, or UUID";
}

// Parse a canonical UUID string back to its INT128 storage value. duckdb's
// UUID::FromString returns false (it does NOT throw) on malformed input, and the
// argument-less convenience overload silently swallows that — so we check the
// bool here and fail loud, mirroring the BIGINT codec's throw-on-unparseable
// contract. The SAM sentinels ("", "*", "=") are not valid UUIDs and must be
// handled by the caller before reaching here; passing one (correctly) throws.
inline hugeint_t ParseUuidOrThrow(const std::string &s) {
	hugeint_t result;
	if (!UUID::FromString(s, result)) {
		throw InvalidInputException("cannot parse value '%s' as UUID", s);
	}
	return result;
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
	} else if (id_type.id() == LogicalTypeId::UUID) {
		// UUID is stored physically as INT128; stringify to the canonical 36-char
		// form so the carrier holds the same lexical id the user sees.
		auto data = UnifiedVectorFormat::GetData<hugeint_t>(fmt);
		for (idx_t i = 0; i < n; i++) {
			auto idx = fmt.sel->get_index(i);
			if (!fmt.validity.RowIsValid(idx)) {
				out_nulls[i] = true;
				continue;
			}
			out_values[i] = UUID::ToString(data[idx]);
		}
	} else {
		throw InternalException("ExtractIdColumnAsStrings: unsupported id type '%s' (must be %s)", id_type.ToString(),
		                        AllowedIdTypeList());
	}
}

// Emit one id cell at `row` into flat vector `out`, dispatching on `id_type`.
// This is the single per-row egress primitive shared by EmitIdColumnFromStrings
// (minimap2/sortmerna output) and the bowtie2 Arrow emit path — keeping the
// VARCHAR/BIGINT dispatch in one place so adding an id type touches one function
// rather than several divergent codecs.
// VARCHAR: writes the bytes verbatim, row marked valid (empty string passes
//          through since VARCHAR ingress carries no NULL encoding here).
// BIGINT:  parses via ParseIdAsInt64; empty string and "*" become NULL; any
//          other unparseable value throws InvalidInputException.
// The caller must resolve SAM sentinels (e.g. "=" for mate_reference) before
// invoking — the codec deliberately fails loud on "=" to surface that contract.
// Validity is always set explicitly (valid or invalid) so reusing the vector is
// safe. Precondition: `out` is a FlatVector sized for at least `row + 1` rows.
inline void EmitIdCell(Vector &out, idx_t row, const std::string &s, const LogicalType &id_type) {
	auto &validity = FlatVector::Validity(out);
	if (id_type.id() == LogicalTypeId::VARCHAR) {
		FlatVector::GetData<string_t>(out)[row] = StringVector::AddString(out, s);
		validity.SetValid(row);
		return;
	}
	if (id_type.id() == LogicalTypeId::BIGINT) {
		auto out_data = FlatVector::GetData<int64_t>(out);
		try {
			auto parsed = ::miint::ParseIdAsInt64(s);
			if (parsed.has_value()) {
				out_data[row] = *parsed;
				validity.SetValid(row);
			} else {
				out_data[row] = 0;
				validity.SetInvalid(row);
			}
		} catch (const std::exception &e) {
			throw InvalidInputException("EmitIdCell: cannot parse value '%s' as BIGINT: %s", s, e.what());
		}
		return;
	}
	if (id_type.id() == LogicalTypeId::UUID) {
		auto out_data = FlatVector::GetData<hugeint_t>(out);
		if (s.empty() || s == "*") {
			out_data[row] = 0;
			validity.SetInvalid(row);
		} else {
			out_data[row] = ParseUuidOrThrow(s);
			validity.SetValid(row);
		}
		return;
	}
	throw InternalException("EmitIdCell: unsupported id type '%s' (must be %s)", id_type.ToString(),
	                        AllowedIdTypeList());
}

// Emit `ids[offset..offset+count)` into `out` via EmitIdCell. See EmitIdCell for
// the per-type contract and the SAM-sentinel caveat.
// Precondition: `out` is a FlatVector sized for at least `count` rows.
inline void EmitIdColumnFromStrings(Vector &out, const std::vector<std::string> &ids, idx_t offset, idx_t count,
                                    const LogicalType &id_type) {
	for (idx_t j = 0; j < count; j++) {
		EmitIdCell(out, j, ids[offset + j], id_type);
	}
}

} // namespace duckdb
