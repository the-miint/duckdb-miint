#include "qc_functions.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/scalar_function.hpp"

#include "qc_algorithms.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>

namespace duckdb {

// ---------------------------------------------------------------------------
// Version stub (Cycle 0)
// ---------------------------------------------------------------------------
static constexpr const char *QC_VERSION_STRING = "qc 0.1.0 (port of fastp algorithms)";

static void QcVersionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto out = ConstantVector::GetData<string_t>(result);
	out[0] = StringVector::AddString(result, QC_VERSION_STRING);
}

// ---------------------------------------------------------------------------
// Shared types and helpers for the trim_quality_* scalars
// ---------------------------------------------------------------------------

static constexpr int32_t DEFAULT_QUAL_WINDOW = 4;
static constexpr int32_t DEFAULT_QUAL_MEAN = 20;

static LogicalType TrimResultStructType() {
	return LogicalType::STRUCT({{"sequence", LogicalType::VARCHAR},
	                            {"quality", LogicalType::LIST(LogicalType::UTINYINT)},
	                            {"trimmed_5p", LogicalType::UINTEGER},
	                            {"trimmed_3p", LogicalType::UINTEGER}});
}

// Locate the child slice for one row of a LIST(UTINYINT) input.
static void GetQualListSlice(Vector &list_vec, UnifiedVectorFormat &list_data, idx_t row_idx, const uint8_t *&out_data,
                             idx_t &out_length) {
	auto list_entries = UnifiedVectorFormat::GetData<list_entry_t>(list_data);
	auto mapped_idx = list_data.sel->get_index(row_idx);
	auto &entry = list_entries[mapped_idx];

	auto &child = ListVector::GetEntry(list_vec);
	auto child_data = FlatVector::GetData<uint8_t>(child);

	out_data = child_data + entry.offset;
	out_length = entry.length;
}

using TrimAlgorithm = miint::qc::TrimResult (*)(const std::uint8_t *, std::size_t, std::size_t, std::uint8_t);

// Execute one of the sliding-window trim algorithms over a DataChunk. Handles
// both the 2-arg form (defaults) and the 4-arg form (explicit window_size +
// mean_quality). Per-row scalar params are evaluated per-row, so they may be
// constants OR vary across rows.
static void TrimQualityExecuteImpl(DataChunk &args, Vector &result, TrimAlgorithm trim_fn, const char *fn_name) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 4;

	UnifiedVectorFormat seq_data, qual_data, window_data, mean_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	if (has_explicit_params) {
		args.data[2].ToUnifiedFormat(row_count, window_data);
		args.data[3].ToUnifiedFormat(row_count, mean_data);
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto window_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(window_data) : nullptr;
	auto mean_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(mean_data) : nullptr;

	auto &entries = StructVector::GetEntries(result);
	auto &seq_out_vec = *entries[0];
	auto &qual_out_vec = *entries[1];
	auto trimmed_5p_data = FlatVector::GetData<uint32_t>(*entries[2]);
	auto trimmed_3p_data = FlatVector::GetData<uint32_t>(*entries[3]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual_out_entries = FlatVector::GetData<list_entry_t>(qual_out_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_out_vec);

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool seq_valid = seq_data.validity.RowIsValid(si);
		bool qual_valid = qual_data.validity.RowIsValid(qi);
		bool window_valid = true;
		bool mean_valid = true;
		if (has_explicit_params) {
			window_valid = window_data.validity.RowIsValid(window_data.sel->get_index(i));
			mean_valid = mean_data.validity.RowIsValid(mean_data.sel->get_index(i));
		}

		if (!seq_valid || !qual_valid || !window_valid || !mean_valid) {
			result_validity.SetInvalid(i);
			qual_out_entries[i].offset = qual_child_offset;
			qual_out_entries[i].length = 0;
			continue;
		}

		const auto &seq = seq_ptr[si];
		const uint8_t *qptr;
		idx_t qlen;
		GetQualListSlice(args.data[1], qual_data, i, qptr, qlen);

		if (seq.GetSize() != qlen) {
			throw InvalidInputException("%s: sequence length (%llu) does not match quality length (%llu)", fn_name,
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		int32_t window = DEFAULT_QUAL_WINDOW;
		int32_t mean = DEFAULT_QUAL_MEAN;
		if (has_explicit_params) {
			window = window_ptr[window_data.sel->get_index(i)];
			mean = mean_ptr[mean_data.sel->get_index(i)];
		}
		if (window <= 0) {
			throw InvalidInputException("%s: window_size must be > 0 (got %d)", fn_name, window);
		}
		if (mean < 0 || mean > 93) {
			throw InvalidInputException("%s: mean_quality must be in 0..93 (got %d)", fn_name, mean);
		}

		miint::qc::TrimResult tr;
		try {
			tr = trim_fn(qptr, qlen, static_cast<std::size_t>(window), static_cast<std::uint8_t>(mean));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("%s: %s", fn_name, e.what());
		}

		// Trimmed sequence
		const idx_t kept_len = tr.end - tr.start;
		FlatVector::GetData<string_t>(seq_out_vec)[i] =
		    StringVector::AddString(seq_out_vec, seq.GetData() + tr.start, kept_len);

		// Trimmed quality (slice of input)
		ListVector::Reserve(qual_out_vec, qual_child_offset + kept_len);
		auto &qual_child = ListVector::GetEntry(qual_out_vec);
		auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
		for (idx_t k = 0; k < kept_len; k++) {
			qual_child_data[qual_child_offset + k] = qptr[tr.start + k];
		}
		qual_out_entries[i].offset = qual_child_offset;
		qual_out_entries[i].length = kept_len;
		qual_child_offset += kept_len;
		ListVector::SetListSize(qual_out_vec, qual_child_offset);

		trimmed_5p_data[i] = static_cast<uint32_t>(tr.start);
		trimmed_3p_data[i] = static_cast<uint32_t>(qlen - tr.end);
	}
}

static void TrimQuality3pFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimQualityExecuteImpl(args, result, &miint::qc::SlidingWindowTrimmer::trim_3p, "trim_quality_3p");
}

static void TrimQuality5pFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimQualityExecuteImpl(args, result, &miint::qc::SlidingWindowTrimmer::trim_5p, "trim_quality_5p");
}

static void TrimQualitySlidingFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimQualityExecuteImpl(args, result, &miint::qc::SlidingWindowTrimmer::trim_sliding, "trim_quality_sliding");
}

// ---------------------------------------------------------------------------
// Shared helpers for the trim_poly* scalars
// ---------------------------------------------------------------------------

static constexpr int32_t DEFAULT_POLY_MIN_LEN = 10;
static constexpr int32_t DEFAULT_POLY_MAX_MISMATCH = 5;
static constexpr int32_t DEFAULT_POLYG_MAX_WINDOW_MEAN_Q = 5;

// Write the trim_result struct fields and copy the kept slice of seq+qual.
// Reused by trim_polyg and trim_polyx scalar implementations.
static void WriteTrimRow(idx_t i, const string_t &seq, const uint8_t *qptr, idx_t qlen, miint::qc::TrimResult tr,
                         Vector &seq_out_vec, Vector &qual_out_vec, list_entry_t *qual_out_entries,
                         idx_t &qual_child_offset, uint32_t *trimmed_5p_data, uint32_t *trimmed_3p_data) {
	const idx_t kept_len = tr.end - tr.start;
	FlatVector::GetData<string_t>(seq_out_vec)[i] =
	    StringVector::AddString(seq_out_vec, seq.GetData() + tr.start, kept_len);

	ListVector::Reserve(qual_out_vec, qual_child_offset + kept_len);
	auto &qual_child = ListVector::GetEntry(qual_out_vec);
	auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
	for (idx_t k = 0; k < kept_len; k++) {
		qual_child_data[qual_child_offset + k] = qptr[tr.start + k];
	}
	qual_out_entries[i].offset = qual_child_offset;
	qual_out_entries[i].length = kept_len;
	qual_child_offset += kept_len;
	ListVector::SetListSize(qual_out_vec, qual_child_offset);

	trimmed_5p_data[i] = static_cast<uint32_t>(tr.start);
	trimmed_3p_data[i] = static_cast<uint32_t>(qlen - tr.end);
}

// ---------------------------------------------------------------------------
// trim_polyg
// ---------------------------------------------------------------------------
static void TrimPolygExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 5; // seq, qual, min_len, max_mm, max_window_q

	UnifiedVectorFormat seq_data, qual_data, p1_data, p2_data, p3_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	if (has_explicit_params) {
		args.data[2].ToUnifiedFormat(row_count, p1_data);
		args.data[3].ToUnifiedFormat(row_count, p2_data);
		args.data[4].ToUnifiedFormat(row_count, p3_data);
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto p1_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p1_data) : nullptr;
	auto p2_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p2_data) : nullptr;
	auto p3_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p3_data) : nullptr;

	auto &entries = StructVector::GetEntries(result);
	auto &seq_out_vec = *entries[0];
	auto &qual_out_vec = *entries[1];
	auto trimmed_5p_data = FlatVector::GetData<uint32_t>(*entries[2]);
	auto trimmed_3p_data = FlatVector::GetData<uint32_t>(*entries[3]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual_out_entries = FlatVector::GetData<list_entry_t>(qual_out_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_out_vec);

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi);
		if (has_explicit_params) {
			all_valid = all_valid && p1_data.validity.RowIsValid(p1_data.sel->get_index(i)) &&
			            p2_data.validity.RowIsValid(p2_data.sel->get_index(i)) &&
			            p3_data.validity.RowIsValid(p3_data.sel->get_index(i));
		}
		if (!all_valid) {
			result_validity.SetInvalid(i);
			qual_out_entries[i].offset = qual_child_offset;
			qual_out_entries[i].length = 0;
			continue;
		}

		const auto &seq = seq_ptr[si];
		const uint8_t *qptr;
		idx_t qlen;
		GetQualListSlice(args.data[1], qual_data, i, qptr, qlen);
		if (seq.GetSize() != qlen) {
			throw InvalidInputException("trim_polyg: sequence length (%llu) does not match quality length (%llu)",
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		int32_t min_len = DEFAULT_POLY_MIN_LEN;
		int32_t max_mismatch = DEFAULT_POLY_MAX_MISMATCH;
		int32_t max_window_q = DEFAULT_POLYG_MAX_WINDOW_MEAN_Q;
		if (has_explicit_params) {
			min_len = p1_ptr[p1_data.sel->get_index(i)];
			max_mismatch = p2_ptr[p2_data.sel->get_index(i)];
			max_window_q = p3_ptr[p3_data.sel->get_index(i)];
		}
		if (min_len < 1) {
			throw InvalidInputException("trim_polyg: min_len must be >= 1 (got %d)", min_len);
		}
		if (max_mismatch < 0) {
			throw InvalidInputException("trim_polyg: max_mismatch must be >= 0 (got %d)", max_mismatch);
		}
		if (max_window_q < 0 || max_window_q > 255) {
			throw InvalidInputException("trim_polyg: max_window_mean_q must be 0..255 (got %d)", max_window_q);
		}

		auto tr = miint::qc::PolyXScanner::scan_polyg(
		    reinterpret_cast<const std::uint8_t *>(seq.GetData()), qptr, qlen, static_cast<std::size_t>(min_len),
		    static_cast<std::uint32_t>(max_mismatch), static_cast<std::uint8_t>(max_window_q));

		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
	}
}

// ---------------------------------------------------------------------------
// trim_polyx
// ---------------------------------------------------------------------------
static void TrimPolyxExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 4; // seq, qual, min_len, max_mm

	UnifiedVectorFormat seq_data, qual_data, p1_data, p2_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	if (has_explicit_params) {
		args.data[2].ToUnifiedFormat(row_count, p1_data);
		args.data[3].ToUnifiedFormat(row_count, p2_data);
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto p1_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p1_data) : nullptr;
	auto p2_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p2_data) : nullptr;

	auto &entries = StructVector::GetEntries(result);
	auto &seq_out_vec = *entries[0];
	auto &qual_out_vec = *entries[1];
	auto trimmed_5p_data = FlatVector::GetData<uint32_t>(*entries[2]);
	auto trimmed_3p_data = FlatVector::GetData<uint32_t>(*entries[3]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual_out_entries = FlatVector::GetData<list_entry_t>(qual_out_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_out_vec);

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi);
		if (has_explicit_params) {
			all_valid = all_valid && p1_data.validity.RowIsValid(p1_data.sel->get_index(i)) &&
			            p2_data.validity.RowIsValid(p2_data.sel->get_index(i));
		}
		if (!all_valid) {
			result_validity.SetInvalid(i);
			qual_out_entries[i].offset = qual_child_offset;
			qual_out_entries[i].length = 0;
			continue;
		}

		const auto &seq = seq_ptr[si];
		const uint8_t *qptr;
		idx_t qlen;
		GetQualListSlice(args.data[1], qual_data, i, qptr, qlen);
		if (seq.GetSize() != qlen) {
			throw InvalidInputException("trim_polyx: sequence length (%llu) does not match quality length (%llu)",
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		int32_t min_len = DEFAULT_POLY_MIN_LEN;
		int32_t max_mismatch = DEFAULT_POLY_MAX_MISMATCH;
		if (has_explicit_params) {
			min_len = p1_ptr[p1_data.sel->get_index(i)];
			max_mismatch = p2_ptr[p2_data.sel->get_index(i)];
		}
		if (min_len < 1) {
			throw InvalidInputException("trim_polyx: min_len must be >= 1 (got %d)", min_len);
		}
		if (max_mismatch < 0) {
			throw InvalidInputException("trim_polyx: max_mismatch must be >= 0 (got %d)", max_mismatch);
		}

		auto tr = miint::qc::PolyXScanner::scan_polyx(reinterpret_cast<const std::uint8_t *>(seq.GetData()), qlen,
		                                              static_cast<std::size_t>(min_len),
		                                              static_cast<std::uint32_t>(max_mismatch));

		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
	}
}

// Register a trim_quality_* function with both the 2-arg (defaults) and 4-arg
// (explicit window_size + mean_quality) overloads.
static void RegisterTrimQualityFamily(ExtensionLoader &loader, const std::string &name, scalar_function_t fn) {
	ScalarFunctionSet set(name);

	ScalarFunction two_arg(name, {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT)},
	                       TrimResultStructType(), fn);
	two_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	set.AddFunction(two_arg);

	ScalarFunction four_arg(
	    name,
	    {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::INTEGER, LogicalType::INTEGER},
	    TrimResultStructType(), fn);
	four_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	set.AddFunction(four_arg);

	loader.RegisterFunction(set);
}

// ---------------------------------------------------------------------------
// Registration entry point
// ---------------------------------------------------------------------------
void QcFunctions::Register(ExtensionLoader &loader) {
	ScalarFunction qc_version("qc_version", {}, LogicalType::VARCHAR, QcVersionFunction);
	loader.RegisterFunction(qc_version);

	RegisterTrimQualityFamily(loader, "trim_quality_3p", TrimQuality3pFunction);
	RegisterTrimQualityFamily(loader, "trim_quality_5p", TrimQuality5pFunction);
	RegisterTrimQualityFamily(loader, "trim_quality_sliding", TrimQualitySlidingFunction);

	// trim_polyg: 2-arg (defaults) + 5-arg (min_len, max_mismatch, max_window_mean_q).
	{
		ScalarFunctionSet set("trim_polyg");
		ScalarFunction two_arg("trim_polyg", {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT)},
		                       TrimResultStructType(), TrimPolygExecute);
		two_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(two_arg);

		ScalarFunction five_arg("trim_polyg",
		                        {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::INTEGER,
		                         LogicalType::INTEGER, LogicalType::INTEGER},
		                        TrimResultStructType(), TrimPolygExecute);
		five_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(five_arg);
		loader.RegisterFunction(set);
	}

	// trim_polyx: 2-arg (defaults) + 4-arg (min_len, max_mismatch).
	{
		ScalarFunctionSet set("trim_polyx");
		ScalarFunction two_arg("trim_polyx", {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT)},
		                       TrimResultStructType(), TrimPolyxExecute);
		two_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(two_arg);

		ScalarFunction four_arg("trim_polyx",
		                        {LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::INTEGER,
		                         LogicalType::INTEGER},
		                        TrimResultStructType(), TrimPolyxExecute);
		four_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(four_arg);
		loader.RegisterFunction(set);
	}
}

} // namespace duckdb
