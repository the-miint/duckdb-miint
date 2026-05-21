#include "qc_functions.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/scalar_function.hpp"

#include "qc_algorithms.hpp"
#include "sequence_utils.hpp"

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

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
// Shared types and helpers for all trim_* scalars
// ---------------------------------------------------------------------------

static constexpr int32_t DEFAULT_QUAL_WINDOW = 4;
static constexpr int32_t DEFAULT_QUAL_MEAN = 20;
// fastp's documented window-size range is 1..1000. The upper bound also keeps
// threshold_sum (mean_quality * window_size) safely within uint32_t.
static constexpr int32_t MAX_QUAL_WINDOW = 1000;

static constexpr int32_t DEFAULT_POLY_MIN_LEN = 10;
static constexpr int32_t DEFAULT_POLY_MAX_MISMATCH = 5;
// Default quality-aware gate for polyG: trim region's mean Phred must be <=
// this threshold. Pass max_window_mean_q=QUAL_GATE_DISABLED (the max valid
// Phred) to make the gate a no-op since any real Phred score satisfies <= 93.
static constexpr int32_t DEFAULT_POLYG_MAX_WINDOW_MEAN_Q = 5;
static constexpr int32_t QUAL_GATE_DISABLED = 93;

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

// Write one row of the trim_result struct: trimmed sequence, trimmed quality
// (as a sliced LIST(UTINYINT)), and the trimmed_5p/trimmed_3p counts. Shared
// by every trim_* scalar.
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
		if (window <= 0 || window > MAX_QUAL_WINDOW) {
			throw InvalidInputException("%s: window_size must be in 1..%d (got %d)", fn_name, MAX_QUAL_WINDOW, window);
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

		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
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
		if (max_window_q < 0 || max_window_q > QUAL_GATE_DISABLED) {
			throw InvalidInputException("trim_polyg: max_window_mean_q must be 0..%d (got %d); pass %d to disable",
			                            QUAL_GATE_DISABLED, max_window_q, QUAL_GATE_DISABLED);
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

		auto tr = miint::qc::PolyXScanner::scan_polyx(reinterpret_cast<const std::uint8_t *>(seq.GetData()),
		                                              seq.GetSize(), static_cast<std::size_t>(min_len),
		                                              static_cast<std::uint32_t>(max_mismatch));

		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
	}
}

// ---------------------------------------------------------------------------
// trim_adapters
// ---------------------------------------------------------------------------

// fastp's auto-scaling: shorter adapter lists get a more lenient min_match.
// Exposed so the caller can override by passing min_match >= 1.
static std::size_t default_min_match(std::size_t n_adapters) {
	if (n_adapters >= 4) {
		return 6;
	}
	if (n_adapters >= 2) {
		return 5;
	}
	return 4;
}

// Reverse-complement a DNA adapter, translating the underlying utility's
// std::runtime_error to InvalidInputException so it surfaces cleanly through
// SQL.
static std::string dna_revcomp_adapter(const std::string &s) {
	try {
		return miint::dna_reverse_complement(s);
	} catch (const std::runtime_error &e) {
		throw InvalidInputException("trim_adapters: invalid DNA base in adapter sequence (%s)", e.what());
	}
}

// Extract adapter strings from one row of the adapter argument. Handles both
// the VARCHAR overload (single adapter, never NULL — caller already guarded)
// and the LIST(VARCHAR) overload (multiple adapters; NULL or empty elements
// are silently skipped).
//
// Returns adapters as owned strings since we may also need to extend with
// reverse complements.
static std::vector<std::string> ExtractAdapters(Vector &arg_vec, UnifiedVectorFormat &arg_data, idx_t row_idx,
                                                bool is_list) {
	std::vector<std::string> out;
	if (is_list) {
		auto list_entries = UnifiedVectorFormat::GetData<list_entry_t>(arg_data);
		auto mapped = arg_data.sel->get_index(row_idx);
		const auto &entry = list_entries[mapped];

		auto &child = ListVector::GetEntry(arg_vec);
		UnifiedVectorFormat child_data;
		child.ToUnifiedFormat(ListVector::GetListSize(arg_vec), child_data);
		auto child_strings = UnifiedVectorFormat::GetData<string_t>(child_data);

		out.reserve(entry.length);
		for (idx_t k = 0; k < entry.length; k++) {
			auto child_idx = child_data.sel->get_index(entry.offset + k);
			if (!child_data.validity.RowIsValid(child_idx)) {
				continue;
			}
			const auto &s = child_strings[child_idx];
			if (s.GetSize() == 0) {
				continue;
			}
			out.emplace_back(s.GetData(), s.GetSize());
		}
	} else {
		auto strings = UnifiedVectorFormat::GetData<string_t>(arg_data);
		auto mapped = arg_data.sel->get_index(row_idx);
		const auto &s = strings[mapped];
		if (s.GetSize() > 0) {
			out.emplace_back(s.GetData(), s.GetSize());
		}
	}
	return out;
}

// Shared trim_adapters execution. arg layout (the function set picks):
//   args[0] seq VARCHAR
//   args[1] qual LIST(UTINYINT)
//   args[2] adapter VARCHAR or LIST(VARCHAR)
//   args[3..5] (optional 6-arg form): match_revcomp BOOLEAN, min_match INTEGER, allow_pre_start BOOLEAN
static void TrimAdaptersExecuteImpl(DataChunk &args, Vector &result, bool adapter_is_list) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 6;

	UnifiedVectorFormat seq_data, qual_data, adapter_data, revcomp_data, minmatch_data, prestart_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	args.data[2].ToUnifiedFormat(row_count, adapter_data);
	if (has_explicit_params) {
		args.data[3].ToUnifiedFormat(row_count, revcomp_data);
		args.data[4].ToUnifiedFormat(row_count, minmatch_data);
		args.data[5].ToUnifiedFormat(row_count, prestart_data);
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto revcomp_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<bool>(revcomp_data) : nullptr;
	auto minmatch_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(minmatch_data) : nullptr;
	auto prestart_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<bool>(prestart_data) : nullptr;

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
		auto ai = adapter_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi) &&
		                 adapter_data.validity.RowIsValid(ai);
		if (has_explicit_params) {
			all_valid = all_valid && revcomp_data.validity.RowIsValid(revcomp_data.sel->get_index(i)) &&
			            minmatch_data.validity.RowIsValid(minmatch_data.sel->get_index(i)) &&
			            prestart_data.validity.RowIsValid(prestart_data.sel->get_index(i));
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
			throw InvalidInputException("trim_adapters: sequence length (%llu) does not match quality length (%llu)",
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		auto adapters = ExtractAdapters(args.data[2], adapter_data, i, adapter_is_list);

		bool match_revcomp = false;
		int32_t min_match_param = 0;
		bool allow_pre_start = false;
		if (has_explicit_params) {
			match_revcomp = revcomp_ptr[revcomp_data.sel->get_index(i)];
			min_match_param = minmatch_ptr[minmatch_data.sel->get_index(i)];
			allow_pre_start = prestart_ptr[prestart_data.sel->get_index(i)];
		}
		if (min_match_param < 0) {
			throw InvalidInputException("trim_adapters: min_match must be >= 0 (got %d; 0 means use default)",
			                            min_match_param);
		}

		// Build full candidate list with optional RCs first, then scale
		// min_match against the total candidate count so RC-enabled searches
		// get the same stringency as if the user had passed adapters+RCs by
		// hand.
		std::vector<std::string> candidates;
		candidates.reserve(adapters.size() * (match_revcomp ? 2 : 1));
		for (const auto &a : adapters) {
			candidates.push_back(a);
		}
		if (match_revcomp) {
			for (const auto &a : adapters) {
				candidates.push_back(dna_revcomp_adapter(a));
			}
		}
		const std::size_t min_match =
		    min_match_param > 0 ? static_cast<std::size_t>(min_match_param) : default_min_match(candidates.size());

		// Run all candidates; take the leftmost trim_start across all matches.
		miint::qc::TrimResult tr {0, seq.GetSize()};
		const auto seq_bytes = reinterpret_cast<const std::uint8_t *>(seq.GetData());
		std::size_t best_trim_start = seq.GetSize();
		for (const auto &cand : candidates) {
			if (cand.size() < min_match) {
				continue;
			}
			auto m = miint::qc::AdapterMatcher::find(seq_bytes, seq.GetSize(),
			                                         reinterpret_cast<const std::uint8_t *>(cand.data()), cand.size(),
			                                         min_match, allow_pre_start);
			if (m.matched && m.trim_start < best_trim_start) {
				best_trim_start = m.trim_start;
			}
		}
		if (best_trim_start < seq.GetSize()) {
			tr.end = best_trim_start;
		}

		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
	}
}

static void TrimAdaptersVarcharExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimAdaptersExecuteImpl(args, result, /*adapter_is_list=*/false);
}

static void TrimAdaptersListExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimAdaptersExecuteImpl(args, result, /*adapter_is_list=*/true);
}

// ---------------------------------------------------------------------------
// filter_read
// ---------------------------------------------------------------------------

static constexpr int32_t DEFAULT_FILTER_MIN_LENGTH = 15;
static constexpr int32_t DEFAULT_FILTER_MAX_LENGTH = 0; // 0 means "off"
static constexpr int32_t DEFAULT_FILTER_QUALIFIED_Q = 15;
static constexpr int32_t DEFAULT_FILTER_MAX_UNQUAL_PCT = 40;
static constexpr int32_t DEFAULT_FILTER_MAX_N = 5;
static constexpr int32_t DEFAULT_FILTER_MIN_AVG_Q = 0; // 0 means "off"

static LogicalType FilterResultStructType() {
	return LogicalType::STRUCT({{"passed", LogicalType::BOOLEAN},
	                            {"fail_reason", LogicalType::VARCHAR},
	                            {"length", LogicalType::UINTEGER},
	                            {"n_bases", LogicalType::UINTEGER},
	                            {"low_qual_bases", LogicalType::UINTEGER},
	                            {"mean_quality", LogicalType::FLOAT}});
}

static void FilterReadExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 8;

	UnifiedVectorFormat seq_data, qual_data, p1_data, p2_data, p3_data, p4_data, p5_data, p6_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	if (has_explicit_params) {
		args.data[2].ToUnifiedFormat(row_count, p1_data);
		args.data[3].ToUnifiedFormat(row_count, p2_data);
		args.data[4].ToUnifiedFormat(row_count, p3_data);
		args.data[5].ToUnifiedFormat(row_count, p4_data);
		args.data[6].ToUnifiedFormat(row_count, p5_data);
		args.data[7].ToUnifiedFormat(row_count, p6_data);
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto p1_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p1_data) : nullptr;
	auto p2_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p2_data) : nullptr;
	auto p3_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p3_data) : nullptr;
	auto p4_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p4_data) : nullptr;
	auto p5_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p5_data) : nullptr;
	auto p6_ptr = has_explicit_params ? UnifiedVectorFormat::GetData<int32_t>(p6_data) : nullptr;

	auto &entries = StructVector::GetEntries(result);
	auto passed_data = FlatVector::GetData<bool>(*entries[0]);
	auto &fail_reason_vec = *entries[1];
	auto length_data = FlatVector::GetData<uint32_t>(*entries[2]);
	auto n_bases_data = FlatVector::GetData<uint32_t>(*entries[3]);
	auto low_qual_data = FlatVector::GetData<uint32_t>(*entries[4]);
	auto mean_q_data = FlatVector::GetData<float>(*entries[5]);
	auto &result_validity = FlatVector::Validity(result);
	auto &fail_reason_validity = FlatVector::Validity(fail_reason_vec);

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi);
		if (has_explicit_params) {
			all_valid = all_valid && p1_data.validity.RowIsValid(p1_data.sel->get_index(i)) &&
			            p2_data.validity.RowIsValid(p2_data.sel->get_index(i)) &&
			            p3_data.validity.RowIsValid(p3_data.sel->get_index(i)) &&
			            p4_data.validity.RowIsValid(p4_data.sel->get_index(i)) &&
			            p5_data.validity.RowIsValid(p5_data.sel->get_index(i)) &&
			            p6_data.validity.RowIsValid(p6_data.sel->get_index(i));
		}
		if (!all_valid) {
			result_validity.SetInvalid(i);
			continue;
		}

		const auto &seq = seq_ptr[si];
		const uint8_t *qptr;
		idx_t qlen;
		GetQualListSlice(args.data[1], qual_data, i, qptr, qlen);
		if (seq.GetSize() != qlen) {
			throw InvalidInputException("filter_read: sequence length (%llu) does not match quality length (%llu)",
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		int32_t min_length = DEFAULT_FILTER_MIN_LENGTH;
		int32_t max_length = DEFAULT_FILTER_MAX_LENGTH;
		int32_t qualified_q = DEFAULT_FILTER_QUALIFIED_Q;
		int32_t max_unqual_pct = DEFAULT_FILTER_MAX_UNQUAL_PCT;
		int32_t max_n = DEFAULT_FILTER_MAX_N;
		int32_t min_avg_q = DEFAULT_FILTER_MIN_AVG_Q;
		if (has_explicit_params) {
			min_length = p1_ptr[p1_data.sel->get_index(i)];
			max_length = p2_ptr[p2_data.sel->get_index(i)];
			qualified_q = p3_ptr[p3_data.sel->get_index(i)];
			max_unqual_pct = p4_ptr[p4_data.sel->get_index(i)];
			max_n = p5_ptr[p5_data.sel->get_index(i)];
			min_avg_q = p6_ptr[p6_data.sel->get_index(i)];
		}
		if (min_length < 0 || max_length < 0 || qualified_q < 0 || qualified_q > 93 || max_unqual_pct < 0 ||
		    max_unqual_pct > 100 || max_n < 0 || min_avg_q < 0 || min_avg_q > 93) {
			throw InvalidInputException(
			    "filter_read: parameter out of range (min_length>=0, max_length>=0, qualified_q in 0..93, "
			    "max_unqualified_pct in 0..100, max_n>=0, min_avg_q in 0..93)");
		}

		// Empty seq is an immediate length failure — skip the metric pass.
		if (seq.GetSize() == 0) {
			passed_data[i] = false;
			FlatVector::GetData<string_t>(fail_reason_vec)[i] = StringVector::AddString(fail_reason_vec, "length");
			length_data[i] = 0;
			n_bases_data[i] = 0;
			low_qual_data[i] = 0;
			mean_q_data[i] = 0.0f;
			continue;
		}

		const auto seq_bytes = reinterpret_cast<const std::uint8_t *>(seq.GetData());
		const auto m =
		    miint::qc::ReadFilter::measure(seq_bytes, qptr, seq.GetSize(), static_cast<std::uint8_t>(qualified_q));
		const float mean_q = static_cast<float>(m.qual_sum) / static_cast<float>(m.length);

		// Apply fastp's fail-reason precedence: low_qual% → avg_qual →
		// n_base → length_min → length_max. First failing check wins.
		const char *fail_reason = nullptr;
		// low_qual %: strictly > limit fails (matches fastp's `>` comparison)
		if (static_cast<double>(m.low_qual_bases) * 100.0 >
		    static_cast<double>(max_unqual_pct) * static_cast<double>(m.length)) {
			fail_reason = "quality";
		} else if (min_avg_q > 0 && mean_q < static_cast<float>(min_avg_q)) {
			fail_reason = "quality";
		} else if (static_cast<int32_t>(m.n_bases) > max_n) {
			fail_reason = "n_base";
		} else if (static_cast<int32_t>(m.length) < min_length) {
			fail_reason = "length";
		} else if (max_length > 0 && static_cast<int32_t>(m.length) > max_length) {
			fail_reason = "too_long";
		}

		passed_data[i] = (fail_reason == nullptr);
		if (fail_reason != nullptr) {
			FlatVector::GetData<string_t>(fail_reason_vec)[i] = StringVector::AddString(fail_reason_vec, fail_reason);
		} else {
			fail_reason_validity.SetInvalid(i);
		}
		length_data[i] = m.length;
		n_bases_data[i] = m.n_bases;
		low_qual_data[i] = m.low_qual_bases;
		mean_q_data[i] = mean_q;
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

	// filter_read: 2-arg (defaults) + 8-arg (min_length, max_length, qualified_q,
	// max_unqualified_pct, max_n, min_avg_q).
	{
		ScalarFunctionSet set("filter_read");
		const auto qual_t = LogicalType::LIST(LogicalType::UTINYINT);
		const auto i = LogicalType::INTEGER;

		ScalarFunction two_arg("filter_read", {LogicalType::VARCHAR, qual_t}, FilterResultStructType(),
		                       FilterReadExecute);
		two_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(two_arg);

		ScalarFunction eight_arg("filter_read", {LogicalType::VARCHAR, qual_t, i, i, i, i, i, i},
		                         FilterResultStructType(), FilterReadExecute);
		eight_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(eight_arg);

		loader.RegisterFunction(set);
	}

	// trim_adapters: 4 overloads — {VARCHAR, LIST(VARCHAR)} adapter × {3-arg defaults,
	// 6-arg full (match_revcomp, min_match, allow_pre_start)}.
	{
		ScalarFunctionSet set("trim_adapters");
		const auto qual_t = LogicalType::LIST(LogicalType::UTINYINT);
		const auto adapter_list_t = LogicalType::LIST(LogicalType::VARCHAR);

		ScalarFunction varchar_3arg("trim_adapters", {LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR},
		                            TrimResultStructType(), TrimAdaptersVarcharExecute);
		varchar_3arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(varchar_3arg);

		ScalarFunction list_3arg("trim_adapters", {LogicalType::VARCHAR, qual_t, adapter_list_t},
		                         TrimResultStructType(), TrimAdaptersListExecute);
		list_3arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(list_3arg);

		ScalarFunction varchar_6arg("trim_adapters",
		                            {LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, LogicalType::BOOLEAN,
		                             LogicalType::INTEGER, LogicalType::BOOLEAN},
		                            TrimResultStructType(), TrimAdaptersVarcharExecute);
		varchar_6arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(varchar_6arg);

		ScalarFunction list_6arg("trim_adapters",
		                         {LogicalType::VARCHAR, qual_t, adapter_list_t, LogicalType::BOOLEAN,
		                          LogicalType::INTEGER, LogicalType::BOOLEAN},
		                         TrimResultStructType(), TrimAdaptersListExecute);
		list_6arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		set.AddFunction(list_6arg);

		loader.RegisterFunction(set);
	}
}

} // namespace duckdb
