#include "qc_functions.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/scalar_function.hpp"

#include "qc_algorithms.hpp"
#include "sequence_utils.hpp"
#include "table_function_common.hpp"

#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

namespace duckdb {

static constexpr const char *QC_VERSION_STRING = "qc 0.1.0 (port of fastp algorithms)";

static void QcVersionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto out = ConstantVector::GetData<string_t>(result);
	out[0] = StringVector::AddString(result, QC_VERSION_STRING);
}

// ---------------------------------------------------------------------------
// Shared types and helpers for all trim_* scalars
// ---------------------------------------------------------------------------

// Phred quality scores top out at 93 (the highest value representable in the
// ASCII-33 FASTQ encoding). Used at validation boundaries.
static constexpr int32_t MAX_PHRED_QUALITY = 93;

static constexpr int32_t DEFAULT_QUAL_WINDOW = 4;
static constexpr int32_t DEFAULT_QUAL_MEAN = 20;
// fastp's documented window-size range is 1..1000. The upper bound also keeps
// threshold_sum (mean_quality * window_size) safely within uint32_t.
static constexpr int32_t MAX_QUAL_WINDOW = 1000;

static constexpr int32_t DEFAULT_POLY_MIN_LEN = 10;
static constexpr int32_t DEFAULT_POLY_MAX_MISMATCH = 5;
// Default quality-aware gate for polyG: trim region's mean Phred must be <=
// this threshold. Pass max_window_mean_q=QUAL_GATE_DISABLED to make the gate
// a no-op since any real Phred score satisfies <= 93.
static constexpr int32_t DEFAULT_POLYG_MAX_WINDOW_MEAN_Q = 5;
// Sentinel for "disable the polyG quality gate" — happens to coincide with
// MAX_PHRED_QUALITY since "any real Phred score" is the threshold that makes
// the gate vacuously true. Kept as a separate named constant because the
// design intent (a disable sentinel) is distinct from the physical maximum.
static constexpr int32_t QUAL_GATE_DISABLED = MAX_PHRED_QUALITY;

static LogicalType TrimResultStructType() {
	return LogicalType::STRUCT({{"sequence", LogicalType::VARCHAR},
	                            {"quality", LogicalType::LIST(LogicalType::UTINYINT)},
	                            {"trimmed_5p", LogicalType::UINTEGER},
	                            {"trimmed_3p", LogicalType::UINTEGER}});
}

// Reserve the trimmed-quality child buffer up front: the trimmed output is at
// most as large as the input qual buffer, so one reserve before the row loop
// eliminates O(log n) reallocations across the chunk.
static void ReserveTrimOutputChildBuffer(DataChunk &args, Vector &qual_out_vec, idx_t qual_child_offset) {
	ListVector::Reserve(qual_out_vec, qual_child_offset + ListVector::GetListSize(args.data[1]));
}

// Write one row of the trim_result struct: trimmed sequence, trimmed quality
// (as a sliced LIST(UTINYINT)), and the trimmed_5p/trimmed_3p counts. The
// caller must have reserved enough child-buffer capacity (see
// ReserveTrimOutputChildBuffer) — this function does no per-row reserve.
static void WriteTrimRow(idx_t i, const string_t &seq, const uint8_t *qptr, idx_t qlen, miint::qc::TrimResult tr,
                         Vector &seq_out_vec, Vector &qual_out_vec, list_entry_t *qual_out_entries,
                         idx_t &qual_child_offset, uint32_t *trimmed_5p_data, uint32_t *trimmed_3p_data) {
	const idx_t kept_len = tr.end - tr.start;
	FlatVector::GetData<string_t>(seq_out_vec)[i] =
	    StringVector::AddString(seq_out_vec, seq.GetData() + tr.start, kept_len);

	auto &qual_child = ListVector::GetEntry(qual_out_vec);
	auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
	std::memcpy(qual_child_data + qual_child_offset, qptr + tr.start, kept_len);
	qual_out_entries[i].offset = qual_child_offset;
	qual_out_entries[i].length = kept_len;
	qual_child_offset += kept_len;
	ListVector::SetListSize(qual_out_vec, qual_child_offset);

	trimmed_5p_data[i] = static_cast<uint32_t>(tr.start);
	trimmed_3p_data[i] = static_cast<uint32_t>(qlen - tr.end);
}

// Shared driver for every per-row trim_* scalar: handles arg unification,
// per-row validity, seq/qual extraction and length check, optional-param
// readout, and the WriteTrimRow call. The trim algorithm itself (and its
// param validation) is provided as a callable.
//
// MakeTrim signature:
//   miint::qc::TrimResult (const string_t &seq, const uint8_t *qptr, idx_t qlen,
//                         bool has_explicit_params, const int32_t opt_params[3],
//                         const char *fn_name)
template <typename MakeTrim>
static void TrimExecuteImpl(DataChunk &args, Vector &result, idx_t n_optional, const char *fn_name,
                            MakeTrim make_trim) {
	const idx_t row_count = args.size();
	const bool has_explicit_params = args.ColumnCount() == 2 + n_optional;
	D_ASSERT(n_optional <= 3);

	UnifiedVectorFormat seq_data, qual_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);

	UnifiedVectorFormat opt_data[3];
	const int32_t *opt_ptr[3] = {nullptr, nullptr, nullptr};
	if (has_explicit_params) {
		for (idx_t k = 0; k < n_optional; k++) {
			args.data[2 + k].ToUnifiedFormat(row_count, opt_data[k]);
			opt_ptr[k] = UnifiedVectorFormat::GetData<int32_t>(opt_data[k]);
		}
	}

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);

	auto &entries = StructVector::GetEntries(result);
	auto &seq_out_vec = *entries[0];
	auto &qual_out_vec = *entries[1];
	auto trimmed_5p_data = FlatVector::GetData<uint32_t>(*entries[2]);
	auto trimmed_3p_data = FlatVector::GetData<uint32_t>(*entries[3]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual_out_entries = FlatVector::GetData<list_entry_t>(qual_out_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_out_vec);
	ReserveTrimOutputChildBuffer(args, qual_out_vec, qual_child_offset);

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi);
		if (has_explicit_params) {
			for (idx_t k = 0; k < n_optional && all_valid; k++) {
				all_valid = opt_data[k].validity.RowIsValid(opt_data[k].sel->get_index(i));
			}
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
		GetListUInt8Slice(args.data[1], qual_data, i, qptr, qlen);

		if (seq.GetSize() != qlen) {
			throw InvalidInputException("%s: sequence length (%llu) does not match quality length (%llu)", fn_name,
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		int32_t p[3] = {0, 0, 0};
		if (has_explicit_params) {
			for (idx_t k = 0; k < n_optional; k++) {
				p[k] = opt_ptr[k][opt_data[k].sel->get_index(i)];
			}
		}

		auto tr = make_trim(seq, qptr, qlen, has_explicit_params, p, fn_name);
		WriteTrimRow(i, seq, qptr, qlen, tr, seq_out_vec, qual_out_vec, qual_out_entries, qual_child_offset,
		             trimmed_5p_data, trimmed_3p_data);
	}
}

// Validate window+mean params (shared by all three sliding-window flavors).
// After this returns, the trim algorithm's own require_window check is
// statically unreachable, so the per-row try/catch is unnecessary.
static void ValidateQualityWindowParams(int32_t window, int32_t mean, const char *fn_name) {
	if (window <= 0 || window > MAX_QUAL_WINDOW) {
		throw InvalidInputException("%s: window_size must be in 1..%d (got %d)", fn_name, MAX_QUAL_WINDOW, window);
	}
	if (mean < 0 || mean > MAX_PHRED_QUALITY) {
		throw InvalidInputException("%s: mean_quality must be in 0..%d (got %d)", fn_name, MAX_PHRED_QUALITY, mean);
	}
}

using SlidingTrimFn = miint::qc::TrimResult (*)(const std::uint8_t *, std::size_t, std::size_t, std::uint8_t);

static void RunTrimQuality(DataChunk &args, Vector &result, SlidingTrimFn trim_fn, const char *fn_name) {
	TrimExecuteImpl(args, result, /*n_optional=*/2, fn_name,
	                [trim_fn](const string_t &, const uint8_t *qptr, idx_t qlen, bool has_explicit_params,
	                          const int32_t p[3], const char *fname) {
		                int32_t window = has_explicit_params ? p[0] : DEFAULT_QUAL_WINDOW;
		                int32_t mean = has_explicit_params ? p[1] : DEFAULT_QUAL_MEAN;
		                ValidateQualityWindowParams(window, mean, fname);
		                return trim_fn(qptr, qlen, static_cast<std::size_t>(window), static_cast<std::uint8_t>(mean));
	                });
}

static void TrimQuality3pFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	RunTrimQuality(args, result, &miint::qc::SlidingWindowTrimmer::trim_3p, "trim_quality_3p");
}

static void TrimQuality5pFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	RunTrimQuality(args, result, &miint::qc::SlidingWindowTrimmer::trim_5p, "trim_quality_5p");
}

static void TrimQualitySlidingFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	RunTrimQuality(args, result, &miint::qc::SlidingWindowTrimmer::trim_sliding, "trim_quality_sliding");
}

static void TrimPolygExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimExecuteImpl(
	    args, result, /*n_optional=*/3, "trim_polyg",
	    [](const string_t &seq, const uint8_t *qptr, idx_t qlen, bool has_explicit_params, const int32_t p[3],
	       const char *fname) {
		    int32_t min_len = has_explicit_params ? p[0] : DEFAULT_POLY_MIN_LEN;
		    int32_t max_mismatch = has_explicit_params ? p[1] : DEFAULT_POLY_MAX_MISMATCH;
		    int32_t max_window_q = has_explicit_params ? p[2] : DEFAULT_POLYG_MAX_WINDOW_MEAN_Q;
		    if (min_len < 1) {
			    throw InvalidInputException("%s: min_len must be >= 1 (got %d)", fname, min_len);
		    }
		    if (max_mismatch < 0) {
			    throw InvalidInputException("%s: max_mismatch must be >= 0 (got %d)", fname, max_mismatch);
		    }
		    if (max_window_q < 0 || max_window_q > QUAL_GATE_DISABLED) {
			    throw InvalidInputException("%s: max_window_mean_q must be 0..%d (got %d); pass %d to disable", fname,
			                                QUAL_GATE_DISABLED, max_window_q, QUAL_GATE_DISABLED);
		    }
		    return miint::qc::PolyXScanner::scan_polyg(
		        reinterpret_cast<const std::uint8_t *>(seq.GetData()), qptr, qlen, static_cast<std::size_t>(min_len),
		        static_cast<std::uint32_t>(max_mismatch), static_cast<std::uint8_t>(max_window_q));
	    });
}

static void TrimPolyxExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	TrimExecuteImpl(args, result, /*n_optional=*/2, "trim_polyx",
	                [](const string_t &seq, const uint8_t *, idx_t, bool has_explicit_params, const int32_t p[3],
	                   const char *fname) {
		                int32_t min_len = has_explicit_params ? p[0] : DEFAULT_POLY_MIN_LEN;
		                int32_t max_mismatch = has_explicit_params ? p[1] : DEFAULT_POLY_MAX_MISMATCH;
		                if (min_len < 1) {
			                throw InvalidInputException("%s: min_len must be >= 1 (got %d)", fname, min_len);
		                }
		                if (max_mismatch < 0) {
			                throw InvalidInputException("%s: max_mismatch must be >= 0 (got %d)", fname, max_mismatch);
		                }
		                return miint::qc::PolyXScanner::scan_polyx(
		                    reinterpret_cast<const std::uint8_t *>(seq.GetData()), seq.GetSize(),
		                    static_cast<std::size_t>(min_len), static_cast<std::uint32_t>(max_mismatch));
	                });
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
// are silently skipped). For the LIST case the caller passes the already-
// unified child format so this function does no per-row ToUnifiedFormat call.
static std::vector<std::string> ExtractAdapters(Vector &arg_vec, UnifiedVectorFormat &arg_data,
                                                const UnifiedVectorFormat *list_child_data, idx_t row_idx,
                                                bool is_list) {
	std::vector<std::string> out;
	if (is_list) {
		auto list_entries = UnifiedVectorFormat::GetData<list_entry_t>(arg_data);
		auto mapped = arg_data.sel->get_index(row_idx);
		const auto &entry = list_entries[mapped];

		auto child_strings = UnifiedVectorFormat::GetData<string_t>(*list_child_data);

		out.reserve(entry.length);
		for (idx_t k = 0; k < entry.length; k++) {
			auto child_idx = list_child_data->sel->get_index(entry.offset + k);
			if (!list_child_data->validity.RowIsValid(child_idx)) {
				continue;
			}
			const auto &s = child_strings[child_idx];
			if (s.GetSize() == 0) {
				continue;
			}
			out.emplace_back(s.GetData(), s.GetSize());
		}
	} else {
		(void)arg_vec; // VARCHAR overload reads directly from arg_data
		auto strings = UnifiedVectorFormat::GetData<string_t>(arg_data);
		auto mapped = arg_data.sel->get_index(row_idx);
		const auto &s = strings[mapped];
		if (s.GetSize() > 0) {
			out.emplace_back(s.GetData(), s.GetSize());
		}
	}
	return out;
}

// Build the full candidate list (with optional RC extension) and resolve the
// effective min_match. min_match scales against the total candidate count so
// RC-enabled searches match the stringency the user would get by passing
// adapters+RCs by hand.
static void BuildAdapterCandidates(std::vector<std::string> adapters, bool match_revcomp, int32_t min_match_param,
                                   std::vector<std::string> &out_candidates, std::size_t &out_min_match) {
	if (match_revcomp) {
		out_candidates.reserve(adapters.size() * 2);
		for (const auto &a : adapters) {
			out_candidates.push_back(a);
		}
		for (const auto &a : adapters) {
			out_candidates.push_back(dna_revcomp_adapter(a));
		}
	} else {
		out_candidates = std::move(adapters);
	}
	out_min_match =
	    min_match_param > 0 ? static_cast<std::size_t>(min_match_param) : default_min_match(out_candidates.size());
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

	// The LIST(VARCHAR) child vector is independent of row — unify it once
	// per chunk, not per row.
	UnifiedVectorFormat list_child_data;
	if (adapter_is_list) {
		auto &child = ListVector::GetEntry(args.data[2]);
		child.ToUnifiedFormat(ListVector::GetListSize(args.data[2]), list_child_data);
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
	ReserveTrimOutputChildBuffer(args, qual_out_vec, qual_child_offset);

	// Hoist adapter preprocessing for the common case where the adapter
	// argument (and all named params) are constant across the chunk — e.g.
	// `trim_adapters(seq, qual, 'AGATCGGAAGAGC')`. Without this, every row
	// would rebuild the same candidate list, allocating per-row std::vector +
	// per-row std::string copies (plus a per-row RC computation under
	// match_revcomp).
	const bool adapter_constant = args.data[2].GetVectorType() == VectorType::CONSTANT_VECTOR;
	const bool params_all_constant =
	    !has_explicit_params || (args.data[3].GetVectorType() == VectorType::CONSTANT_VECTOR &&
	                             args.data[4].GetVectorType() == VectorType::CONSTANT_VECTOR &&
	                             args.data[5].GetVectorType() == VectorType::CONSTANT_VECTOR);
	std::vector<std::string> hoisted_candidates;
	std::size_t hoisted_min_match = 0;
	bool hoisted_allow_pre_start = false;
	bool hoisted_null = false; // a constant value is itself NULL → all rows produce a NULL result
	const bool hoisted = adapter_constant && params_all_constant;
	if (hoisted) {
		bool any_null = !adapter_data.validity.RowIsValid(adapter_data.sel->get_index(0));
		if (has_explicit_params && !any_null) {
			any_null = !revcomp_data.validity.RowIsValid(revcomp_data.sel->get_index(0)) ||
			           !minmatch_data.validity.RowIsValid(minmatch_data.sel->get_index(0)) ||
			           !prestart_data.validity.RowIsValid(prestart_data.sel->get_index(0));
		}
		if (any_null) {
			hoisted_null = true;
		} else {
			auto adapters = ExtractAdapters(args.data[2], adapter_data, adapter_is_list ? &list_child_data : nullptr, 0,
			                                adapter_is_list);
			const bool mr = has_explicit_params ? revcomp_ptr[revcomp_data.sel->get_index(0)] : false;
			const int32_t mm_param = has_explicit_params ? minmatch_ptr[minmatch_data.sel->get_index(0)] : 0;
			hoisted_allow_pre_start = has_explicit_params ? prestart_ptr[prestart_data.sel->get_index(0)] : false;
			if (mm_param < 0) {
				throw InvalidInputException("trim_adapters: min_match must be >= 0 (got %d; 0 means use default)",
				                            mm_param);
			}
			BuildAdapterCandidates(std::move(adapters), mr, mm_param, hoisted_candidates, hoisted_min_match);
		}
	}

	for (idx_t i = 0; i < row_count; i++) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi);
		if (hoisted) {
			all_valid = all_valid && !hoisted_null;
		} else {
			all_valid = all_valid && adapter_data.validity.RowIsValid(adapter_data.sel->get_index(i));
			if (has_explicit_params) {
				all_valid = all_valid && revcomp_data.validity.RowIsValid(revcomp_data.sel->get_index(i)) &&
				            minmatch_data.validity.RowIsValid(minmatch_data.sel->get_index(i)) &&
				            prestart_data.validity.RowIsValid(prestart_data.sel->get_index(i));
			}
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
		GetListUInt8Slice(args.data[1], qual_data, i, qptr, qlen);
		if (seq.GetSize() != qlen) {
			throw InvalidInputException("trim_adapters: sequence length (%llu) does not match quality length (%llu)",
			                            (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		const std::vector<std::string> *candidates_ptr;
		std::size_t min_match;
		bool allow_pre_start;
		std::vector<std::string> row_candidates; // lifetime extends through end of row iteration
		if (hoisted) {
			candidates_ptr = &hoisted_candidates;
			min_match = hoisted_min_match;
			allow_pre_start = hoisted_allow_pre_start;
		} else {
			auto adapters = ExtractAdapters(args.data[2], adapter_data, adapter_is_list ? &list_child_data : nullptr, i,
			                                adapter_is_list);
			const bool mr = has_explicit_params ? revcomp_ptr[revcomp_data.sel->get_index(i)] : false;
			const int32_t mm_param = has_explicit_params ? minmatch_ptr[minmatch_data.sel->get_index(i)] : 0;
			allow_pre_start = has_explicit_params ? prestart_ptr[prestart_data.sel->get_index(i)] : false;
			if (mm_param < 0) {
				throw InvalidInputException("trim_adapters: min_match must be >= 0 (got %d; 0 means use default)",
				                            mm_param);
			}
			BuildAdapterCandidates(std::move(adapters), mr, mm_param, row_candidates, min_match);
			candidates_ptr = &row_candidates;
		}

		// Leftmost trim_start across all candidates wins.
		miint::qc::TrimResult tr {0, seq.GetSize()};
		const auto seq_bytes = reinterpret_cast<const std::uint8_t *>(seq.GetData());
		std::size_t best_trim_start = seq.GetSize();
		for (const auto &cand : *candidates_ptr) {
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

// fail_reason values returned by filter_read. SQL callers compare against
// these directly (e.g., WHERE (f).fail_reason = 'quality'), so they are part
// of the user-facing contract; centralizing them here avoids drift between
// code paths and tests.
static constexpr const char *FAIL_REASON_QUALITY = "quality";
static constexpr const char *FAIL_REASON_N_BASE = "n_base";
static constexpr const char *FAIL_REASON_LENGTH = "length";
static constexpr const char *FAIL_REASON_TOO_LONG = "too_long";

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
		GetListUInt8Slice(args.data[1], qual_data, i, qptr, qlen);
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
		if (min_length < 0 || max_length < 0 || qualified_q < 0 || qualified_q > MAX_PHRED_QUALITY ||
		    max_unqual_pct < 0 || max_unqual_pct > 100 || max_n < 0 || min_avg_q < 0 || min_avg_q > MAX_PHRED_QUALITY) {
			throw InvalidInputException(
			    "filter_read: parameter out of range (min_length>=0, max_length>=0, qualified_q in 0..%d, "
			    "max_unqualified_pct in 0..100, max_n>=0, min_avg_q in 0..%d)",
			    MAX_PHRED_QUALITY, MAX_PHRED_QUALITY);
		}

		// Empty seq is an immediate length failure — skip the metric pass.
		// (The seq/qual length-mismatch check above already ran, so this is
		// reached only when qlen is also 0.)
		if (seq.GetSize() == 0) {
			passed_data[i] = false;
			FlatVector::GetData<string_t>(fail_reason_vec)[i] =
			    StringVector::AddString(fail_reason_vec, FAIL_REASON_LENGTH);
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
			fail_reason = FAIL_REASON_QUALITY;
		} else if (min_avg_q > 0 && mean_q < static_cast<float>(min_avg_q)) {
			fail_reason = FAIL_REASON_QUALITY;
		} else if (static_cast<int32_t>(m.n_bases) > max_n) {
			fail_reason = FAIL_REASON_N_BASE;
		} else if (static_cast<int32_t>(m.length) < min_length) {
			fail_reason = FAIL_REASON_LENGTH;
		} else if (max_length > 0 && static_cast<int32_t>(m.length) > max_length) {
			fail_reason = FAIL_REASON_TOO_LONG;
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
