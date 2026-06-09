#include "qc_functions.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"

#include "qc_algorithms.hpp"
#include "sequence_utils.hpp"
#include "table_function_common.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <unordered_set>
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

// Remove exact-duplicate candidates in place, preserving first-occurrence order.
// Production adapter sets are highly redundant (exact dups, and reverse-complement
// pairs once RC expansion runs), so collapsing them cuts matching cost without
// changing the leftmost-match result. Must run AFTER min_match is fixed (see
// BuildAdapterCandidates) so that list size still drives sensitivity.
static void DedupAdapterCandidates(std::vector<std::string> &candidates) {
	std::unordered_set<std::string> seen;
	seen.reserve(candidates.size());
	std::size_t write = 0;
	for (std::size_t read = 0; read < candidates.size(); read++) {
		if (seen.insert(candidates[read]).second) {
			if (write != read) {
				candidates[write] = std::move(candidates[read]);
			}
			write++;
		}
	}
	candidates.resize(write);
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
	// Collapse redundant candidates AFTER fixing min_match: dedup then only
	// removes wasted matching work, never changes the auto-scaled sensitivity.
	DedupAdapterCandidates(out_candidates);
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

		// Leftmost trim_start across all candidates wins. find_leftmost applies
		// the position-0 early-exit and returns the kept end directly (= seq size
		// when nothing matches).
		const auto seq_bytes = reinterpret_cast<const std::uint8_t *>(seq.GetData());
		const std::size_t kept_end = miint::qc::AdapterMatcher::find_leftmost(seq_bytes, seq.GetSize(), *candidates_ptr,
		                                                                      min_match, allow_pre_start);
		const miint::qc::TrimResult tr {0, kept_end};

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
// trim_adapters_pe
// ---------------------------------------------------------------------------

// fastp's default paired-end overlap-analysis parameters
// (--overlap_len_require / --overlap_diff_limit / --overlap_diff_percent_limit).
static constexpr int32_t DEFAULT_PE_OVERLAP_REQUIRE = 30;
static constexpr int32_t DEFAULT_PE_OVERLAP_DIFF_LIMIT = 5;
static constexpr int32_t DEFAULT_PE_OVERLAP_DIFF_PERCENT = 20;

// Bind-time configuration: the overlap parameters and the (RC-expanded, deduped)
// adapter-by-sequence fallback candidate list. The 4-arg form leaves `candidates`
// empty (overlap-only); the 11-arg form populates it (fastp step 8 -> step 9).
struct TrimAdaptersPeBindData : public FunctionData {
	int32_t overlap_require = DEFAULT_PE_OVERLAP_REQUIRE;
	int32_t overlap_diff_limit = DEFAULT_PE_OVERLAP_DIFF_LIMIT;
	int32_t overlap_diff_percent = DEFAULT_PE_OVERLAP_DIFF_PERCENT;
	std::vector<std::string> candidates; // empty => no adapter-by-sequence fallback
	std::size_t min_match = 0;
	bool allow_pre_start = false;

	unique_ptr<FunctionData> Copy() const override {
		auto c = make_uniq<TrimAdaptersPeBindData>();
		c->overlap_require = overlap_require;
		c->overlap_diff_limit = overlap_diff_limit;
		c->overlap_diff_percent = overlap_diff_percent;
		c->candidates = candidates;
		c->min_match = min_match;
		c->allow_pre_start = allow_pre_start;
		return std::move(c);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &o = other_p.Cast<TrimAdaptersPeBindData>();
		return overlap_require == o.overlap_require && overlap_diff_limit == o.overlap_diff_limit &&
		       overlap_diff_percent == o.overlap_diff_percent && min_match == o.min_match &&
		       allow_pre_start == o.allow_pre_start && candidates == o.candidates;
	}
};

// Per-thread copy of the bind-time config so Execute reads plain members without
// touching bind_info. The candidate list is small (protocol adapter set) and
// copied once per thread in InitLocalState, not per row.
struct TrimAdaptersPeLocalState : public FunctionLocalState {
	int32_t overlap_require;
	int32_t overlap_diff_limit;
	int32_t overlap_diff_percent;
	std::vector<std::string> candidates;
	std::size_t min_match;
	bool allow_pre_start;

	explicit TrimAdaptersPeLocalState(const TrimAdaptersPeBindData &bd)
	    : overlap_require(bd.overlap_require), overlap_diff_limit(bd.overlap_diff_limit),
	      overlap_diff_percent(bd.overlap_diff_percent), candidates(bd.candidates), min_match(bd.min_match),
	      allow_pre_start(bd.allow_pre_start) {
	}
};

static unique_ptr<FunctionLocalState>
TrimAdaptersPeInitLocalState(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data) {
	(void)state;
	(void)expr;
	return make_uniq<TrimAdaptersPeLocalState>(bind_data->Cast<TrimAdaptersPeBindData>());
}

static LogicalType TrimAdaptersPeResultStructType() {
	return LogicalType::STRUCT({{"sequence1", LogicalType::VARCHAR},
	                            {"quality1", LogicalType::LIST(LogicalType::UTINYINT)},
	                            {"sequence2", LogicalType::VARCHAR},
	                            {"quality2", LogicalType::LIST(LogicalType::UTINYINT)},
	                            // overlap_len is the overlap-analysis result (detected overlap width,
	                            // 0 if none), independent of adapter_trimmed and of the trim source: it
	                            // is nonzero for a full-insert overlap at offset >= 0 (adapter_trimmed
	                            // false), and in the 11-arg form a row trimmed by the adapter-by-sequence
	                            // fallback still reports the overlap width found by analysis.
	                            {"overlap_len", LogicalType::INTEGER},
	                            {"adapter_trimmed", LogicalType::BOOLEAN},
	                            {"trimmed1_3p", LogicalType::UINTEGER},
	                            {"trimmed2_3p", LogicalType::UINTEGER}});
}

// Write one trimmed mate into its (sequence, quality) output pair. The kept
// range is [0, keep_len) of both seq and qual (overlap trimming is 3'-only).
// The caller must have reserved enough quality child-buffer capacity up front.
static void WritePeMate(idx_t i, const string_t &seq, const uint8_t *qptr, idx_t keep_len, Vector &seq_out_vec,
                        Vector &qual_out_vec, list_entry_t *qual_entries, uint8_t *qual_child,
                        idx_t &qual_child_offset) {
	FlatVector::GetData<string_t>(seq_out_vec)[i] = StringVector::AddString(seq_out_vec, seq.GetData(), keep_len);
	std::memcpy(qual_child + qual_child_offset, qptr, keep_len);
	qual_entries[i].offset = qual_child_offset;
	qual_entries[i].length = keep_len;
	qual_child_offset += keep_len;
	ListVector::SetListSize(qual_out_vec, qual_child_offset);
}

// 4-arg form: overlap-only with fastp defaults; no adapter-by-sequence fallback.
static unique_ptr<FunctionData> TrimAdaptersPeBind4(ClientContext &ctx, ScalarFunction &fn,
                                                    vector<unique_ptr<Expression>> &arguments) {
	(void)ctx;
	(void)fn;
	(void)arguments;
	return make_uniq<TrimAdaptersPeBindData>(); // all defaults; candidates empty
}

// 11-arg form: trim_adapters_pe(seq1, qual1, seq2, qual2, adapters, overlap_require,
// overlap_diff_limit, overlap_diff_percent_limit, match_revcomp, min_match,
// allow_pre_start). The adapter list + tuning params must be constant; they are
// evaluated once here and the fallback candidate set (RC-expanded, deduped) is
// built with min_match fixed from the pre-dedup count.
static unique_ptr<FunctionData> TrimAdaptersPeBind11(ClientContext &ctx, ScalarFunction &fn,
                                                     vector<unique_ptr<Expression>> &arguments) {
	(void)fn;
	for (idx_t k = 4; k < 11; k++) {
		if (!arguments[k]->IsFoldable()) {
			throw InvalidInputException("trim_adapters_pe: the adapter list and tuning parameters must be constant "
			                            "values, not column references");
		}
	}
	// Evaluate each constant param, rejecting NULL with a clean error (GetValue on
	// a NULL Value would otherwise throw an opaque InternalException).
	auto eval_int = [&](idx_t idx, const char *name) {
		const Value v = ExpressionExecutor::EvaluateScalar(ctx, *arguments[idx]);
		if (v.IsNull()) {
			throw InvalidInputException("trim_adapters_pe: %s must not be NULL", name);
		}
		return v.GetValue<int32_t>();
	};
	auto eval_bool = [&](idx_t idx, const char *name) {
		const Value v = ExpressionExecutor::EvaluateScalar(ctx, *arguments[idx]);
		if (v.IsNull()) {
			throw InvalidInputException("trim_adapters_pe: %s must not be NULL", name);
		}
		return v.GetValue<bool>();
	};
	const int32_t overlap_require = eval_int(5, "overlap_require");
	const int32_t overlap_diff_limit = eval_int(6, "overlap_diff_limit");
	const int32_t overlap_diff_percent = eval_int(7, "overlap_diff_percent_limit");
	const bool match_revcomp = eval_bool(8, "match_revcomp");
	const int32_t min_match_param = eval_int(9, "min_match");
	const bool allow_pre_start = eval_bool(10, "allow_pre_start");

	if (overlap_require < 1) {
		throw InvalidInputException("trim_adapters_pe: overlap_require must be >= 1 (got %d)", overlap_require);
	}
	if (overlap_diff_limit < 0) {
		throw InvalidInputException("trim_adapters_pe: overlap_diff_limit must be >= 0 (got %d)", overlap_diff_limit);
	}
	if (overlap_diff_percent < 0 || overlap_diff_percent > 100) {
		throw InvalidInputException("trim_adapters_pe: overlap_diff_percent_limit must be 0..100 (got %d)",
		                            overlap_diff_percent);
	}
	if (min_match_param < 0) {
		throw InvalidInputException("trim_adapters_pe: min_match must be >= 0 (got %d; 0 means use default)",
		                            min_match_param);
	}

	// Extract the adapter list (LIST(VARCHAR); NULL list or NULL/empty elements
	// are skipped, mirroring trim_adapters).
	std::vector<std::string> adapters;
	const Value adapters_val = ExpressionExecutor::EvaluateScalar(ctx, *arguments[4]);
	if (!adapters_val.IsNull()) {
		for (const auto &child : ListValue::GetChildren(adapters_val)) {
			if (child.IsNull()) {
				continue;
			}
			auto s = child.GetValue<std::string>();
			if (!s.empty()) {
				adapters.push_back(std::move(s));
			}
		}
	}

	auto bd = make_uniq<TrimAdaptersPeBindData>();
	bd->overlap_require = overlap_require;
	bd->overlap_diff_limit = overlap_diff_limit;
	bd->overlap_diff_percent = overlap_diff_percent;
	bd->allow_pre_start = allow_pre_start;
	// Reuse trim_adapters' candidate builder (RC expansion, min_match auto-scale
	// from the pre-dedup count, then exact-duplicate dedup). Invalid DNA bases in
	// an adapter throw here at bind time.
	BuildAdapterCandidates(std::move(adapters), match_revcomp, min_match_param, bd->candidates, bd->min_match);
	return std::move(bd);
}

// trim_adapters_pe: infer the insert boundary from the R1 / revcomp(R2) overlap
// (fastp's primary PE adapter mechanism) and trim both mates to the insert when
// each reads through into adapter. When overlap analysis finds no adapter and a
// fallback adapter list was supplied (11-arg form), fall back to adapter-by-
// sequence on each mate independently (fastp step 8 -> step 9).
static void TrimAdaptersPeExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<TrimAdaptersPeLocalState>();
	const idx_t row_count = args.size();

	UnifiedVectorFormat seq1_data, qual1_data, seq2_data, qual2_data;
	args.data[0].ToUnifiedFormat(row_count, seq1_data);
	args.data[1].ToUnifiedFormat(row_count, qual1_data);
	args.data[2].ToUnifiedFormat(row_count, seq2_data);
	args.data[3].ToUnifiedFormat(row_count, qual2_data);

	auto seq1_ptr = UnifiedVectorFormat::GetData<string_t>(seq1_data);
	auto seq2_ptr = UnifiedVectorFormat::GetData<string_t>(seq2_data);

	auto &entries = StructVector::GetEntries(result);
	auto &seq1_out = *entries[0];
	auto &qual1_out = *entries[1];
	auto &seq2_out = *entries[2];
	auto &qual2_out = *entries[3];
	auto overlap_len_data = FlatVector::GetData<int32_t>(*entries[4]);
	auto adapter_trimmed_data = FlatVector::GetData<bool>(*entries[5]);
	auto trimmed1_data = FlatVector::GetData<uint32_t>(*entries[6]);
	auto trimmed2_data = FlatVector::GetData<uint32_t>(*entries[7]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual1_entries = FlatVector::GetData<list_entry_t>(qual1_out);
	auto qual2_entries = FlatVector::GetData<list_entry_t>(qual2_out);
	idx_t qual1_off = ListVector::GetListSize(qual1_out);
	idx_t qual2_off = ListVector::GetListSize(qual2_out);
	// Trimmed output is at most as long as the input quality, so one reserve per
	// mate before the row loop eliminates per-row reallocations.
	ListVector::Reserve(qual1_out, qual1_off + ListVector::GetListSize(args.data[1]));
	ListVector::Reserve(qual2_out, qual2_off + ListVector::GetListSize(args.data[3]));
	// Both Reserves happen before either child pointer is taken, and no Reserve
	// runs inside the row loop (WritePeMate only SetListSize, which never
	// reallocates), so these cached child pointers stay valid for the whole chunk.
	auto qual1_child = FlatVector::GetData<uint8_t>(ListVector::GetEntry(qual1_out));
	auto qual2_child = FlatVector::GetData<uint8_t>(ListVector::GetEntry(qual2_out));

	for (idx_t i = 0; i < row_count; i++) {
		auto s1i = seq1_data.sel->get_index(i);
		auto q1i = qual1_data.sel->get_index(i);
		auto s2i = seq2_data.sel->get_index(i);
		auto q2i = qual2_data.sel->get_index(i);

		if (!seq1_data.validity.RowIsValid(s1i) || !qual1_data.validity.RowIsValid(q1i) ||
		    !seq2_data.validity.RowIsValid(s2i) || !qual2_data.validity.RowIsValid(q2i)) {
			result_validity.SetInvalid(i);
			qual1_entries[i].offset = qual1_off;
			qual1_entries[i].length = 0;
			qual2_entries[i].offset = qual2_off;
			qual2_entries[i].length = 0;
			continue;
		}

		const auto &seq1 = seq1_ptr[s1i];
		const auto &seq2 = seq2_ptr[s2i];
		const uint8_t *q1ptr;
		const uint8_t *q2ptr;
		idx_t q1len, q2len;
		GetListUInt8Slice(args.data[1], qual1_data, i, q1ptr, q1len);
		GetListUInt8Slice(args.data[3], qual2_data, i, q2ptr, q2len);
		if (seq1.GetSize() != q1len) {
			throw InvalidInputException(
			    "trim_adapters_pe: read 1 sequence length (%llu) does not match quality length (%llu)",
			    (unsigned long long)seq1.GetSize(), (unsigned long long)q1len);
		}
		if (seq2.GetSize() != q2len) {
			throw InvalidInputException(
			    "trim_adapters_pe: read 2 sequence length (%llu) does not match quality length (%llu)",
			    (unsigned long long)seq2.GetSize(), (unsigned long long)q2len);
		}

		const auto ov = miint::qc::OverlapAnalyzer::analyze(seq1.GetData(), seq1.GetSize(), seq2.GetData(),
		                                                    seq2.GetSize(), lstate.overlap_diff_limit,
		                                                    lstate.overlap_require, lstate.overlap_diff_percent);

		idx_t keep1 = seq1.GetSize();
		idx_t keep2 = seq2.GetSize();
		bool adapter_trimmed = false;
		// offset < 0 means the insert is shorter than the read: each mate reads
		// through into adapter past the overlap, so trim both to the insert.
		// (overlap_len > 0 holds for any offset<0 result from OverlapAnalyzer; the
		// explicit check guards the unsigned cast against a degenerate result.)
		if (ov.overlapped && ov.offset < 0 && ov.overlap_len > 0) {
			const idx_t ol = static_cast<idx_t>(ov.overlap_len);
			keep1 = std::min(keep1, ol);
			keep2 = std::min(keep2, ol);
			adapter_trimmed = true;
		} else if (!lstate.candidates.empty()) {
			// fastp step 9: overlap analysis found no adapter, so fall back to
			// adapter-by-sequence on each mate independently (leftmost match wins).
			keep1 = miint::qc::AdapterMatcher::find_leftmost(reinterpret_cast<const uint8_t *>(seq1.GetData()),
			                                                 seq1.GetSize(), lstate.candidates, lstate.min_match,
			                                                 lstate.allow_pre_start);
			keep2 = miint::qc::AdapterMatcher::find_leftmost(reinterpret_cast<const uint8_t *>(seq2.GetData()),
			                                                 seq2.GetSize(), lstate.candidates, lstate.min_match,
			                                                 lstate.allow_pre_start);
			adapter_trimmed = keep1 < seq1.GetSize() || keep2 < seq2.GetSize();
		}

		overlap_len_data[i] = ov.overlapped ? ov.overlap_len : 0;
		adapter_trimmed_data[i] = adapter_trimmed;
		trimmed1_data[i] = static_cast<uint32_t>(seq1.GetSize() - keep1);
		trimmed2_data[i] = static_cast<uint32_t>(seq2.GetSize() - keep2);

		WritePeMate(i, seq1, q1ptr, keep1, seq1_out, qual1_out, qual1_entries, qual1_child, qual1_off);
		WritePeMate(i, seq2, q2ptr, keep2, seq2_out, qual2_out, qual2_entries, qual2_child, qual2_off);
	}
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

	// trim_adapters_pe: 4-arg overlap-only (fastp PE defaults) + 11-arg full form
	// (custom overlap params + adapter-by-sequence fallback). The bind captures
	// tuning params; init_local_state caches them per thread.
	{
		ScalarFunctionSet set("trim_adapters_pe");
		const auto qual_t = LogicalType::LIST(LogicalType::UTINYINT);
		const auto adapter_list_t = LogicalType::LIST(LogicalType::VARCHAR);
		const auto i = LogicalType::INTEGER;
		const auto b = LogicalType::BOOLEAN;

		ScalarFunction four_arg("trim_adapters_pe", {LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, qual_t},
		                        TrimAdaptersPeResultStructType(), TrimAdaptersPeExecute, TrimAdaptersPeBind4);
		four_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		four_arg.init_local_state = TrimAdaptersPeInitLocalState;
		set.AddFunction(four_arg);

		ScalarFunction eleven_arg(
		    "trim_adapters_pe",
		    {LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, qual_t, adapter_list_t, i, i, i, b, i, b},
		    TrimAdaptersPeResultStructType(), TrimAdaptersPeExecute, TrimAdaptersPeBind11);
		eleven_arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
		eleven_arg.init_local_state = TrimAdaptersPeInitLocalState;
		set.AddFunction(eleven_arg);

		loader.RegisterFunction(set);
	}
}

} // namespace duckdb
