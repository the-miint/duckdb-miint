#include "extract_linked_amplicon.hpp"

#include "WFA2Aligner.hpp"
#include "alignment_functions_internal.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"

#include <cmath>
#include <cstring>
#include <limits>
#include <optional>
#include <string>

namespace duckdb {

// ---------------------------------------------------------------------------
// Defaults
// ---------------------------------------------------------------------------
static constexpr double DEFAULT_ERROR_RATE = 0.10;
static constexpr int64_t DEFAULT_MIN_LEN = 0;
static constexpr int64_t DEFAULT_MAX_LEN = std::numeric_limits<int32_t>::max();
// Matches cutadapt's `-O` default. Floors how many anchor bases must align
// in partial (non-anchored) mode; prevents spurious short partial overlaps.
static constexpr int64_t DEFAULT_MIN_OVERLAP = 3;

// ---------------------------------------------------------------------------
// Output struct shape
// ---------------------------------------------------------------------------
static LogicalType LinkedAmpliconResultType() {
	return LogicalType::STRUCT({{"sequence", LogicalType::VARCHAR},
	                            {"qual", LogicalType::LIST(LogicalType::UTINYINT)},
	                            {"start", LogicalType::INTEGER},
	                            {"stop", LogicalType::INTEGER}});
}

// ---------------------------------------------------------------------------
// Bind data — three constants captured at bind time
// ---------------------------------------------------------------------------
struct LinkedAmpliconBindData : public FunctionData {
	int64_t min_len;
	int64_t max_len;
	double error_rate;
	int64_t min_overlap;

	LinkedAmpliconBindData(int64_t mn, int64_t mx, double er, int64_t mo)
	    : min_len(mn), max_len(mx), error_rate(er), min_overlap(mo) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<LinkedAmpliconBindData>(min_len, max_len, error_rate, min_overlap);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<LinkedAmpliconBindData>();
		return min_len == other.min_len && max_len == other.max_len && error_rate == other.error_rate &&
		       min_overlap == other.min_overlap;
	}

	static void Validate(int64_t min_len, int64_t max_len, double error_rate, int64_t min_overlap) {
		if (min_len < 0) {
			throw InvalidInputException("extract_linked_amplicon: min_len must be >= 0 (got %lld)",
			                            static_cast<long long>(min_len));
		}
		if (max_len < 0) {
			throw InvalidInputException("extract_linked_amplicon: max_len must be >= 0 (got %lld)",
			                            static_cast<long long>(max_len));
		}
		if (min_len > max_len) {
			throw InvalidInputException("extract_linked_amplicon: min_len (%lld) must be <= max_len (%lld)",
			                            static_cast<long long>(min_len), static_cast<long long>(max_len));
		}
		if (!(error_rate >= 0.0 && error_rate <= 1.0)) {
			throw InvalidInputException("extract_linked_amplicon: error_rate must be in [0.0, 1.0] (got %f)",
			                            error_rate);
		}
		if (min_overlap < 1) {
			throw InvalidInputException("extract_linked_amplicon: min_overlap must be >= 1 (got %lld)",
			                            static_cast<long long>(min_overlap));
		}
	}
};

// ---------------------------------------------------------------------------
// Per-thread local state — reuses one aligner and string buffers across rows
// ---------------------------------------------------------------------------
struct LinkedAmpliconLocalState : public FunctionLocalState {
	miint::WFA2Aligner aligner;
	std::string anchor_buf;
	std::string window_buf;
	std::string suffix_buf;
	// Bind-time params cached here so Execute doesn't need to dig into bind_info
	int64_t min_len;
	int64_t max_len;
	double error_rate;
	int64_t min_overlap;

	LinkedAmpliconLocalState(int64_t mn, int64_t mx, double er, int64_t mo)
	    : aligner(), min_len(mn), max_len(mx), error_rate(er), min_overlap(mo) {
		anchor_buf.reserve(64);
		window_buf.reserve(1024);
		suffix_buf.reserve(1024);
	}
};

static unique_ptr<FunctionLocalState> InitLocalState(ExpressionState &state, const BoundFunctionExpression &expr,
                                                     FunctionData *bind_data) {
	(void)state;
	(void)expr;
	auto &bd = bind_data->Cast<LinkedAmpliconBindData>();
	return make_uniq<LinkedAmpliconLocalState>(bd.min_len, bd.max_len, bd.error_rate, bd.min_overlap);
}

// ---------------------------------------------------------------------------
// Anchor matching against a window via WFA2 semi-global.
//
// WFA2 semi-global convention (see WFA2Aligner::align_full_semiglobal):
//   pattern = subject (window) with free begin/end gaps
//   text    = query   (anchor); text_begin_free / text_end_free allow free
//   trim of the anchor's prefix/suffix (partial / non-anchored mode)
//
// CIGAR interpretation for our purposes:
//   - leading 'D' runs       = window bases skipped at start (free begin)
//   - trailing 'D' runs      = window bases skipped at end (free end)
//   - leading 'I' runs       = anchor bases trimmed at start, only free when
//                              text_begin_free > 0 (cutadapt's non-anchored
//                              5' adapter partial-overlap-at-read-start)
//   - trailing 'I' runs      = anchor bases trimmed at end, only free when
//                              text_end_free > 0 (mirror case at read 3')
//   - inner ops (=, X, I, D) form the anchor↔window alignment
//   - errors = sum of inner X + I + D lengths
//
// Returns std::nullopt if the alignment fails OR inner errors > budget.
// ---------------------------------------------------------------------------
struct AnchorMatch {
	size_t start;
	size_t stop;
	int errors;
};

static std::optional<AnchorMatch> FindAnchor(miint::WFA2Aligner &aligner, const std::string &anchor,
                                             const std::string &window, int error_budget, int text_begin_free,
                                             int text_end_free) {
	// An empty window can never carry useful match information. Cap checks below
	// would also reject an all-I CIGAR (since trail_i == anchor.size() > any
	// in-budget text_end_free), but guarding here is cheaper and keeps the
	// safety invariant local rather than emergent from the strip+cap interplay.
	if (window.empty()) {
		return std::nullopt;
	}
	// Short-circuit: a fully-anchored anchor longer than the window cannot fit.
	// In partial mode the anchor *can* legitimately exceed the window length —
	// e.g., long primer with most of its 3' tail off the read end — so only
	// short-circuit when both ends are anchored.
	if (anchor.size() > window.size() && text_begin_free == 0 && text_end_free == 0) {
		return std::nullopt;
	}
	auto aln =
	    aligner.align_cigar_semiglobal_iupac(/*query=*/anchor, /*subject=*/window, text_begin_free, text_end_free);
	if (!aln.has_value()) {
		return std::nullopt;
	}

	std::vector<miint::CigarOperation> ops;
	try {
		ops = miint::ParseCigarOperations(aln->cigar);
	} catch (const miint::InvalidInputException &) {
		return std::nullopt;
	}
	if (ops.empty()) {
		return std::nullopt;
	}

	// Strip leading free ops: D always (window-skip), I when partial-5'.
	// Either order is valid; consume both kinds until we hit a non-free op.
	size_t lead_d = 0;
	size_t lead_i = 0;
	size_t i = 0;
	while (i < ops.size()) {
		if (ops[i].op == 'D') {
			lead_d += static_cast<size_t>(ops[i].length);
			++i;
		} else if (ops[i].op == 'I' && text_begin_free > 0) {
			lead_i += static_cast<size_t>(ops[i].length);
			++i;
		} else {
			break;
		}
	}

	// Strip trailing free ops: D always, I when partial-3'.
	size_t trail_d = 0;
	size_t trail_i = 0;
	size_t j = ops.size();
	while (j > i) {
		if (ops[j - 1].op == 'D') {
			trail_d += static_cast<size_t>(ops[j - 1].length);
			--j;
		} else if (ops[j - 1].op == 'I' && text_end_free > 0) {
			trail_i += static_cast<size_t>(ops[j - 1].length);
			--j;
		} else {
			break;
		}
	}

	// Enforce the per-side free-trim caps. WFA2 will use as much free trim as
	// it likes; we cap at the configured budget so min_overlap is respected
	// even if a sub-min_overlap match would have been cheaper for WFA2.
	if (static_cast<int>(lead_i) > text_begin_free) {
		return std::nullopt;
	}
	if (static_cast<int>(trail_i) > text_end_free) {
		return std::nullopt;
	}

	// Count errors in inner [i, j)
	int errors = 0;
	for (size_t k = i; k < j; ++k) {
		switch (ops[k].op) {
		case 'X': // mismatch
		case 'I': // insertion in anchor (anchor has extra base)
		case 'D': // deletion in anchor (window has extra base mid-alignment)
			errors += static_cast<int>(ops[k].length);
			break;
		default:
			break;
		}
	}
	if (errors > error_budget) {
		return std::nullopt;
	}

	AnchorMatch m;
	m.start = lead_d;
	m.stop = window.size() - trail_d;
	m.errors = errors;
	return m;
}

// ---------------------------------------------------------------------------
// Execute
// ---------------------------------------------------------------------------
static void Execute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<LinkedAmpliconLocalState>();

	const idx_t row_count = args.size();

	UnifiedVectorFormat seq_data, qual_data, a5_data, a3_data;
	args.data[0].ToUnifiedFormat(row_count, seq_data);
	args.data[1].ToUnifiedFormat(row_count, qual_data);
	args.data[2].ToUnifiedFormat(row_count, a5_data);
	args.data[3].ToUnifiedFormat(row_count, a3_data);

	auto seq_ptr = UnifiedVectorFormat::GetData<string_t>(seq_data);
	auto a5_ptr = UnifiedVectorFormat::GetData<string_t>(a5_data);
	auto a3_ptr = UnifiedVectorFormat::GetData<string_t>(a3_data);

	auto &entries = StructVector::GetEntries(result);
	auto &seq_out_vec = *entries[0];
	auto &qual_out_vec = *entries[1];
	auto start_data = FlatVector::GetData<int32_t>(*entries[2]);
	auto stop_data = FlatVector::GetData<int32_t>(*entries[3]);
	auto &result_validity = FlatVector::Validity(result);

	auto qual_out_entries = FlatVector::GetData<list_entry_t>(qual_out_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_out_vec);
	// Worst case: every output row keeps the full input qual list.
	ListVector::Reserve(qual_out_vec, qual_child_offset + ListVector::GetListSize(args.data[1]));

	auto invalidate = [&](idx_t i) {
		result_validity.SetInvalid(i);
		FlatVector::Validity(seq_out_vec).SetInvalid(i);
		FlatVector::Validity(qual_out_vec).SetInvalid(i);
		FlatVector::Validity(*entries[2]).SetInvalid(i);
		FlatVector::Validity(*entries[3]).SetInvalid(i);
		qual_out_entries[i].offset = qual_child_offset;
		qual_out_entries[i].length = 0;
		start_data[i] = 0;
		stop_data[i] = 0;
	};

	for (idx_t i = 0; i < row_count; ++i) {
		auto si = seq_data.sel->get_index(i);
		auto qi = qual_data.sel->get_index(i);
		auto a5i = a5_data.sel->get_index(i);
		auto a3i = a3_data.sel->get_index(i);

		bool all_valid = seq_data.validity.RowIsValid(si) && qual_data.validity.RowIsValid(qi) &&
		                 a5_data.validity.RowIsValid(a5i) && a3_data.validity.RowIsValid(a3i);
		if (!all_valid) {
			invalidate(i);
			continue;
		}

		const auto &seq = seq_ptr[si];
		const auto &a5 = a5_ptr[a5i];
		const auto &a3 = a3_ptr[a3i];

		// Strip cutadapt-style sigils: '^' prefix on anchor5 anchors its 5' edge
		// to the read start; '$' suffix on anchor3 anchors its 3' edge to the
		// read end. Default (no sigil) is cutadapt's non-anchored mode: the
		// anchor's far end may hang off the read edge subject to min_overlap.
		// Local raw/_size names avoid shadowing the function-scope
		// UnifiedVectorFormat a5_data / a3_data above.
		const char *a5_raw = a5.GetData();
		idx_t a5_raw_size = a5.GetSize();
		bool anchor5_anchored = false;
		if (a5_raw_size > 0 && a5_raw[0] == '^') {
			anchor5_anchored = true;
			++a5_raw;
			--a5_raw_size;
		}
		const char *a3_raw = a3.GetData();
		idx_t a3_raw_size = a3.GetSize();
		bool anchor3_anchored = false;
		if (a3_raw_size > 0 && a3_raw[a3_raw_size - 1] == '$') {
			anchor3_anchored = true;
			--a3_raw_size;
		}

		if (a5_raw_size == 0 || a3_raw_size == 0) {
			throw InvalidInputException("extract_linked_amplicon: anchor sequences must be non-empty");
		}

		const uint8_t *qptr;
		idx_t qlen;
		GetListUInt8Slice(args.data[1], qual_data, i, qptr, qlen);
		if (seq.GetSize() != qlen) {
			throw InvalidInputException(
			    "extract_linked_amplicon: sequence length (%llu) does not match quality length (%llu)",
			    (unsigned long long)seq.GetSize(), (unsigned long long)qlen);
		}

		lstate.window_buf.assign(seq.GetData(), seq.GetSize());
		lstate.anchor_buf.assign(a5_raw, a5_raw_size);

		// Partial-5' free-trim cap: in non-anchored mode the anchor's 5' end
		// may hang off the window's 5' edge by up to len(anchor) - min_overlap.
		// If min_overlap >= len(anchor) the cap is 0, i.e., partial mode
		// silently degenerates to anchored for too-short anchors rather than
		// throwing (anchors may be column-driven).
		int tbf5_cap = 0;
		if (!anchor5_anchored) {
			int64_t cap = static_cast<int64_t>(a5_raw_size) - lstate.min_overlap;
			tbf5_cap = cap > 0 ? static_cast<int>(cap) : 0;
		}

		// Find anchor5: two-pass to match cutadapt's preference for full
		// anchor matches over partial-overlap-at-edge. Pass 1 calls WFA2 with
		// text_begin_free = 0, which *structurally disables* the partial path
		// (not merely deprioritizes it) — WFA2 has no representation for "trim
		// the anchor prefix for free" and must align the anchor end-to-end.
		// Pass 2 runs only when pass 1 found nothing, lifting the cap so the
		// edge-overlap case can match. This avoids WFA2's internal tie-break
		// (which favors maximal free-trim) producing a partial-at-edge match
		// when a clean internal full match also exists at zero cost.
		int budget5 = static_cast<int>(std::ceil(static_cast<double>(lstate.anchor_buf.size()) * lstate.error_rate));
		auto m5 = FindAnchor(lstate.aligner, lstate.anchor_buf, lstate.window_buf, budget5, /*tbf=*/0, /*tef=*/0);
		if (!m5.has_value() && tbf5_cap > 0) {
			m5 = FindAnchor(lstate.aligner, lstate.anchor_buf, lstate.window_buf, budget5, tbf5_cap, /*tef=*/0);
		}
		if (!m5.has_value()) {
			invalidate(i);
			continue;
		}

		// Find anchor3 in seq[m5.stop:]
		lstate.suffix_buf.assign(seq.GetData() + m5->stop, seq.GetSize() - m5->stop);
		lstate.anchor_buf.assign(a3_raw, a3_raw_size);
		int budget3 = static_cast<int>(std::ceil(static_cast<double>(lstate.anchor_buf.size()) * lstate.error_rate));

		// Partial-3' free-trim cap (symmetric: anchor's 3' end may hang off
		// the window's 3' edge by up to len(anchor) - min_overlap).
		int tef3_cap = 0;
		if (!anchor3_anchored) {
			int64_t cap = static_cast<int64_t>(a3_raw_size) - lstate.min_overlap;
			tef3_cap = cap > 0 ? static_cast<int>(cap) : 0;
		}

		// Same two-pass policy for anchor3 (anchored first structurally
		// disables the partial path; partial fallback runs only on a miss).
		auto m3 = FindAnchor(lstate.aligner, lstate.anchor_buf, lstate.suffix_buf, budget3, /*tbf=*/0, /*tef=*/0);
		if (!m3.has_value() && tef3_cap > 0) {
			m3 = FindAnchor(lstate.aligner, lstate.anchor_buf, lstate.suffix_buf, budget3, /*tbf=*/0, tef3_cap);
		}
		if (!m3.has_value()) {
			invalidate(i);
			continue;
		}

		// Translate anchor3 match back to global window coordinates
		size_t global_start = m5->stop;
		size_t global_stop = m5->stop + m3->start;
		size_t extracted_len = global_stop - global_start;

		// Length bounds
		if (static_cast<int64_t>(extracted_len) < lstate.min_len ||
		    static_cast<int64_t>(extracted_len) > lstate.max_len) {
			invalidate(i);
			continue;
		}

		// Write extracted sequence
		FlatVector::GetData<string_t>(seq_out_vec)[i] =
		    StringVector::AddString(seq_out_vec, seq.GetData() + global_start, extracted_len);

		// Write extracted qual list
		auto &qual_child = ListVector::GetEntry(qual_out_vec);
		auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
		if (extracted_len > 0) {
			std::memcpy(qual_child_data + qual_child_offset, qptr + global_start, extracted_len);
		}
		qual_out_entries[i].offset = qual_child_offset;
		qual_out_entries[i].length = extracted_len;
		qual_child_offset += extracted_len;
		ListVector::SetListSize(qual_out_vec, qual_child_offset);

		start_data[i] = static_cast<int32_t>(global_start);
		stop_data[i] = static_cast<int32_t>(global_stop);
	}
}

// ---------------------------------------------------------------------------
// Bind: 4-arg (all defaults), 7-arg (min_len/max_len/error_rate), 8-arg
// (also min_overlap). All explicit params must be foldable constants.
// ---------------------------------------------------------------------------
static unique_ptr<FunctionData> Bind4Arg(ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
	(void)ctx;
	(void)fn;
	(void)args;
	LinkedAmpliconBindData::Validate(DEFAULT_MIN_LEN, DEFAULT_MAX_LEN, DEFAULT_ERROR_RATE, DEFAULT_MIN_OVERLAP);
	return make_uniq<LinkedAmpliconBindData>(DEFAULT_MIN_LEN, DEFAULT_MAX_LEN, DEFAULT_ERROR_RATE,
	                                         DEFAULT_MIN_OVERLAP);
}

static unique_ptr<FunctionData> Bind7Arg(ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
	(void)fn;
	for (idx_t i = 4; i < 7; ++i) {
		if (!args[i]->IsFoldable()) {
			throw InvalidInputException(
			    "extract_linked_amplicon: min_len, max_len, error_rate must be constants, not column references");
		}
	}
	int64_t min_len = ExpressionExecutor::EvaluateScalar(ctx, *args[4]).GetValue<int64_t>();
	int64_t max_len = ExpressionExecutor::EvaluateScalar(ctx, *args[5]).GetValue<int64_t>();
	double error_rate = ExpressionExecutor::EvaluateScalar(ctx, *args[6]).GetValue<double>();
	LinkedAmpliconBindData::Validate(min_len, max_len, error_rate, DEFAULT_MIN_OVERLAP);
	return make_uniq<LinkedAmpliconBindData>(min_len, max_len, error_rate, DEFAULT_MIN_OVERLAP);
}

static unique_ptr<FunctionData> Bind8Arg(ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
	(void)fn;
	for (idx_t i = 4; i < 8; ++i) {
		if (!args[i]->IsFoldable()) {
			throw InvalidInputException("extract_linked_amplicon: min_len, max_len, error_rate, min_overlap must be "
			                            "constants, not column references");
		}
	}
	int64_t min_len = ExpressionExecutor::EvaluateScalar(ctx, *args[4]).GetValue<int64_t>();
	int64_t max_len = ExpressionExecutor::EvaluateScalar(ctx, *args[5]).GetValue<int64_t>();
	double error_rate = ExpressionExecutor::EvaluateScalar(ctx, *args[6]).GetValue<double>();
	int64_t min_overlap = ExpressionExecutor::EvaluateScalar(ctx, *args[7]).GetValue<int64_t>();
	LinkedAmpliconBindData::Validate(min_len, max_len, error_rate, min_overlap);
	return make_uniq<LinkedAmpliconBindData>(min_len, max_len, error_rate, min_overlap);
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------
void ExtractLinkedAmpliconFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet set("extract_linked_amplicon");

	auto qual_t = LogicalType::LIST(LogicalType::UTINYINT);
	auto ret_t = LinkedAmpliconResultType();

	// 4-arg: seq, qual, anchor5, anchor3
	ScalarFunction four({LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, LogicalType::VARCHAR}, ret_t, Execute,
	                    Bind4Arg);
	four.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	four.init_local_state = InitLocalState;
	set.AddFunction(four);

	// 7-arg: seq, qual, anchor5, anchor3, min_len, max_len, error_rate
	ScalarFunction seven({LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BIGINT,
	                      LogicalType::BIGINT, LogicalType::DOUBLE},
	                     ret_t, Execute, Bind7Arg);
	seven.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	seven.init_local_state = InitLocalState;
	set.AddFunction(seven);

	// 8-arg: seq, qual, anchor5, anchor3, min_len, max_len, error_rate, min_overlap
	ScalarFunction eight({LogicalType::VARCHAR, qual_t, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BIGINT,
	                      LogicalType::BIGINT, LogicalType::DOUBLE, LogicalType::BIGINT},
	                     ret_t, Execute, Bind8Arg);
	eight.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	eight.init_local_state = InitLocalState;
	set.AddFunction(eight);

	loader.RegisterFunction(set);
}

} // namespace duckdb
