#include "align_pairwise_ksw2_dual_affine_functions.hpp"

#include "KSW2Aligner.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "pairwise_align_shared.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// Bind data: dual-affine KSW2 penalty constants validated at bind time
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2DualAffineBindData : public FunctionData {
	int match;
	int mismatch;
	int gap_open;
	int gap_extend;
	int gap_open2;
	int gap_extend2;
	int w;
	int zdrop;

	AlignPairwiseKsw2DualAffineBindData(int match_p, int mismatch_p, int gap_open_p, int gap_extend_p, int gap_open2_p,
	                                    int gap_extend2_p, int w_p, int zdrop_p)
	    : match(match_p), mismatch(mismatch_p), gap_open(gap_open_p), gap_extend(gap_extend_p), gap_open2(gap_open2_p),
	      gap_extend2(gap_extend2_p), w(w_p), zdrop(zdrop_p) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<AlignPairwiseKsw2DualAffineBindData>(match, mismatch, gap_open, gap_extend, gap_open2,
		                                                      gap_extend2, w, zdrop);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<AlignPairwiseKsw2DualAffineBindData>();
		return match == other.match && mismatch == other.mismatch && gap_open == other.gap_open &&
		       gap_extend == other.gap_extend && gap_open2 == other.gap_open2 && gap_extend2 == other.gap_extend2 &&
		       w == other.w && zdrop == other.zdrop;
	}

	// Penalty validation delegated to KSW2Aligner dual-affine ctor (single source of truth)
	// via the shared exception-translation helper.
	static void ValidatePenalties(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int gap_extend2,
	                              int w, int zdrop) {
		TranslateCtorInvalidArgument(
		    [&] { miint::KSW2Aligner test(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, w, zdrop); });
	}

	// 8-arg: (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2)
	// w and zdrop default to -1 (no band, no z-drop).
	static unique_ptr<AlignPairwiseKsw2DualAffineBindData> FromArgs8(ClientContext &context,
	                                                                 vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 8, "align_pairwise_ksw2_dual_affine_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		auto gap_open2 = EvalIntArg(context, arguments, 6);
		auto gap_extend2 = EvalIntArg(context, arguments, 7);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, -1, -1);
		return make_uniq<AlignPairwiseKsw2DualAffineBindData>(match, mismatch, gap_open, gap_extend, gap_open2,
		                                                      gap_extend2, -1, -1);
	}

	// 10-arg: (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, w, zdrop)
	static unique_ptr<AlignPairwiseKsw2DualAffineBindData> FromArgs10(ClientContext &context,
	                                                                  vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 10, "align_pairwise_ksw2_dual_affine_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		auto gap_open2 = EvalIntArg(context, arguments, 6);
		auto gap_extend2 = EvalIntArg(context, arguments, 7);
		auto w = EvalIntArg(context, arguments, 8);
		auto zdrop = EvalIntArg(context, arguments, 9);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, w, zdrop);
		return make_uniq<AlignPairwiseKsw2DualAffineBindData>(match, mismatch, gap_open, gap_extend, gap_open2,
		                                                      gap_extend2, w, zdrop);
	}

	// Defaults: match=2, mismatch=4, gap_open=6, gap_extend=2, gap_open2=24, gap_extend2=1.
	// Minimap2-style long-read defaults: first affine cheap-open for SNVs/micro-indels;
	// second affine expensive-open + cheap-extend for long indels.
	static unique_ptr<AlignPairwiseKsw2DualAffineBindData> Defaults() {
		return make_uniq<AlignPairwiseKsw2DualAffineBindData>(2, 4, 6, 2, 24, 1, -1, -1);
	}
};

// ---------------------------------------------------------------------------
// Local state: per-thread KSW2Aligner constructed with dual-affine params
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2DualAffineLocalState : public FunctionLocalState {
	miint::KSW2Aligner aligner;
	std::string query_buf;
	std::string subject_buf;

	AlignPairwiseKsw2DualAffineLocalState(int match, int mismatch, int gap_open, int gap_extend, int gap_open2,
	                                      int gap_extend2, int w, int zdrop)
	    : aligner(match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, w, zdrop) {
		query_buf.reserve(256);
		subject_buf.reserve(256);
	}
};

static unique_ptr<FunctionLocalState> AlignPairwiseKsw2DualAffineInitLocalState(ExpressionState &state,
                                                                                const BoundFunctionExpression &expr,
                                                                                FunctionData *bind_data) {
	auto &data = bind_data->Cast<AlignPairwiseKsw2DualAffineBindData>();
	return make_uniq<AlignPairwiseKsw2DualAffineLocalState>(data.match, data.mismatch, data.gap_open, data.gap_extend,
	                                                        data.gap_open2, data.gap_extend2, data.w, data.zdrop);
}

// ---------------------------------------------------------------------------
// Argument-type vectors
// ---------------------------------------------------------------------------
// (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2)
static vector<LogicalType> EightArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER};
}

// (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, gap_extend2, w, zdrop)
static vector<LogicalType> TenArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER};
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_dual_affine_score → INTEGER
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2DualAffineScoreExecute =
    RunPairwiseAlignScoreExecute<AlignPairwiseKsw2DualAffineLocalState, &miint::KSW2Aligner::align_extd_score>;

void AlignPairwiseKsw2DualAffineScoreFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_dual_affine_score");

	ScalarFunction score_2arg("align_pairwise_ksw2_dual_affine_score", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          LogicalType::INTEGER, AlignPairwiseKsw2DualAffineScoreExecute);
	score_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::Defaults().release());
	};
	score_2arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(score_2arg);

	ScalarFunction score_8arg("align_pairwise_ksw2_dual_affine_score", EightArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseKsw2DualAffineScoreExecute);
	score_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs8(ctx, args).release());
	};
	score_8arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(score_8arg);

	ScalarFunction score_10arg("align_pairwise_ksw2_dual_affine_score", TenArgTypes(), LogicalType::INTEGER,
	                           AlignPairwiseKsw2DualAffineScoreExecute);
	score_10arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_10arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs10(ctx, args).release());
	};
	score_10arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(score_10arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_dual_affine_cigar → STRUCT(score INTEGER, cigar VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2DualAffineCigarExecute =
    RunPairwiseAlignCigarExecute<AlignPairwiseKsw2DualAffineLocalState, &miint::KSW2Aligner::align_extd_cigar>;

void AlignPairwiseKsw2DualAffineCigarFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_dual_affine_cigar");

	ScalarFunction cigar_2arg("align_pairwise_ksw2_dual_affine_cigar", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          PairwiseCigarReturnType(), AlignPairwiseKsw2DualAffineCigarExecute);
	cigar_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::Defaults().release());
	};
	cigar_2arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(cigar_2arg);

	ScalarFunction cigar_8arg("align_pairwise_ksw2_dual_affine_cigar", EightArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseKsw2DualAffineCigarExecute);
	cigar_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs8(ctx, args).release());
	};
	cigar_8arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(cigar_8arg);

	ScalarFunction cigar_10arg("align_pairwise_ksw2_dual_affine_cigar", TenArgTypes(), PairwiseCigarReturnType(),
	                           AlignPairwiseKsw2DualAffineCigarExecute);
	cigar_10arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_10arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs10(ctx, args).release());
	};
	cigar_10arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(cigar_10arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_dual_affine_full → STRUCT(score INTEGER, cigar VARCHAR,
//                                                query_aligned VARCHAR, subject_aligned VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2DualAffineFullExecute =
    RunPairwiseAlignFullExecute<AlignPairwiseKsw2DualAffineLocalState, &miint::KSW2Aligner::align_extd_full>;

void AlignPairwiseKsw2DualAffineFullFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_dual_affine_full");

	ScalarFunction full_2arg("align_pairwise_ksw2_dual_affine_full", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         PairwiseFullReturnType(), AlignPairwiseKsw2DualAffineFullExecute);
	full_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::Defaults().release());
	};
	full_2arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(full_2arg);

	ScalarFunction full_8arg("align_pairwise_ksw2_dual_affine_full", EightArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseKsw2DualAffineFullExecute);
	full_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs8(ctx, args).release());
	};
	full_8arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(full_8arg);

	ScalarFunction full_10arg("align_pairwise_ksw2_dual_affine_full", TenArgTypes(), PairwiseFullReturnType(),
	                          AlignPairwiseKsw2DualAffineFullExecute);
	full_10arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_10arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2DualAffineBindData::FromArgs10(ctx, args).release());
	};
	full_10arg.init_local_state = AlignPairwiseKsw2DualAffineInitLocalState;
	function_set.AddFunction(full_10arg);

	loader.RegisterFunction(function_set);
}

} // namespace duckdb
