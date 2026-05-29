#include "align_pairwise_ksw2_functions.hpp"

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
// Bind data: ksw2 (extz) penalty constants validated at bind time
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2BindData : public FunctionData {
	int match;
	int mismatch;
	int gap_open;
	int gap_extend;
	int w;
	int zdrop;

	AlignPairwiseKsw2BindData(int match_p, int mismatch_p, int gap_open_p, int gap_extend_p, int w_p, int zdrop_p)
	    : match(match_p), mismatch(mismatch_p), gap_open(gap_open_p), gap_extend(gap_extend_p), w(w_p), zdrop(zdrop_p) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<AlignPairwiseKsw2BindData>(match, mismatch, gap_open, gap_extend, w, zdrop);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<AlignPairwiseKsw2BindData>();
		return match == other.match && mismatch == other.mismatch && gap_open == other.gap_open &&
		       gap_extend == other.gap_extend && w == other.w && zdrop == other.zdrop;
	}

	// Penalty validation delegated to KSW2Aligner constructor (single source of truth)
	// via the shared exception-translation helper.
	static void ValidatePenalties(int match, int mismatch, int gap_open, int gap_extend, int w, int zdrop) {
		TranslateCtorInvalidArgument([&] { miint::KSW2Aligner test(match, mismatch, gap_open, gap_extend, w, zdrop); });
	}

	// 6-arg overload: (query, subject, match, mismatch, gap_open, gap_extend) — w=-1, zdrop=-1.
	static unique_ptr<AlignPairwiseKsw2BindData> FromArgs6(ClientContext &context,
	                                                       vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 6, "align_pairwise_ksw2_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, -1, -1);
		return make_uniq<AlignPairwiseKsw2BindData>(match, mismatch, gap_open, gap_extend, -1, -1);
	}

	// 8-arg overload: (query, subject, match, mismatch, gap_open, gap_extend, w, zdrop).
	static unique_ptr<AlignPairwiseKsw2BindData> FromArgs8(ClientContext &context,
	                                                       vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 8, "align_pairwise_ksw2_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		auto w = EvalIntArg(context, arguments, 6);
		auto zdrop = EvalIntArg(context, arguments, 7);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, w, zdrop);
		return make_uniq<AlignPairwiseKsw2BindData>(match, mismatch, gap_open, gap_extend, w, zdrop);
	}

	static unique_ptr<AlignPairwiseKsw2BindData> Defaults() {
		return make_uniq<AlignPairwiseKsw2BindData>(2, 4, 6, 2, -1, -1);
	}
};

// ---------------------------------------------------------------------------
// Local state: per-thread KSW2Aligner reused across all rows
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2LocalState : public FunctionLocalState {
	miint::KSW2Aligner aligner;
	std::string query_buf;
	std::string subject_buf;

	AlignPairwiseKsw2LocalState(int match, int mismatch, int gap_open, int gap_extend, int w, int zdrop)
	    : aligner(match, mismatch, gap_open, gap_extend, w, zdrop) {
		query_buf.reserve(256);
		subject_buf.reserve(256);
	}
};

static unique_ptr<FunctionLocalState>
AlignPairwiseKsw2InitLocalState(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data) {
	auto &data = bind_data->Cast<AlignPairwiseKsw2BindData>();
	return make_uniq<AlignPairwiseKsw2LocalState>(data.match, data.mismatch, data.gap_open, data.gap_extend, data.w,
	                                              data.zdrop);
}

// ---------------------------------------------------------------------------
// Argument-type vectors
// ---------------------------------------------------------------------------
// (query, subject, match, mismatch, gap_open, gap_extend)
static vector<LogicalType> SixArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER};
}

// (query, subject, match, mismatch, gap_open, gap_extend, w, zdrop)
static vector<LogicalType> EightArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER};
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_score → INTEGER
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2ScoreExecute =
    RunPairwiseAlignScoreExecute<AlignPairwiseKsw2LocalState, &miint::KSW2Aligner::align_extz_score>;

void AlignPairwiseKsw2ScoreFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_score");

	ScalarFunction score_2arg("align_pairwise_ksw2_score", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          LogicalType::INTEGER, AlignPairwiseKsw2ScoreExecute);
	score_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::Defaults().release());
	};
	score_2arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(score_2arg);

	ScalarFunction score_6arg("align_pairwise_ksw2_score", SixArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseKsw2ScoreExecute);
	score_6arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_6arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs6(ctx, args).release());
	};
	score_6arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(score_6arg);

	ScalarFunction score_8arg("align_pairwise_ksw2_score", EightArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseKsw2ScoreExecute);
	score_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs8(ctx, args).release());
	};
	score_8arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(score_8arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_cigar → STRUCT(score INTEGER, cigar VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2CigarExecute =
    RunPairwiseAlignCigarExecute<AlignPairwiseKsw2LocalState, &miint::KSW2Aligner::align_extz_cigar>;

void AlignPairwiseKsw2CigarFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_cigar");

	ScalarFunction cigar_2arg("align_pairwise_ksw2_cigar", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          PairwiseCigarReturnType(), AlignPairwiseKsw2CigarExecute);
	cigar_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::Defaults().release());
	};
	cigar_2arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(cigar_2arg);

	ScalarFunction cigar_6arg("align_pairwise_ksw2_cigar", SixArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseKsw2CigarExecute);
	cigar_6arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_6arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs6(ctx, args).release());
	};
	cigar_6arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(cigar_6arg);

	ScalarFunction cigar_8arg("align_pairwise_ksw2_cigar", EightArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseKsw2CigarExecute);
	cigar_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs8(ctx, args).release());
	};
	cigar_8arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(cigar_8arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_full → STRUCT(score INTEGER, cigar VARCHAR,
//                                    query_aligned VARCHAR, subject_aligned VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2FullExecute =
    RunPairwiseAlignFullExecute<AlignPairwiseKsw2LocalState, &miint::KSW2Aligner::align_extz_full>;

void AlignPairwiseKsw2FullFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_full");

	ScalarFunction full_2arg("align_pairwise_ksw2_full", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         PairwiseFullReturnType(), AlignPairwiseKsw2FullExecute);
	full_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::Defaults().release());
	};
	full_2arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(full_2arg);

	ScalarFunction full_6arg("align_pairwise_ksw2_full", SixArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseKsw2FullExecute);
	full_6arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_6arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs6(ctx, args).release());
	};
	full_6arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(full_6arg);

	ScalarFunction full_8arg("align_pairwise_ksw2_full", EightArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseKsw2FullExecute);
	full_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2BindData::FromArgs8(ctx, args).release());
	};
	full_8arg.init_local_state = AlignPairwiseKsw2InitLocalState;
	function_set.AddFunction(full_8arg);

	loader.RegisterFunction(function_set);
}

} // namespace duckdb
