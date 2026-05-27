#include "align_pairwise_ksw2_splice_functions.hpp"

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
// Bind data: splice-mode penalty constants validated at bind time
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2SpliceBindData : public FunctionData {
	int match;
	int mismatch;
	int gap_open;
	int gap_extend;
	int gap_open2;
	int noncan;
	int zdrop;

	AlignPairwiseKsw2SpliceBindData(int match_p, int mismatch_p, int gap_open_p, int gap_extend_p, int gap_open2_p,
	                                int noncan_p, int zdrop_p)
	    : match(match_p), mismatch(mismatch_p), gap_open(gap_open_p), gap_extend(gap_extend_p), gap_open2(gap_open2_p),
	      noncan(noncan_p), zdrop(zdrop_p) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<AlignPairwiseKsw2SpliceBindData>(match, mismatch, gap_open, gap_extend, gap_open2, noncan,
		                                                  zdrop);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<AlignPairwiseKsw2SpliceBindData>();
		return match == other.match && mismatch == other.mismatch && gap_open == other.gap_open &&
		       gap_extend == other.gap_extend && gap_open2 == other.gap_open2 && noncan == other.noncan &&
		       zdrop == other.zdrop;
	}

	// Penalty validation delegated to KSW2Aligner splice ctor via the shared exception helper.
	static void ValidatePenalties(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int noncan,
	                              int zdrop) {
		TranslateCtorInvalidArgument(
		    [&] { miint::KSW2Aligner test(match, mismatch, gap_open, gap_extend, gap_open2, noncan, zdrop); });
	}

	// 8-arg: (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, noncan)
	// zdrop defaults to -1.
	static unique_ptr<AlignPairwiseKsw2SpliceBindData> FromArgs8(ClientContext &context,
	                                                             vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 8, "align_pairwise_ksw2_splice_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		auto gap_open2 = EvalIntArg(context, arguments, 6);
		auto noncan = EvalIntArg(context, arguments, 7);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, gap_open2, noncan, -1);
		return make_uniq<AlignPairwiseKsw2SpliceBindData>(match, mismatch, gap_open, gap_extend, gap_open2, noncan, -1);
	}

	// 9-arg: (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, noncan, zdrop)
	static unique_ptr<AlignPairwiseKsw2SpliceBindData> FromArgs9(ClientContext &context,
	                                                             vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 9, "align_pairwise_ksw2_splice_*");
		auto match = EvalIntArg(context, arguments, 2);
		auto mismatch = EvalIntArg(context, arguments, 3);
		auto gap_open = EvalIntArg(context, arguments, 4);
		auto gap_extend = EvalIntArg(context, arguments, 5);
		auto gap_open2 = EvalIntArg(context, arguments, 6);
		auto noncan = EvalIntArg(context, arguments, 7);
		auto zdrop = EvalIntArg(context, arguments, 8);
		ValidatePenalties(match, mismatch, gap_open, gap_extend, gap_open2, noncan, zdrop);
		return make_uniq<AlignPairwiseKsw2SpliceBindData>(match, mismatch, gap_open, gap_extend, gap_open2, noncan,
		                                                  zdrop);
	}

	// Defaults match minimap2's --splice preset shape: match=2, mismatch=4, gap_open=6,
	// gap_extend=2, gap_open2=24 (intron open), noncan=9, zdrop=-1 (disabled).
	static unique_ptr<AlignPairwiseKsw2SpliceBindData> Defaults() {
		return make_uniq<AlignPairwiseKsw2SpliceBindData>(2, 4, 6, 2, 24, 9, -1);
	}
};

// ---------------------------------------------------------------------------
// Local state: per-thread KSW2Aligner constructed with splice params
// ---------------------------------------------------------------------------
struct AlignPairwiseKsw2SpliceLocalState : public FunctionLocalState {
	miint::KSW2Aligner aligner;
	std::string query_buf;
	std::string subject_buf;

	AlignPairwiseKsw2SpliceLocalState(int match, int mismatch, int gap_open, int gap_extend, int gap_open2, int noncan,
	                                  int zdrop)
	    : aligner(match, mismatch, gap_open, gap_extend, gap_open2, noncan, zdrop) {
		query_buf.reserve(256);
		subject_buf.reserve(256);
	}
};

static unique_ptr<FunctionLocalState> AlignPairwiseKsw2SpliceInitLocalState(ExpressionState &state,
                                                                            const BoundFunctionExpression &expr,
                                                                            FunctionData *bind_data) {
	auto &data = bind_data->Cast<AlignPairwiseKsw2SpliceBindData>();
	return make_uniq<AlignPairwiseKsw2SpliceLocalState>(data.match, data.mismatch, data.gap_open, data.gap_extend,
	                                                    data.gap_open2, data.noncan, data.zdrop);
}

// ---------------------------------------------------------------------------
// Argument-type vectors
// ---------------------------------------------------------------------------
// (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, noncan)
static vector<LogicalType> EightArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER};
}

// (query, subject, match, mismatch, gap_open, gap_extend, gap_open2, noncan, zdrop)
static vector<LogicalType> NineArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER};
}

// ---------------------------------------------------------------------------
// Execute bodies (shared templates from pairwise_align_shared.hpp)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseKsw2SpliceScoreExecute =
    RunPairwiseAlignScoreExecute<AlignPairwiseKsw2SpliceLocalState, &miint::KSW2Aligner::align_exts_score>;
static constexpr auto AlignPairwiseKsw2SpliceCigarExecute =
    RunPairwiseAlignCigarExecute<AlignPairwiseKsw2SpliceLocalState, &miint::KSW2Aligner::align_exts_cigar>;
static constexpr auto AlignPairwiseKsw2SpliceFullExecute =
    RunPairwiseAlignFullExecute<AlignPairwiseKsw2SpliceLocalState, &miint::KSW2Aligner::align_exts_full>;

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_splice_score → INTEGER
// ---------------------------------------------------------------------------
void AlignPairwiseKsw2SpliceScoreFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_splice_score");

	ScalarFunction score_2arg("align_pairwise_ksw2_splice_score", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          LogicalType::INTEGER, AlignPairwiseKsw2SpliceScoreExecute);
	score_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::Defaults().release());
	};
	score_2arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(score_2arg);

	ScalarFunction score_8arg("align_pairwise_ksw2_splice_score", EightArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseKsw2SpliceScoreExecute);
	score_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs8(ctx, args).release());
	};
	score_8arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(score_8arg);

	ScalarFunction score_9arg("align_pairwise_ksw2_splice_score", NineArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseKsw2SpliceScoreExecute);
	score_9arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_9arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs9(ctx, args).release());
	};
	score_9arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(score_9arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_splice_cigar → STRUCT(score INTEGER, cigar VARCHAR)
// ---------------------------------------------------------------------------
void AlignPairwiseKsw2SpliceCigarFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_splice_cigar");

	ScalarFunction cigar_2arg("align_pairwise_ksw2_splice_cigar", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          PairwiseCigarReturnType(), AlignPairwiseKsw2SpliceCigarExecute);
	cigar_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::Defaults().release());
	};
	cigar_2arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(cigar_2arg);

	ScalarFunction cigar_8arg("align_pairwise_ksw2_splice_cigar", EightArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseKsw2SpliceCigarExecute);
	cigar_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs8(ctx, args).release());
	};
	cigar_8arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(cigar_8arg);

	ScalarFunction cigar_9arg("align_pairwise_ksw2_splice_cigar", NineArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseKsw2SpliceCigarExecute);
	cigar_9arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_9arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs9(ctx, args).release());
	};
	cigar_9arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(cigar_9arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_ksw2_splice_full → STRUCT(score INTEGER, cigar VARCHAR,
//                                            query_aligned VARCHAR, subject_aligned VARCHAR)
// ---------------------------------------------------------------------------
void AlignPairwiseKsw2SpliceFullFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_ksw2_splice_full");

	ScalarFunction full_2arg("align_pairwise_ksw2_splice_full", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         PairwiseFullReturnType(), AlignPairwiseKsw2SpliceFullExecute);
	full_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::Defaults().release());
	};
	full_2arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(full_2arg);

	ScalarFunction full_8arg("align_pairwise_ksw2_splice_full", EightArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseKsw2SpliceFullExecute);
	full_8arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_8arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs8(ctx, args).release());
	};
	full_8arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(full_8arg);

	ScalarFunction full_9arg("align_pairwise_ksw2_splice_full", NineArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseKsw2SpliceFullExecute);
	full_9arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_9arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseKsw2SpliceBindData::FromArgs9(ctx, args).release());
	};
	full_9arg.init_local_state = AlignPairwiseKsw2SpliceInitLocalState;
	function_set.AddFunction(full_9arg);

	loader.RegisterFunction(function_set);
}

} // namespace duckdb
