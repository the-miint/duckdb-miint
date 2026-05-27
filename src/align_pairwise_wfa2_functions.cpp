#include "align_pairwise_wfa2_functions.hpp"

#include "WFA2Aligner.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "pairwise_align_shared.hpp"

namespace duckdb {

// ---------------------------------------------------------------------------
// Bind data: penalty constants validated at bind time
// ---------------------------------------------------------------------------
struct AlignPairwiseWfa2BindData : public FunctionData {
	int mismatch;
	int gap_open;
	int gap_extend;

	AlignPairwiseWfa2BindData(int mismatch_p, int gap_open_p, int gap_extend_p)
	    : mismatch(mismatch_p), gap_open(gap_open_p), gap_extend(gap_extend_p) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<AlignPairwiseWfa2BindData>(mismatch, gap_open, gap_extend);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<AlignPairwiseWfa2BindData>();
		return mismatch == other.mismatch && gap_open == other.gap_open && gap_extend == other.gap_extend;
	}

	// Penalty validation is delegated to WFA2Aligner constructor (single source of truth)
	// via the shared exception-translation helper.
	static void ValidatePenalties(int mismatch, int gap_open, int gap_extend) {
		TranslateCtorInvalidArgument([&] { miint::WFA2Aligner test(mismatch, gap_open, gap_extend); });
	}

	// Extract and validate penalty parameters from 5-arg overload bind arguments.
	// Arg layout: (query, subject, mismatch, gap_open, gap_extend)
	static unique_ptr<AlignPairwiseWfa2BindData> FromArgs5(ClientContext &context,
	                                                       vector<unique_ptr<Expression>> &arguments) {
		EnsurePenaltyArgsFoldable(arguments, 2, 5, "align_pairwise_wfa2_*");
		auto mismatch = EvalIntArg(context, arguments, 2);
		auto gap_open = EvalIntArg(context, arguments, 3);
		auto gap_extend = EvalIntArg(context, arguments, 4);
		ValidatePenalties(mismatch, gap_open, gap_extend);
		return make_uniq<AlignPairwiseWfa2BindData>(mismatch, gap_open, gap_extend);
	}

	static unique_ptr<AlignPairwiseWfa2BindData> Defaults() {
		return make_uniq<AlignPairwiseWfa2BindData>(4, 6, 2);
	}
};

// ---------------------------------------------------------------------------
// Local state: per-thread WFA2Aligner reused across all rows
// ---------------------------------------------------------------------------
struct AlignPairwiseWfa2LocalState : public FunctionLocalState {
	miint::WFA2Aligner aligner;
	// Reusable buffers to avoid per-row heap allocations in the execute loop.
	// Pre-reserved for typical Illumina read length (150bp).
	std::string query_buf;
	std::string subject_buf;

	AlignPairwiseWfa2LocalState(int mismatch, int gap_open, int gap_extend) : aligner(mismatch, gap_open, gap_extend) {
		query_buf.reserve(256);
		subject_buf.reserve(256);
	}
};

static unique_ptr<FunctionLocalState>
AlignPairwiseWfa2InitLocalState(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data) {
	auto &data = bind_data->Cast<AlignPairwiseWfa2BindData>();
	return make_uniq<AlignPairwiseWfa2LocalState>(data.mismatch, data.gap_open, data.gap_extend);
}

// 5-arg input types shared by all three functions: (query, subject, mismatch, gap_open, gap_extend)
static vector<LogicalType> FiveArgTypes() {
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER};
}

// ---------------------------------------------------------------------------
// align_pairwise_wfa2_score → INTEGER
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseWfa2ScoreExecute =
    RunPairwiseAlignScoreExecute<AlignPairwiseWfa2LocalState, &miint::WFA2Aligner::align_score>;

void AlignPairwiseWfa2ScoreFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_wfa2_score");

	ScalarFunction score_2arg("align_pairwise_wfa2_score", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          LogicalType::INTEGER, AlignPairwiseWfa2ScoreExecute);
	score_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::Defaults().release());
	};
	score_2arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(score_2arg);

	ScalarFunction score_5arg("align_pairwise_wfa2_score", FiveArgTypes(), LogicalType::INTEGER,
	                          AlignPairwiseWfa2ScoreExecute);
	score_5arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	score_5arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::FromArgs5(ctx, args).release());
	};
	score_5arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(score_5arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_wfa2_cigar → STRUCT(score INTEGER, cigar VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseWfa2CigarExecute =
    RunPairwiseAlignCigarExecute<AlignPairwiseWfa2LocalState, &miint::WFA2Aligner::align_cigar>;

void AlignPairwiseWfa2CigarFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_wfa2_cigar");

	ScalarFunction cigar_2arg("align_pairwise_wfa2_cigar", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                          PairwiseCigarReturnType(), AlignPairwiseWfa2CigarExecute);
	cigar_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::Defaults().release());
	};
	cigar_2arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(cigar_2arg);

	ScalarFunction cigar_5arg("align_pairwise_wfa2_cigar", FiveArgTypes(), PairwiseCigarReturnType(),
	                          AlignPairwiseWfa2CigarExecute);
	cigar_5arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	cigar_5arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseCigarReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::FromArgs5(ctx, args).release());
	};
	cigar_5arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(cigar_5arg);

	loader.RegisterFunction(function_set);
}

// ---------------------------------------------------------------------------
// align_pairwise_wfa2_full → STRUCT(score INTEGER, cigar VARCHAR,
//                                   query_aligned VARCHAR, subject_aligned VARCHAR)
// ---------------------------------------------------------------------------
static constexpr auto AlignPairwiseWfa2FullExecute =
    RunPairwiseAlignFullExecute<AlignPairwiseWfa2LocalState, &miint::WFA2Aligner::align_full>;

void AlignPairwiseWfa2FullFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("align_pairwise_wfa2_full");

	ScalarFunction full_2arg("align_pairwise_wfa2_full", {LogicalType::VARCHAR, LogicalType::VARCHAR},
	                         PairwiseFullReturnType(), AlignPairwiseWfa2FullExecute);
	full_2arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_2arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::Defaults().release());
	};
	full_2arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(full_2arg);

	ScalarFunction full_5arg("align_pairwise_wfa2_full", FiveArgTypes(), PairwiseFullReturnType(),
	                         AlignPairwiseWfa2FullExecute);
	full_5arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	full_5arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = PairwiseFullReturnType();
		return unique_ptr<FunctionData>(AlignPairwiseWfa2BindData::FromArgs5(ctx, args).release());
	};
	full_5arg.init_local_state = AlignPairwiseWfa2InitLocalState;
	function_set.AddFunction(full_5arg);

	loader.RegisterFunction(function_set);
}

} // namespace duckdb
