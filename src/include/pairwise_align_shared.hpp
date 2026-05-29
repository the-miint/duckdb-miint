#pragma once

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/planner/expression.hpp"

#include <stdexcept>
#include <string>

namespace duckdb {

// Return types shared by every align_pairwise_* family. The _cigar and _full STRUCT
// shapes are identical across WFA2 + KSW2 (extz/extd/exts) — promoted here to keep the
// LogicalType definition in one place. Defined as inline functions (not constants) because
// LogicalType construction is not constexpr.
inline LogicalType PairwiseCigarReturnType() {
	return LogicalType::STRUCT({{"score", LogicalType::INTEGER}, {"cigar", LogicalType::VARCHAR}});
}

inline LogicalType PairwiseFullReturnType() {
	return LogicalType::STRUCT({{"score", LogicalType::INTEGER},
	                            {"cigar", LogicalType::VARCHAR},
	                            {"query_aligned", LogicalType::VARCHAR},
	                            {"subject_aligned", LogicalType::VARCHAR}});
}

// Throw InvalidInputException if any argument in [start, end) is not foldable (i.e. not a
// constant). Used by every align_pairwise_* FromArgsN to enforce that penalty params are
// bind-time constants, not column references. `family_name` is the user-visible function
// prefix (e.g. "align_pairwise_wfa2_*") used in the error message.
inline void EnsurePenaltyArgsFoldable(const vector<unique_ptr<Expression>> &arguments, idx_t start, idx_t end,
                                      const char *family_name) {
	for (idx_t i = start; i < end; i++) {
		if (!arguments[i]->IsFoldable()) {
			throw InvalidInputException(std::string(family_name) +
			                            " penalty parameters must be constant values, not column references");
		}
	}
}

// Evaluate one bind-time integer argument. Wraps the EvaluateScalar+GetValue<int>() pattern
// so each FromArgsN body reads as a flat list of `auto x = EvalIntArg(ctx, args, k);` lines.
inline int EvalIntArg(ClientContext &context, vector<unique_ptr<Expression>> &arguments, idx_t i) {
	return ExpressionExecutor::EvaluateScalar(context, *arguments[i]).GetValue<int>();
}

// Run an Aligner constructor and translate its std::invalid_argument (the C++-API contract
// for invalid penalty values) into DuckDB's InvalidInputException. Keeps validation as a
// single source of truth (the Aligner ctor) while surfacing the right exception type to SQL.
template <typename Func>
inline void TranslateCtorInvalidArgument(Func &&ctor_call) {
	try {
		ctor_call();
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException(e.what());
	}
}

// Prepared input vectors for a (query VARCHAR, subject VARCHAR) scalar-function execute loop.
// No raw pointers held: row-level access derives `string_t *` from the UnifiedVectorFormat
// on demand inside GetAlignInput, so the lifetime of every pointer is bounded by a single row.
struct AlignInputVectors {
	UnifiedVectorFormat query_data;
	UnifiedVectorFormat subject_data;
};

// Unify the first two args of the chunk into AlignInputVectors. Caller's responsibility to
// ensure args.data[0] and args.data[1] are VARCHAR.
inline AlignInputVectors PrepareAlignInputs(DataChunk &args) {
	AlignInputVectors v;
	args.data[0].ToUnifiedFormat(args.size(), v.query_data);
	args.data[1].ToUnifiedFormat(args.size(), v.subject_data);
	return v;
}

// Populate reusable buffers from input vectors at row i. Returns false if either input is NULL.
inline bool GetAlignInput(AlignInputVectors &v, idx_t i, std::string &query_buf, std::string &subject_buf) {
	auto qi = v.query_data.sel->get_index(i);
	auto si = v.subject_data.sel->get_index(i);
	if (!v.query_data.validity.RowIsValid(qi) || !v.subject_data.validity.RowIsValid(si)) {
		return false;
	}
	auto q_ptr = UnifiedVectorFormat::GetData<string_t>(v.query_data);
	auto s_ptr = UnifiedVectorFormat::GetData<string_t>(v.subject_data);
	query_buf.assign(q_ptr[qi].GetData(), q_ptr[qi].GetSize());
	subject_buf.assign(s_ptr[si].GetData(), s_ptr[si].GetSize());
	return true;
}

// Templated scalar-function execute bodies shared across all pairwise-alignment families
// (WFA2 + KSW2 extz/extd/exts). Each family supplies a LocalState type (which must expose
// `aligner`, `query_buf`, `subject_buf` members) and a pointer-to-member-function on the aligner.
// The output-shape differs per detail level (score / cigar / full), so three templates are
// supplied rather than one.

template <typename LocalState, auto Method>
inline void RunPairwiseAlignScoreExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->template Cast<LocalState>();
	auto inputs = PrepareAlignInputs(args);
	auto result_data = FlatVector::GetData<int32_t>(result);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!GetAlignInput(inputs, i, lstate.query_buf, lstate.subject_buf)) {
			result_validity.SetInvalid(i);
			continue;
		}
		auto score = (lstate.aligner.*Method)(lstate.query_buf, lstate.subject_buf);
		if (!score.has_value()) {
			result_validity.SetInvalid(i);
			continue;
		}
		result_data[i] = score.value();
	}
}

template <typename LocalState, auto Method>
inline void RunPairwiseAlignCigarExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->template Cast<LocalState>();
	auto inputs = PrepareAlignInputs(args);

	auto &entries = StructVector::GetEntries(result);
	auto score_data = FlatVector::GetData<int32_t>(*entries[0]);
	auto &cigar_vec = *entries[1];
	auto cigar_data = FlatVector::GetData<string_t>(cigar_vec);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!GetAlignInput(inputs, i, lstate.query_buf, lstate.subject_buf)) {
			result_validity.SetInvalid(i);
			continue;
		}
		auto cigar_result = (lstate.aligner.*Method)(lstate.query_buf, lstate.subject_buf);
		if (!cigar_result.has_value()) {
			result_validity.SetInvalid(i);
			continue;
		}
		score_data[i] = cigar_result->score;
		cigar_data[i] = StringVector::AddString(cigar_vec, cigar_result->cigar);
	}
}

template <typename LocalState, auto Method>
inline void RunPairwiseAlignFullExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->template Cast<LocalState>();
	auto inputs = PrepareAlignInputs(args);

	auto &entries = StructVector::GetEntries(result);
	auto score_data = FlatVector::GetData<int32_t>(*entries[0]);
	auto &cigar_vec = *entries[1];
	auto &query_aligned_vec = *entries[2];
	auto &subject_aligned_vec = *entries[3];
	auto cigar_data = FlatVector::GetData<string_t>(cigar_vec);
	auto query_aligned_data = FlatVector::GetData<string_t>(query_aligned_vec);
	auto subject_aligned_data = FlatVector::GetData<string_t>(subject_aligned_vec);
	auto &result_validity = FlatVector::Validity(result);

	for (idx_t i = 0; i < args.size(); i++) {
		if (!GetAlignInput(inputs, i, lstate.query_buf, lstate.subject_buf)) {
			result_validity.SetInvalid(i);
			continue;
		}
		auto full_result = (lstate.aligner.*Method)(lstate.query_buf, lstate.subject_buf);
		if (!full_result.has_value()) {
			result_validity.SetInvalid(i);
			continue;
		}
		score_data[i] = full_result->score;
		cigar_data[i] = StringVector::AddString(cigar_vec, full_result->cigar);
		query_aligned_data[i] = StringVector::AddString(query_aligned_vec, full_result->query_aligned);
		subject_aligned_data[i] = StringVector::AddString(subject_aligned_vec, full_result->subject_aligned);
	}
}

} // namespace duckdb
