#include "compute_coverage_depth.hpp"
#include "alignment_functions_internal.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/aggregate_function.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"

namespace duckdb {

struct CoverageDepthState {
	miint::CoverageDepthCalculator *calculator;
	int64_t reference_length;
	bool initialized;

	CoverageDepthState() : calculator(nullptr), reference_length(0), initialized(false) {
	}

	const std::vector<uint32_t> &GetDepths() const {
		return calculator->GetDepths();
	}

	bool Empty() const {
		return !initialized || calculator->Empty();
	}

	int64_t Size() const {
		return reference_length;
	}
};

struct CoverageDepthBindData : public FunctionData {
	bool include_deletions;

	explicit CoverageDepthBindData(bool include_deletions) : include_deletions(include_deletions) {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<CoverageDepthBindData>(include_deletions);
	}

	bool Equals(const FunctionData &other) const override {
		auto &other_data = other.Cast<CoverageDepthBindData>();
		return include_deletions == other_data.include_deletions;
	}
};

struct CoverageDepthOperation {
	template <class STATE>
	static void Initialize(STATE &state) {
		state.calculator = nullptr;
		state.reference_length = 0;
		state.initialized = false;
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &aggr_input_data) {
		delete state.calculator;
	}

	static void Operation(Vector inputs[], AggregateInputData &aggr_input_data, idx_t input_count, Vector &states,
	                      idx_t count) {
		auto &bind_data = aggr_input_data.bind_data->Cast<CoverageDepthBindData>();

		auto &position_vector = inputs[0];
		auto &stop_position_vector = inputs[1];
		auto &cigar_vector = inputs[2];
		auto &ref_length_vector = inputs[3];

		UnifiedVectorFormat position_data, stop_data, cigar_data, ref_length_data, state_data;
		position_vector.ToUnifiedFormat(count, position_data);
		stop_position_vector.ToUnifiedFormat(count, stop_data);
		cigar_vector.ToUnifiedFormat(count, cigar_data);
		ref_length_vector.ToUnifiedFormat(count, ref_length_data);
		states.ToUnifiedFormat(count, state_data);

		auto position_ptr = UnifiedVectorFormat::GetData<int64_t>(position_data);
		auto stop_ptr = UnifiedVectorFormat::GetData<int64_t>(stop_data);
		auto cigar_ptr = UnifiedVectorFormat::GetData<string_t>(cigar_data);
		auto ref_length_ptr = UnifiedVectorFormat::GetData<int64_t>(ref_length_data);
		auto state_ptr = UnifiedVectorFormat::GetData<CoverageDepthState *>(state_data);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto pos_idx = position_data.sel->get_index(i);
			auto stop_idx = stop_data.sel->get_index(i);
			auto cigar_idx = cigar_data.sel->get_index(i);
			auto ref_len_idx = ref_length_data.sel->get_index(i);

			if (!position_data.validity.RowIsValid(pos_idx) || !stop_data.validity.RowIsValid(stop_idx) ||
			    !cigar_data.validity.RowIsValid(cigar_idx) || !ref_length_data.validity.RowIsValid(ref_len_idx)) {
				continue;
			}

			auto *state = state_ptr[state_idx];
			auto ref_length = ref_length_ptr[ref_len_idx];

			if (ref_length <= 0) {
				throw InvalidInputException("compute_coverage_depth: reference_length must be positive, got %lld",
				                            ref_length);
			}
			if (ref_length > 2'000'000'000LL) {
				throw InvalidInputException(
				    "compute_coverage_depth: reference_length %lld exceeds maximum (2,000,000,000). "
				    "Consider using compress_intervals for large references.",
				    ref_length);
			}

			if (!state->initialized) {
				state->reference_length = ref_length;
				state->calculator =
				    new miint::CoverageDepthCalculator(ref_length, bind_data.include_deletions);
				state->initialized = true;
			} else if (state->reference_length != ref_length) {
				throw InvalidInputException(
				    "compute_coverage_depth: reference_length must be the same for all rows in a group, "
				    "got %lld and %lld",
				    state->reference_length, ref_length);
			}

			auto position = position_ptr[pos_idx];
			auto stop_position = stop_ptr[stop_idx];
			auto cigar = cigar_ptr[cigar_idx].GetString();

			try {
				state->calculator->AddRead(position, stop_position, cigar);
			} catch (const miint::InvalidInputException &e) {
				throw InvalidInputException("%s", e.what());
			}
		}
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &aggr_input_data) {
		if (!source.initialized) {
			return;
		}
		if (!target.initialized) {
			target.reference_length = source.reference_length;
			target.calculator = new miint::CoverageDepthCalculator(
			    source.reference_length, aggr_input_data.bind_data->Cast<CoverageDepthBindData>().include_deletions);
			target.initialized = true;
		} else if (target.reference_length != source.reference_length) {
			throw InvalidInputException(
			    "compute_coverage_depth: reference_length must be the same for all rows in a group, "
			    "got %lld and %lld",
			    target.reference_length, source.reference_length);
		}
		try {
			target.calculator->Combine(*source.calculator);
		} catch (const miint::InvalidInputException &e) {
			throw InvalidInputException("%s", e.what());
		}
	}

	static void Finalize(Vector &state_vector, AggregateInputData &aggr_input_data, Vector &result, idx_t count,
	                     idx_t offset) {
		UnifiedVectorFormat state_data;
		state_vector.ToUnifiedFormat(count, state_data);
		auto states = UnifiedVectorFormat::GetData<CoverageDepthState *>(state_data);

		auto &result_validity = FlatVector::Validity(result);
		auto result_data = FlatVector::GetData<list_entry_t>(result);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto &state = *states[state_idx];

			if (state.Empty()) {
				result_validity.SetInvalid(i + offset);
				continue;
			}

			auto &depths = state.GetDepths();
			auto list_size = static_cast<idx_t>(depths.size());

			auto &list_entry = ListVector::GetEntry(result);
			auto list_offset = ListVector::GetListSize(result);
			ListVector::Reserve(result, list_offset + list_size);

			auto child_ptr = FlatVector::GetData<uint32_t>(list_entry);

			for (idx_t j = 0; j < list_size; j++) {
				child_ptr[list_offset + j] = depths[j];
			}

			ListVector::SetListSize(result, list_offset + list_size);

			result_data[i + offset].offset = list_offset;
			result_data[i + offset].length = list_size;
		}
	}

	static bool IgnoreNull() {
		return true;
	}
};

void ComputeCoverageDepthFunction::Register(ExtensionLoader &loader) {
	auto fun = AggregateFunction(
	    "compute_coverage_depth",
	    {LogicalType::BIGINT, LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::VARCHAR},
	    LogicalType::LIST(LogicalType::UINTEGER), AggregateFunction::StateSize<CoverageDepthState>,
	    AggregateFunction::StateInitialize<CoverageDepthState, CoverageDepthOperation>,
	    CoverageDepthOperation::Operation,
	    AggregateFunction::StateCombine<CoverageDepthState, CoverageDepthOperation>,
	    CoverageDepthOperation::Finalize, nullptr, nullptr,
	    AggregateFunction::StateDestroy<CoverageDepthState, CoverageDepthOperation>);

	// Validate mode parameter at bind time
	fun.bind = [](ClientContext &context, AggregateFunction &function,
	              vector<unique_ptr<Expression>> &arguments) -> unique_ptr<FunctionData> {
		if (arguments[4]->IsFoldable()) {
			auto mode_val = ExpressionExecutor::EvaluateScalar(context, *arguments[4]);
			if (!mode_val.IsNull()) {
				auto mode = mode_val.ToString();
				if (mode == "include_deletions") {
					return make_uniq<CoverageDepthBindData>(true);
				} else if (mode == "exclude_deletions") {
					return make_uniq<CoverageDepthBindData>(false);
				} else {
					throw InvalidInputException(
					    "Invalid mode '%s' for compute_coverage_depth. Must be 'include_deletions' or "
					    "'exclude_deletions'.",
					    mode);
				}
			}
		}
		throw InvalidInputException(
		    "compute_coverage_depth: mode parameter must be a constant string "
		    "('include_deletions' or 'exclude_deletions')");
	};

	loader.RegisterFunction(fun);
}

} // namespace duckdb
