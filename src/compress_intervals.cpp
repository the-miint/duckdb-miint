#include "compress_intervals.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/aggregate_function.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"
#include <memory>

namespace duckdb {

struct IntervalState {
	std::unique_ptr<miint::IntervalCompressor> compressor;

	IntervalState() : compressor(nullptr) {
	}

	void Compress() {
		compressor->Compress();
	}

	void Add(int64_t start, int64_t stop) {
		compressor->Add(start, stop);
	}

	bool Empty() const {
		return compressor->Empty();
	}

	size_t Size() const {
		return compressor->Size();
	}

	const std::vector<int64_t> &Starts() const {
		return compressor->starts;
	}

	const std::vector<int64_t> &Stops() const {
		return compressor->stops;
	}
};

struct CompressIntervalsBindData : public FunctionData {
	explicit CompressIntervalsBindData() {
	}

	unique_ptr<FunctionData> Copy() const override {
		return make_uniq<CompressIntervalsBindData>();
	}

	bool Equals(const FunctionData &other) const override {
		return true;
	}
};

struct CompressIntervalsOperation {
	template <class STATE>
	static void Initialize(STATE &state) {
		new (&state) STATE();
		state.compressor = std::make_unique<miint::IntervalCompressor>();
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &aggr_input_data) {
		state.~STATE();
	}

	static void Operation(Vector inputs[], AggregateInputData &aggr_input_data, idx_t input_count, Vector &states,
	                      idx_t count) {
		auto &start_vector = inputs[0];
		auto &stop_vector = inputs[1];

		UnifiedVectorFormat start_data;
		UnifiedVectorFormat stop_data;
		start_vector.ToUnifiedFormat(count, start_data);
		stop_vector.ToUnifiedFormat(count, stop_data);

		auto start_ptr = UnifiedVectorFormat::GetData<int64_t>(start_data);
		auto stop_ptr = UnifiedVectorFormat::GetData<int64_t>(stop_data);

		UnifiedVectorFormat state_data;
		states.ToUnifiedFormat(count, state_data);
		auto state_ptr = UnifiedVectorFormat::GetData<IntervalState *>(state_data);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto start_idx = start_data.sel->get_index(i);
			auto stop_idx = stop_data.sel->get_index(i);

			if (!start_data.validity.RowIsValid(start_idx) || !stop_data.validity.RowIsValid(stop_idx)) {
				continue;
			}

			auto *state = state_ptr[state_idx];
			state->Add(start_ptr[start_idx], stop_ptr[stop_idx]);
		}
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &aggr_input_data) {
		for (idx_t i = 0; i < source.Size(); i++) {
			target.Add(source.Starts()[i], source.Stops()[i]);
		}
	}

	static void Finalize(Vector &state_vector, AggregateInputData &aggr_input_data, Vector &result, idx_t count,
	                     idx_t offset) {
		UnifiedVectorFormat state_data;
		state_vector.ToUnifiedFormat(count, state_data);
		auto states = UnifiedVectorFormat::GetData<IntervalState *>(state_data);

		auto &result_validity = FlatVector::Validity(result);
		auto result_data = FlatVector::GetData<list_entry_t>(result);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto &state = *states[state_idx];

			state.Compress();

			if (state.Empty()) {
				result_validity.SetInvalid(i + offset);
				continue;
			}

			auto &list_entry = ListVector::GetEntry(result);
			auto list_offset = ListVector::GetListSize(result);
			ListVector::Reserve(result, list_offset + state.Size());

			auto &struct_children = StructVector::GetEntries(list_entry);
			auto &start_child = struct_children[0];
			auto &stop_child = struct_children[1];

			auto start_ptr = FlatVector::GetData<int64_t>(*start_child);
			auto stop_ptr = FlatVector::GetData<int64_t>(*stop_child);

			for (idx_t j = 0; j < state.Size(); j++) {
				start_ptr[list_offset + j] = state.Starts()[j];
				stop_ptr[list_offset + j] = state.Stops()[j];
			}

			ListVector::SetListSize(result, list_offset + state.Size());

			result_data[i + offset].offset = list_offset;
			result_data[i + offset].length = state.Size();
		}
	}

	static bool IgnoreNull() {
		return true;
	}
};

void CompressIntervalsFunction::Register(ExtensionLoader &loader) {
	auto fun = AggregateFunction(
	    "compress_intervals", {LogicalType::BIGINT, LogicalType::BIGINT},
	    LogicalType::LIST(LogicalType::STRUCT({{"start", LogicalType::BIGINT}, {"stop", LogicalType::BIGINT}})),
	    AggregateFunction::StateSize<IntervalState>,
	    AggregateFunction::StateInitialize<IntervalState, CompressIntervalsOperation>,
	    CompressIntervalsOperation::Operation,
	    AggregateFunction::StateCombine<IntervalState, CompressIntervalsOperation>,
	    CompressIntervalsOperation::Finalize, nullptr, nullptr,
	    AggregateFunction::StateDestroy<IntervalState, CompressIntervalsOperation>);

	loader.RegisterFunction(fun);
}

} // namespace duckdb
