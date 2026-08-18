#include "compress_intervals.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/aggregate_function.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"

namespace duckdb {

struct IntervalState {
	miint::IntervalCompressor *compressor;

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
		state.compressor = new miint::IntervalCompressor();
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &aggr_input_data) {
		delete state.compressor;
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
		// Deliberately NO Compress() here. Add() checks IntervalCompressor's growing
		// compression floor on every push, so the state is already bounded when this loop
		// ends -- and Finalize() compresses each state before reading Size()/Empty(), so the
		// output is unaffected either way. An unconditional compression here therefore
		// reclaimed nothing on disjoint input and cost O(m log m) PER COMBINE. It could also
		// compress twice in a row, when the Add loop above had just crossed the floor itself.
		//
		// Disjoint intervals are exactly what this aggregate is normally fed, and Compress()
		// merges only overlaps, so "compaction reclaims nothing" is the common case rather
		// than a corner. See compress_floor_ in IntervalCompressor.
		//
		// Measured on 500k rows / 50 groups of disjoint intervals, results bit-identical
		// (500000 intervals, 2500000 total width at both thread counts):
		//
		//   threads=1    0.016 s -> 0.012 s
		//   threads=16   0.025 s -> 0.019 s
		//
		// So about 25% off this shape at either thread count. Note what this does NOT fix:
		// threads=16 remains slower than threads=1 here, and still is at 4M rows / 200
		// groups (0.112 s vs 0.120 s). That residual is not the merge policy -- Combine runs
		// at threads=1 too, since the hash aggregate builds partials per radix partition --
		// so do not read this as making the aggregate scale with cores.
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
