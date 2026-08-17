#include "cumulative_coverage.hpp"
#include "alignment_functions_internal.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/function/aggregate_function.hpp"
#include "duckdb/planner/expression/bound_aggregate_expression.hpp"

namespace duckdb {

// cumulative_coverage(rank, start, stop)
//   -> LIST<STRUCT(rank INTEGER, covered BIGINT)>
//
// Cumulative breadth of coverage in rank order (issue #214): element k reports how many
// bases are covered by the union of every interval with rank <= k. The algorithm lives
// in miint::CumulativeCoverageAccumulator; this is only the DuckDB glue.
//
// The rank is computed by the CALLER in SQL, not here, which is the main departure from
// the shape #214 proposed. That proposal folded sample_id, the breadth ranking, the tie
// policy and a genome_length constant-check into the aggregate. Pushing the ranking out
// to a window function removes the type-mirroring bind, the heap-allocated id keys in
// the aggregate state, the tiebreak policy and the constant-check -- and lets callers
// rank by anything, not only breadth. cumulative_coverage_curve() in miint_macros.hpp
// is the wrapper that does the ranking and reattaches sample_id.
//
// A row with a NULL start or stop registers its rank with no coverage. That is the
// zero-coverage backfill, and it is load-bearing rather than a convenience: a sample in
// the group with no coverage of the target produces no interval rows at all, so without
// an explicit signal it would vanish, understating the group size and making curves
// from groups with different detection rates incomparable. The caller produces those
// rows with an ordinary LEFT JOIN from the roster.
//
// That is also why this aggregate must see NULL rows. It supplies its own Operation
// rather than going through the AggregateExecutor::Unary/Binary templates, and
// IgnoreNull() is only consulted inside those templates -- so NULL rows arrive here
// regardless. compress_intervals relies on the same fact in the opposite direction, and
// has an explicit RowIsValid check because of it. IgnoreNull() below returns false to
// state the requirement rather than to change behaviour.
//
// A NULL rank, by contrast, raises. ROW_NUMBER() cannot produce one, so a NULL means
// the rank column came from somewhere else and the curve's x-axis is not what the
// caller thinks it is.

struct CumulativeCoverageState {
	miint::CumulativeCoverageAccumulator *accumulator;

	CumulativeCoverageState() : accumulator(nullptr) {
	}
};

struct CumulativeCoverageOperation {
	template <class STATE>
	static void Initialize(STATE &state) {
		state.accumulator = new miint::CumulativeCoverageAccumulator();
	}

	template <class STATE>
	static void Destroy(STATE &state, AggregateInputData &aggr_input_data) {
		delete state.accumulator;
	}

	static void Operation(Vector inputs[], AggregateInputData &aggr_input_data, idx_t input_count, Vector &states,
	                      idx_t count) {
		auto &rank_vector = inputs[0];
		auto &start_vector = inputs[1];
		auto &stop_vector = inputs[2];

		UnifiedVectorFormat rank_data;
		UnifiedVectorFormat start_data;
		UnifiedVectorFormat stop_data;
		rank_vector.ToUnifiedFormat(count, rank_data);
		start_vector.ToUnifiedFormat(count, start_data);
		stop_vector.ToUnifiedFormat(count, stop_data);

		auto rank_ptr = UnifiedVectorFormat::GetData<int32_t>(rank_data);
		auto start_ptr = UnifiedVectorFormat::GetData<int64_t>(start_data);
		auto stop_ptr = UnifiedVectorFormat::GetData<int64_t>(stop_data);

		UnifiedVectorFormat state_data;
		states.ToUnifiedFormat(count, state_data);
		auto state_ptr = UnifiedVectorFormat::GetData<CumulativeCoverageState *>(state_data);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto rank_idx = rank_data.sel->get_index(i);
			auto start_idx = start_data.sel->get_index(i);
			auto stop_idx = stop_data.sel->get_index(i);

			if (!rank_data.validity.RowIsValid(rank_idx)) {
				throw InvalidInputException(
				    "cumulative_coverage: rank must not be NULL. Produce ranks with ROW_NUMBER() OVER (...) - 1; a "
				    "NULL rank usually means the ranking window returned no row for this sample.");
			}

			auto *state = state_ptr[state_idx];
			const int32_t rank = rank_ptr[rank_idx];

			const bool start_null = !start_data.validity.RowIsValid(start_idx);
			const bool stop_null = !stop_data.validity.RowIsValid(stop_idx);

			// Exactly one NULL is a data error, not the backfill signal. A row like
			// (rank, NULL, 100) comes from a join that matched one coordinate column but
			// not the other, or a blank field in one column; treating it as "covers
			// nothing" would discard real coverage and understate this rank and every
			// rank above it. Both-NULL is the documented signal and is handled below.
			if (start_null != stop_null) {
				throw InvalidInputException(
				    "cumulative_coverage: start and stop must both be NULL (a rank with no coverage) or both be "
				    "present; rank %d has %s. One NULL is a broken join or a blank field.",
				    rank, start_null ? "start NULL and stop present" : "start present and stop NULL");
			}

			// Both NULL == "this rank exists and covers nothing".
			if (start_null) {
				state->accumulator->AddEmptyRank(rank);
				continue;
			}

			state->accumulator->Add(rank, start_ptr[start_idx], stop_ptr[stop_idx]);
		}
	}

	template <class STATE, class OP>
	static void Combine(const STATE &source, STATE &target, AggregateInputData &aggr_input_data) {
		target.accumulator->Absorb(*source.accumulator);
	}

	static void Finalize(Vector &state_vector, AggregateInputData &aggr_input_data, Vector &result, idx_t count,
	                     idx_t offset) {
		UnifiedVectorFormat state_data;
		state_vector.ToUnifiedFormat(count, state_data);
		auto states = UnifiedVectorFormat::GetData<CumulativeCoverageState *>(state_data);

		auto &result_validity = FlatVector::Validity(result);
		auto result_data = FlatVector::GetData<list_entry_t>(result);

		for (idx_t i = 0; i < count; i++) {
			auto state_idx = state_data.sel->get_index(i);
			auto &state = *states[state_idx];

			if (state.accumulator->Empty()) {
				result_validity.SetInvalid(i + offset);
				continue;
			}

			// Rank contiguity is validated here rather than in Operation: under a
			// parallel GROUP BY each thread sees an arbitrary subset, so no single
			// Operation call can tell a genuine gap from a rank another thread holds.
			std::vector<miint::CumulativePoint> curve;
			try {
				curve = state.accumulator->Curve();
			} catch (const miint::InvalidInputException &e) {
				throw InvalidInputException("%s", e.what());
			}

			auto &list_entry = ListVector::GetEntry(result);
			auto list_offset = ListVector::GetListSize(result);
			ListVector::Reserve(result, list_offset + curve.size());

			auto &struct_children = StructVector::GetEntries(list_entry);
			auto rank_ptr = FlatVector::GetData<int32_t>(*struct_children[0]);
			auto covered_ptr = FlatVector::GetData<int64_t>(*struct_children[1]);

			for (idx_t j = 0; j < curve.size(); j++) {
				rank_ptr[list_offset + j] = curve[j].rank;
				covered_ptr[list_offset + j] = curve[j].covered;
			}

			ListVector::SetListSize(result, list_offset + curve.size());

			result_data[i + offset].offset = list_offset;
			result_data[i + offset].length = curve.size();
		}
	}

	// See the header note: NULL start/stop rows are DATA here, not absence, so they
	// must not be filtered out before Operation sees them.
	static bool IgnoreNull() {
		return false;
	}
};

void CumulativeCoverageFunction::Register(ExtensionLoader &loader) {
	auto fun = AggregateFunction(
	    "cumulative_coverage", {LogicalType::INTEGER, LogicalType::BIGINT, LogicalType::BIGINT},
	    LogicalType::LIST(LogicalType::STRUCT({{"rank", LogicalType::INTEGER}, {"covered", LogicalType::BIGINT}})),
	    AggregateFunction::StateSize<CumulativeCoverageState>,
	    AggregateFunction::StateInitialize<CumulativeCoverageState, CumulativeCoverageOperation>,
	    CumulativeCoverageOperation::Operation,
	    AggregateFunction::StateCombine<CumulativeCoverageState, CumulativeCoverageOperation>,
	    CumulativeCoverageOperation::Finalize, nullptr, nullptr,
	    AggregateFunction::StateDestroy<CumulativeCoverageState, CumulativeCoverageOperation>);

	loader.RegisterFunction(fun);
}

} // namespace duckdb
