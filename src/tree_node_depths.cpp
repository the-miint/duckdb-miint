#include "tree_node_depths.hpp"

#include "TreeIndex.hpp"
#include "catalog_utils.hpp"
#include "tree_function_utils.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <algorithm>
#include <atomic>

namespace duckdb {

namespace {

constexpr const char *kFn = "tree_node_depths";

struct BindData : public TableFunctionData {
	std::string tree_table;
};

struct GlobalState : public GlobalTableFunctionState {
	IndexedTree tree;
	std::atomic<idx_t> next_node {0};
	idx_t max_threads = 1;

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<BindData>();
	data->tree_table = input.inputs[0].ToString();
	RejectCTERelationName(input, data->tree_table);
	// Schemas only: no rows are read until InitGlobal (docs/internals/architecture.md).
	ValidateTreeTableSchema(context, data->tree_table);

	names.emplace_back("node_index");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("depth");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("distance_to_root");
	return_types.push_back(LogicalType::DOUBLE);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<BindData>();
	auto state = make_uniq<GlobalState>();
	// No names are needed, so none are read, and there is no name type to check.
	TipMatches no_names;
	state->tree = ReadIndexedTree(context, kFn, bind.tree_table, no_names, LogicalType::INVALID,
	                              miint::LengthPolicy::UnknownAllowed, "read");

	const idx_t blocks = (state->tree.index.size() + STANDARD_VECTOR_SIZE - 1) / STANDARD_VECTOR_SIZE;
	const idx_t threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	state->max_threads = std::max<idx_t>(1, std::min(blocks, threads));
	return std::move(state);
}

// Each call answers the next block of nodes; threads claim blocks independently, so the
// rows come out in no particular order.
void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &state = data_p.global_state->Cast<GlobalState>();
	const auto &index = state.tree.index;
	const idx_t start = state.next_node.fetch_add(STANDARD_VECTOR_SIZE);
	if (start >= index.size()) {
		output.SetChildCardinality(0);
		return;
	}
	const idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, index.size() - start);

	auto node = FlatVector::GetDataMutable<int64_t>(output.data[0]);
	auto depth = FlatVector::GetDataMutable<int64_t>(output.data[1]);
	auto distance = FlatVector::GetDataMutable<double>(output.data[2]);
	auto &distance_valid = FlatVector::ValidityMutable(output.data[2]);
	for (idx_t i = 0; i < count; ++i) {
		const auto v = static_cast<uint32_t>(start + i);
		node[i] = state.tree.node_ids[v];
		depth[i] = index.depth(v);
		const auto d = index.distance_to_ancestor(v, index.root());
		if (d) {
			distance[i] = *d;
			distance_valid.SetValid(i);
		} else {
			distance_valid.SetInvalid(i);
		}
	}
	output.SetChildCardinality(count);
}

} // namespace

TableFunction TreeNodeDepthsTableFunction::GetFunction() {
	TableFunction tf("tree_node_depths", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void TreeNodeDepthsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
