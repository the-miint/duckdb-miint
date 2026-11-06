#include "tree_ancestors.hpp"

#include "TreeIndex.hpp"
#include "catalog_utils.hpp"
#include "miint_named_parameter.hpp"
#include "tree_function_utils.hpp"
#include "tree_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <algorithm>
#include <atomic>

namespace duckdb {

namespace {

constexpr const char *kFn = "tree_ancestors";

// Nodes are claimed a few at a time: a deep node has thousands of ancestors, and small
// claims spread those nodes over the threads.
constexpr idx_t kNodesPerClaim = 64;

struct BindData : public TableFunctionData {
	std::string tree_table;
	bool has_nodes_table = false;
	std::string nodes_table;
};

struct GlobalState : public GlobalTableFunctionState {
	IndexedTree tree;
	// The nodes whose ancestors are listed: nodes_table's, or every node when there is none.
	bool all_nodes = true;
	std::vector<uint32_t> nodes;
	idx_t node_count = 0;
	std::atomic<idx_t> next_claim {0};
	idx_t max_threads = 1;

	idx_t MaxThreads() const override {
		return max_threads;
	}
	uint32_t NodeAt(idx_t i) const {
		return all_nodes ? static_cast<uint32_t>(i) : nodes[i];
	}
};

// Where a thread is: the positions [next, end) of its claim still to start, and the node
// whose ancestors it is listing, the next of them and how many edges up that one is. A
// node's ancestors can run past an output chunk; the next call carries on from here.
struct LocalState : public LocalTableFunctionState {
	idx_t next = 0;
	idx_t end = 0;
	uint32_t node = 0;
	uint32_t ancestor = miint::TreeIndex::NO_PARENT;
	int64_t edges = 0;
};

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<Identifier> &names) {
	auto data = make_uniq<BindData>();
	data->tree_table = input.inputs[0].ToString();
	RejectCTERelationName(input, data->tree_table);
	auto nodes = input.named_parameters.find(Identifier("nodes_table"));
	if (nodes != input.named_parameters.end() && !nodes->second.IsNull()) {
		data->has_nodes_table = true;
		data->nodes_table = StringValue::Get(nodes->second);
		RejectCTERelationName(input, data->nodes_table);
	}

	// Schemas only: no rows are read until InitGlobal (docs/internals/architecture.md).
	ValidateTreeTableSchema(context, data->tree_table);
	if (data->has_nodes_table) {
		// The tree is read in its own projection, so passing it as nodes_table would read it twice.
		if (StringUtil::CIEquals(data->nodes_table, data->tree_table)) {
			throw BinderException("%s: '%s' is the tree table; leave nodes_table out for every node, or pass the nodes "
			                      "as a relation of their own",
			                      kFn, data->nodes_table);
		}
		auto columns = GetTableOrViewColumns(context, data->nodes_table, "Nodes table");
		LogicalType type;
		if (!FindColumnType(columns, "node_index", type)) {
			throw BinderException("%s: nodes table '%s' has no 'node_index' column", kFn, data->nodes_table);
		}
		if (!type.IsIntegral()) {
			throw BinderException("%s: node_index in '%s' must be an integer type, but it is %s", kFn,
			                      data->nodes_table, type.ToString());
		}
	}

	names.emplace_back("node_index");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("ancestor_index");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("edges");
	return_types.push_back(LogicalType::BIGINT);
	names.emplace_back("distance");
	return_types.push_back(LogicalType::DOUBLE);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bind = input.bind_data->Cast<BindData>();
	auto state = make_uniq<GlobalState>();

	// nodes_table first, once: its distinct ids, ascending, NULLs dropped. Any integer type
	// is read as BIGINT; a value out of BIGINT's range fails the read.
	std::vector<int64_t> ids;
	if (bind.has_nodes_table) {
		ReadRelation(context, kFn, bind.nodes_table, {{"CAST(node_index AS BIGINT)", LogicalType::BIGINT}},
		             [&](DataChunk &chunk) {
			             UnifiedVectorFormat format;
			             chunk.data[0].ToUnifiedFormat(format);
			             auto values = UnifiedVectorFormat::GetData<int64_t>(format);
			             for (idx_t i = 0; i < chunk.size(); ++i) {
				             const auto k = format.sel->get_index(i);
				             if (format.validity.RowIsValid(k)) {
					             ids.push_back(values[k]);
				             }
			             }
		             });
		std::sort(ids.begin(), ids.end());
		ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
	}

	// No names are needed, so none are read, and there is no name type to check.
	TipMatches no_names;
	state->tree = ReadIndexedTree(context, kFn, bind.tree_table, no_names, LogicalType::INVALID,
	                              miint::LengthPolicy::UnknownAllowed, "read");

	if (bind.has_nodes_table) {
		std::vector<std::string> missing;
		state->nodes.reserve(ids.size());
		for (const int64_t id : ids) {
			const uint32_t v = miint::DenseTree::Find(state->tree.node_ids, id);
			if (v == miint::DenseTree::NOT_FOUND) {
				missing.push_back(std::to_string(id));
			} else {
				state->nodes.push_back(v);
			}
		}
		if (!missing.empty()) {
			throw InvalidInputException("%s: %d node_index value(s) in '%s' are not in the tree '%s': %s", kFn,
			                            static_cast<int64_t>(missing.size()), bind.nodes_table, bind.tree_table,
			                            QuoteSome(missing));
		}
		state->all_nodes = false;
		state->node_count = state->nodes.size();
	} else {
		state->node_count = state->tree.index.size();
	}

	const idx_t claims = (state->node_count + kNodesPerClaim - 1) / kNodesPerClaim;
	const idx_t threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	state->max_threads = std::max<idx_t>(1, std::min(claims, threads));
	return std::move(state);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &, TableFunctionInitInput &,
                                              GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

// Fills a chunk with (node, ancestor) rows, a node's ancestors bottom up, claiming more
// nodes as its claim runs out; threads claim independently, so the rows come out in no
// particular order.
void Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &state = data_p.global_state->Cast<GlobalState>();
	auto &local = data_p.local_state->Cast<LocalState>();
	const auto &index = state.tree.index;
	const auto &ids = state.tree.node_ids;

	auto node_out = FlatVector::GetDataMutable<int64_t>(output.data[0]);
	auto ancestor_out = FlatVector::GetDataMutable<int64_t>(output.data[1]);
	auto edges_out = FlatVector::GetDataMutable<int64_t>(output.data[2]);
	auto distance_out = FlatVector::GetDataMutable<double>(output.data[3]);
	auto &distance_valid = FlatVector::ValidityMutable(output.data[3]);
	idx_t count = 0;
	while (count < STANDARD_VECTOR_SIZE) {
		if (local.ancestor == miint::TreeIndex::NO_PARENT) {
			// The last node is done (past the root), or none has started.
			if (local.next == local.end) {
				const idx_t start = state.next_claim.fetch_add(kNodesPerClaim);
				if (start >= state.node_count) {
					break;
				}
				local.next = start;
				local.end = std::min(start + kNodesPerClaim, state.node_count);
			}
			local.node = state.NodeAt(local.next++);
			local.ancestor = local.node;
			local.edges = 0;
		}
		node_out[count] = ids[local.node];
		ancestor_out[count] = ids[local.ancestor];
		edges_out[count] = local.edges;
		const auto d = index.distance_to_ancestor(local.node, local.ancestor);
		if (d) {
			distance_out[count] = *d;
			distance_valid.SetValid(count);
		} else {
			distance_valid.SetInvalid(count);
		}
		++count;
		local.ancestor = index.parent(local.ancestor);
		++local.edges;
	}
	output.SetChildCardinality(count);
}

} // namespace

TableFunction TreeAncestorsTableFunction::GetFunction() {
	TableFunction tf("tree_ancestors", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	AddNamedParameter(tf, "nodes_table", LogicalType::VARCHAR);
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void TreeAncestorsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
