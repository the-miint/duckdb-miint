#include "cluster_upgma.hpp"

#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <stdexcept>

namespace duckdb {

namespace {

using unifrac_internal::ReadDistanceTable;

struct ClusterUpgmaBindData : public TableFunctionData {
	std::string table_name;
};

// One emitted tree-table row (read_newick shape).
struct UpgmaRow {
	int64_t node_index;
	std::string name; // sample id for tips, "" for internal nodes (matches read_newick)
	double branch_length;
	int64_t edge_id;
	int64_t parent_index;
	bool has_parent; // false only for the root
	bool is_tip;
};

struct ClusterUpgmaGlobalState : public GlobalTableFunctionState {
	std::vector<UpgmaRow> rows;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> ClusterUpgmaBind(ClientContext &, TableFunctionBindInput &input,
                                          vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ClusterUpgmaBindData>();
	data->table_name = input.inputs[0].GetValue<string>();

	// Emit the read_newick tree-table schema verbatim, so cluster_upgma output is
	// a valid tree table for every downstream tree consumer.
	names = {"node_index", "name", "branch_length", "edge_id", "parent_index", "is_tip"};
	return_types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::DOUBLE,
	                LogicalType::BIGINT, LogicalType::BIGINT,  LogicalType::BOOLEAN};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ClusterUpgmaInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<ClusterUpgmaBindData>();
	auto gstate = make_uniq<ClusterUpgmaGlobalState>();

	// ReadDistanceTable validates the (sample_a, sample_b, distance) schema,
	// completeness, n >= 2, and non-negative/finite distances (throwing a
	// prefixed InvalidInputException), then hands back a dense fp32 matrix with
	// lexicographically-sorted sample ids.
	auto dist = ReadDistanceTable(context, data.table_name, "cluster_upgma");
	const uint32_t n = dist.n_samples;

	std::vector<double> dmat(static_cast<size_t>(n) * n);
	for (size_t i = 0; i < dmat.size(); ++i) {
		dmat[i] = static_cast<double>(dist.matrix[i]);
	}

	std::vector<miint::UpgmaNode> nodes;
	try {
		nodes = miint::UpgmaAverageLinkage(dmat, n);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	gstate->rows.reserve(nodes.size());
	for (size_t i = 0; i < nodes.size(); ++i) {
		const auto &nd = nodes[i];
		UpgmaRow row;
		row.node_index = static_cast<int64_t>(i);
		row.name = nd.leaf_index >= 0 ? dist.sample_ids[static_cast<size_t>(nd.leaf_index)] : std::string();
		row.branch_length = nd.branch_length;
		row.edge_id = static_cast<int64_t>(i);
		row.has_parent = nd.parent >= 0;
		row.parent_index = nd.parent >= 0 ? static_cast<int64_t>(nd.parent) : 0;
		row.is_tip = nd.leaf_index >= 0;
		gstate->rows.push_back(std::move(row));
	}
	return std::move(gstate);
}

void ClusterUpgmaExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<ClusterUpgmaGlobalState>();
	const idx_t total = g.rows.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto node_index = FlatVector::GetData<int64_t>(output.data[0]);
	auto &name_vec = output.data[1];
	auto name_data = FlatVector::GetData<string_t>(name_vec);
	auto branch_length = FlatVector::GetData<double>(output.data[2]);
	auto edge_id = FlatVector::GetData<int64_t>(output.data[3]);
	auto parent_index = FlatVector::GetData<int64_t>(output.data[4]);
	auto &parent_validity = FlatVector::Validity(output.data[4]);
	auto is_tip = FlatVector::GetData<bool>(output.data[5]);

	for (idx_t r = 0; r < count; ++r) {
		const auto &row = g.rows[g.cursor + r];
		node_index[r] = row.node_index;
		name_data[r] = StringVector::AddString(name_vec, row.name);
		branch_length[r] = row.branch_length;
		edge_id[r] = row.edge_id;
		if (row.has_parent) {
			parent_index[r] = row.parent_index;
			parent_validity.SetValid(r);
		} else {
			parent_index[r] = 0;
			parent_validity.SetInvalid(r); // root has no parent
		}
		is_tip[r] = row.is_tip;
	}
	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterClusterUpgma(ExtensionLoader &loader) {
	TableFunction fn("cluster_upgma", {LogicalType::VARCHAR}, ClusterUpgmaExecute, ClusterUpgmaBind,
	                 ClusterUpgmaInitGlobal);
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
