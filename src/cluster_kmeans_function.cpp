#include "cluster_kmeans.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <algorithm>
#include <stdexcept>

namespace duckdb {

namespace {

using unifrac_internal::ReadCoordinateTable;
using unifrac_internal::ResolveSampleIdOutputType;

struct ClusterKmeansBindData : public TableFunctionData {
	std::string table_name;
	int32_t k = 0;
	int64_t seed = 0;
	int32_t max_iter = 100;
	int32_t n_init = 10;
	int32_t n_dims = 0; // 0 => use every axis present
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

struct ClusterKmeansGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sample_ids; // sorted, aligned to assignments
	std::vector<int32_t> assignments;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> ClusterKmeansBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<ClusterKmeansBindData>();
	data->table_name = input.inputs[0].GetValue<string>();

	bool has_k = false;
	for (auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "k") {
			data->k = kv.second.GetValue<int32_t>();
			has_k = true;
		} else if (key == "seed") {
			data->seed = kv.second.GetValue<int64_t>();
		} else if (key == "max_iter") {
			data->max_iter = kv.second.GetValue<int32_t>();
		} else if (key == "n_init") {
			data->n_init = kv.second.GetValue<int32_t>();
		} else if (key == "n_dims") {
			data->n_dims = kv.second.GetValue<int32_t>();
		}
	}
	if (!has_k) {
		throw BinderException("cluster_kmeans: the number of clusters 'k' is required (e.g. k := 3)");
	}
	if (data->k < 1) {
		throw BinderException("cluster_kmeans: k must be >= 1 (got %d)", data->k);
	}
	if (data->max_iter < 1) {
		throw BinderException("cluster_kmeans: max_iter must be >= 1 (got %d)", data->max_iter);
	}
	if (data->n_init < 1) {
		throw BinderException("cluster_kmeans: n_init must be >= 1 (got %d)", data->n_init);
	}
	if (data->n_dims < 0) {
		throw BinderException("cluster_kmeans: n_dims must be >= 0 (0 = all axes; got %d)", data->n_dims);
	}

	// Mirror the coordinate table's sample_id type onto the output.
	auto cols = GetTableOrViewColumns(context, data->table_name, "coordinate-table");
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == "sample_id") {
			data->sample_id_type = ResolveSampleIdOutputType(cols.types[i]);
			break;
		}
	}

	names = {"sample_id", "cluster"};
	return_types = {data->sample_id_type, LogicalType::INTEGER};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ClusterKmeansInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<ClusterKmeansBindData>();
	auto gstate = make_uniq<ClusterKmeansGlobalState>();
	gstate->sample_id_type = data.sample_id_type;

	// ReadCoordinateTable validates the (sample_id, axis, coordinate) schema, caps
	// the axes to the leading n_dims, rejects duplicate (sample_id, axis) pairs and
	// ragged samples, and hands back a dense row-major cloud in
	// lexicographically-sorted sample-id order.
	auto coords = ReadCoordinateTable(context, data.table_name, "cluster_kmeans", data.n_dims);
	if (coords.n_samples < static_cast<uint32_t>(data.k)) {
		throw InvalidInputException("cluster_kmeans: k (%d) must be <= number of samples (%u)", data.k,
		                            coords.n_samples);
	}

	miint::KMeansResult res;
	try {
		res = miint::KMeans(coords.points, coords.n_samples, coords.n_dims, data.k, data.seed, data.max_iter,
		                    data.n_init);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	gstate->sample_ids = std::move(coords.sample_ids);
	gstate->assignments = std::move(res.assignments);
	return std::move(gstate);
}

void ClusterKmeansExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<ClusterKmeansGlobalState>();
	const idx_t total = g.sample_ids.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &sid_vec = output.data[0];
	auto cluster = FlatVector::GetData<int32_t>(output.data[1]);
	for (idx_t r = 0; r < count; ++r) {
		const idx_t k = g.cursor + r;
		EmitIdCell(sid_vec, r, g.sample_ids[k], g.sample_id_type);
		cluster[r] = g.assignments[k];
	}
	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterClusterKmeans(ExtensionLoader &loader) {
	TableFunction fn("cluster_kmeans", {LogicalType::VARCHAR}, ClusterKmeansExecute, ClusterKmeansBind,
	                 ClusterKmeansInitGlobal);
	fn.named_parameters["k"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::BIGINT;
	fn.named_parameters["max_iter"] = LogicalType::INTEGER;
	fn.named_parameters["n_init"] = LogicalType::INTEGER;
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
