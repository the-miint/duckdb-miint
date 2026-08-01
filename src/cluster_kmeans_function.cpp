#include "cluster_kmeans.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/query_result.hpp"

#include <algorithm>
#include <map>
#include <stdexcept>

namespace duckdb {

namespace {

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

	auto conn = MakeReadOnlyHelperConnection(context);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(data.table_name);
	auto probe = conn.Query("SELECT sample_id::VARCHAR, axis::INTEGER, coordinate::DOUBLE FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException(
		    "cluster_kmeans: coordinate-table '%s' must expose (sample_id, axis INTEGER, coordinate DOUBLE): %s",
		    data.table_name, probe->GetError());
	}
	auto result = conn.Query("SELECT sample_id::VARCHAR, axis::INTEGER, coordinate::DOUBLE FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("cluster_kmeans: failed to read coordinate-table '%s': %s", data.table_name,
		                            result->GetError());
	}

	// Collect (sample_id, axis) -> coordinate, dropping NULLs. std::map keeps the
	// axis set sorted ascending; sample ids sorted lexicographically below.
	std::map<std::string, std::map<int32_t, double>> bysample;
	std::map<int32_t, bool> axis_seen;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat s_u, a_u, c_u;
		chunk->data[0].ToUnifiedFormat(rn, s_u);
		chunk->data[1].ToUnifiedFormat(rn, a_u);
		chunk->data[2].ToUnifiedFormat(rn, c_u);
		auto s_d = UnifiedVectorFormat::GetData<string_t>(s_u);
		auto a_d = UnifiedVectorFormat::GetData<int32_t>(a_u);
		auto c_d = UnifiedVectorFormat::GetData<double>(c_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto si = s_u.sel->get_index(i);
			const auto ai = a_u.sel->get_index(i);
			const auto ci = c_u.sel->get_index(i);
			if (!s_u.validity.RowIsValid(si) || !a_u.validity.RowIsValid(ai) || !c_u.validity.RowIsValid(ci)) {
				continue;
			}
			const std::string sid = s_d[si].GetString();
			const int32_t axis = a_d[ai];
			auto &am = bysample[sid];
			if (!am.emplace(axis, c_d[ci]).second) {
				throw InvalidInputException(
				    "cluster_kmeans: coordinate-table '%s' has duplicate (sample_id='%s', axis=%d); "
				    "pass a single-iteration coordinate table",
				    data.table_name, sid, axis);
			}
			axis_seen[axis] = true;
		}
	}

	const uint32_t n = static_cast<uint32_t>(bysample.size());
	if (n < static_cast<uint32_t>(data.k)) {
		throw InvalidInputException("cluster_kmeans: k (%d) must be <= number of samples (%u)", data.k, n);
	}

	// Axis order: ascending; optionally capped to the first n_dims axes.
	std::vector<int32_t> axes;
	axes.reserve(axis_seen.size());
	for (auto &kv : axis_seen) {
		axes.push_back(kv.first);
	}
	if (data.n_dims > 0 && static_cast<size_t>(data.n_dims) < axes.size()) {
		axes.resize(static_cast<size_t>(data.n_dims));
	}
	const uint32_t d = static_cast<uint32_t>(axes.size());
	if (d == 0) {
		throw InvalidInputException("cluster_kmeans: coordinate-table '%s' has no coordinates", data.table_name);
	}

	// Dense n x d point matrix in sorted-sample-id order; every used axis must be
	// present for every sample.
	std::vector<double> points(static_cast<size_t>(n) * d);
	std::vector<std::string> sample_ids;
	sample_ids.reserve(n);
	uint32_t si = 0;
	for (auto &sp : bysample) {
		sample_ids.push_back(sp.first);
		for (uint32_t j = 0; j < d; ++j) {
			auto it = sp.second.find(axes[j]);
			if (it == sp.second.end()) {
				throw InvalidInputException("cluster_kmeans: sample '%s' is missing a coordinate for axis %d", sp.first,
				                            axes[j]);
			}
			points[static_cast<size_t>(si) * d + j] = it->second;
		}
		++si;
	}

	miint::KMeansResult res;
	try {
		res = miint::KMeans(points, n, d, data.k, data.seed, data.max_iter, data.n_init);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	gstate->sample_ids = std::move(sample_ids);
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
