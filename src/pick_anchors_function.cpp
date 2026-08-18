#include "pick_anchors.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {
namespace {

using unifrac_internal::ReadCoordinateTable;
using unifrac_internal::ResolveSampleIdOutputType;

// Selection rules. Defaults to `stratified` because that is the measured winner
// for progressive-PCoA anchors; see pick_anchors.hpp for the bake-off and the
// mechanism.
enum class SelectionMethod { STRATIFIED, FARTHEST_POINT };

struct PickAnchorsData : public TableFunctionData {
	std::string table_name;
	int32_t n_anchors = 0;
	SelectionMethod method = SelectionMethod::STRATIFIED;
	int64_t seed = 0;
	int32_t n_dims = 3;
	int32_t n_bins = 4;
	LogicalType sample_id_type = LogicalType::VARCHAR;
};

struct PickAnchorsGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> anchors; // in selection order
	idx_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<FunctionData> PickAnchorsBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<PickAnchorsData>();
	data->table_name = input.inputs[0].GetValue<string>();
	if (data->table_name.empty()) {
		throw BinderException("pick_anchors: coordinate-table name must not be empty");
	}

	bool has_n = false;
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "n_anchors") {
			data->n_anchors = kv.second.GetValue<int32_t>();
			has_n = true;
		} else if (key == "method") {
			const auto method = StringUtil::Lower(kv.second.GetValue<string>());
			if (method == "stratified") {
				data->method = SelectionMethod::STRATIFIED;
			} else if (method == "farthest_point") {
				data->method = SelectionMethod::FARTHEST_POINT;
			} else {
				throw BinderException("pick_anchors: unknown method '%s'; expected 'stratified' or 'farthest_point'",
				                      kv.second.GetValue<string>());
			}
		} else if (key == "seed") {
			data->seed = kv.second.GetValue<int64_t>();
		} else if (key == "n_dims") {
			data->n_dims = kv.second.GetValue<int32_t>();
		} else if (key == "n_bins") {
			data->n_bins = kv.second.GetValue<int32_t>();
		}
	}
	if (!has_n) {
		throw BinderException("pick_anchors: n_anchors is required (there is no default anchor count); "
		                      "pass n_anchors := N");
	}
	if (data->n_anchors < 1) {
		throw BinderException("pick_anchors: n_anchors must be >= 1 (got %d)", data->n_anchors);
	}
	if (data->n_dims < 0) {
		throw BinderException("pick_anchors: n_dims must be >= 0 (0 = every axis; got %d)", data->n_dims);
	}
	if (data->n_bins < 1) {
		throw BinderException("pick_anchors: n_bins must be >= 1 (got %d)", data->n_bins);
	}

	// Mirror the coordinate table's sample_id type onto the output (same as
	// cluster_kmeans): a metadata lookup at bind, so the table itself is read once,
	// in InitGlobal.
	auto cols = GetTableOrViewColumns(context, data->table_name, "coordinate-table");
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == "sample_id") {
			data->sample_id_type = ResolveSampleIdOutputType(cols.types[i]);
			break;
		}
	}

	names.emplace_back("anchor_rank");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("sample_id");
	return_types.emplace_back(data->sample_id_type);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> PickAnchorsInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<PickAnchorsData>();
	auto gstate = make_uniq<PickAnchorsGlobalState>();
	gstate->sample_id_type = data.sample_id_type;

	auto coords = ReadCoordinateTable(context, data.table_name, "pick_anchors", data.n_dims);
	if (static_cast<uint32_t>(data.n_anchors) > coords.n_samples) {
		throw InvalidInputException(
		    "pick_anchors: n_anchors (%d) exceeds the %u distinct sample(s) in the coordinate-table", data.n_anchors,
		    coords.n_samples);
	}

	std::vector<uint32_t> picked;
	try {
		const auto k = static_cast<uint32_t>(data.n_anchors);
		if (data.method == SelectionMethod::STRATIFIED) {
			picked = miint::SelectStratified(coords.points, coords.sample_ids, coords.n_samples, coords.n_dims, k,
			                                 static_cast<uint32_t>(data.n_bins), data.seed);
		} else {
			picked = miint::SelectFarthestPoint(coords.points, coords.n_samples, coords.n_dims, k);
		}
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	gstate->anchors.reserve(picked.size());
	for (const auto index : picked) {
		gstate->anchors.push_back(coords.sample_ids[index]);
	}
	return std::move(gstate);
}

void PickAnchorsExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<PickAnchorsGlobalState>();
	const idx_t total = gstate.anchors.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t n = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto rank_data = FlatVector::GetData<int32_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	for (idx_t i = 0; i < n; ++i) {
		rank_data[i] = static_cast<int32_t>(gstate.cursor + i);
		EmitIdCell(sample_id_vec, i, gstate.anchors[gstate.cursor + i], gstate.sample_id_type);
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterPickAnchors(ExtensionLoader &loader) {
	TableFunction fn("pick_anchors", {LogicalType::VARCHAR}, PickAnchorsExecute, PickAnchorsBind,
	                 PickAnchorsInitGlobal);
	fn.named_parameters["n_anchors"] = LogicalType::INTEGER;
	fn.named_parameters["method"] = LogicalType::VARCHAR;
	fn.named_parameters["seed"] = LogicalType::BIGINT;
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["n_bins"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
