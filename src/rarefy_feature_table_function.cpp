#include "unifrac_table_functions.hpp"

#include <algorithm>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"
#include "api.hpp"
#include "unifrac_omp_scope.hpp"
#include "unifrac_subsample_bridge.hpp"
#include "unifrac_support_biom.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {
namespace {

using unifrac_internal::ReadFeatureTable;
using unifrac_internal::ResolveSampleIdOutputType;
using unifrac_internal::ResolveThreadsParameter;

struct RarefyData : public TableFunctionData {
	std::vector<miint::unifrac::CooRow> rows;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
};

struct RarefyGlobalState : public GlobalTableFunctionState {
	std::vector<miint::unifrac::CooRow> rows;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Flatten a rarefied biom view back into COO rows.
//
// The view is CSR obs-major (api.hpp: "CSR-encoded input table"): indptr has
// n_obs + 1 entries, and obs i owns indices/data in [indptr[i], indptr[i+1]),
// where `indices` holds SAMPLE column indices. Misreading this as CSC would
// transpose the table silently, which is why the layout is spelled out here and
// pinned by the depth-oracle cases in rarefy_feature_table.test.
//
// BridgeSubsample has already applied the sparse-storage invariant (cells that
// subsample to zero are gone, and so are observations left empty across every
// surviving sample), so no filtering is needed here.
std::vector<miint::unifrac::CooRow> FlattenToCoo(const miint::unifrac::UnifracSupportBiomView &view) {
	const support_biom_t *b = view.support_biom();
	std::vector<miint::unifrac::CooRow> out;
	const auto n_obs = static_cast<uint32_t>(b->n_obs);
	out.reserve(n_obs > 0 ? b->indptr[n_obs] : 0);
	for (uint32_t i = 0; i < n_obs; ++i) {
		for (uint32_t j = b->indptr[i]; j < b->indptr[i + 1]; ++j) {
			out.push_back({std::string(b->sample_ids[b->indices[j]]), std::string(b->obs_ids[i]), b->data[j]});
		}
	}
	return out;
}

unique_ptr<FunctionData> RarefyBind(ClientContext &context, TableFunctionBindInput &input,
                                    vector<LogicalType> &return_types, vector<string> &names) {
	const std::string table_name = input.inputs[0].GetValue<string>();
	if (table_name.empty()) {
		throw BinderException("rarefy_feature_table: feature-table name must not be empty");
	}

	// INTEGER-typed named params match the project convention (int32_t).
	// `depth` is deliberately required and has no default: there is no defensible
	// universal per-sample depth, and silently picking one would change every
	// downstream diversity result.
	bool has_depth = false;
	int32_t depth = 0;
	bool with_replacement = false;
	int32_t seed = -1;
	int32_t threads = 0; // 0 = follow DuckDB's TaskScheduler::NumberOfThreads()
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "depth") {
			depth = kv.second.GetValue<int32_t>();
			has_depth = true;
		} else if (key == "with_replacement") {
			with_replacement = kv.second.GetValue<bool>();
		} else if (key == "seed") {
			seed = kv.second.GetValue<int32_t>();
		} else if (key == "threads") {
			threads = kv.second.GetValue<int32_t>();
		}
	}
	if (!has_depth) {
		throw BinderException("rarefy_feature_table: depth is required (there is no default per-sample depth); "
		                      "pass depth := N");
	}
	if (depth < 1) {
		throw BinderException("rarefy_feature_table: depth must be >= 1 (got %d)", depth);
	}
	const int n_threads = ResolveThreadsParameter(context, threads, "rarefy_feature_table");

	LogicalType sample_id_col_type = LogicalType::VARCHAR;
	LogicalType feature_id_col_type = LogicalType::VARCHAR;
	auto coo_rows =
	    ReadFeatureTable(context, table_name, "rarefy_feature_table", &sample_id_col_type, &feature_id_col_type);
	if (coo_rows.empty()) {
		throw InvalidInputException("rarefy_feature_table: feature-table '%s' is empty after dropping NULL/zero rows",
		                            table_name);
	}
	miint::unifrac::UnifracSupportBiomView biom_view = [&]() {
		try {
			return miint::unifrac::UnifracSupportBiomView::FromCoo(std::move(coo_rows));
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("rarefy_feature_table: %s", e.what());
		}
	}();

	auto data = make_uniq<RarefyData>();
	data->sample_id_type = ResolveSampleIdOutputType(sample_id_col_type);
	data->feature_id_type = ResolveSampleIdOutputType(feature_id_col_type);

	// The draw itself. ComputeCallScope pins this thread's OpenMP width and
	// guarantees a non-negative seed, so the draw uses a per-call generator rather
	// than libssu's process-global one — required for concurrent callers, and the
	// reason an unseeded call here is still safe (it gets a derived per-call seed).
	// NOTE the width matters to the RESULT, not just the speed: libssu distributes
	// the draw across the OpenMP team, so a seeded rarefaction reproduces per
	// thread count, not across thread counts (docs/diversity.md).
	miint::unifrac::UnifracSupportBiomView rarefied = [&]() {
		try {
			miint::unifrac::ComputeCallScope scope(n_threads, seed);
			return miint::unifrac::BridgeSubsample(biom_view, static_cast<uint32_t>(depth), with_replacement,
			                                       scope.seed());
		} catch (const std::invalid_argument &e) {
			// Every sample fell below `depth`. Fail loudly: an empty result would
			// read as a successful analysis of nothing.
			throw InvalidInputException("rarefy_feature_table: no sample in '%s' has at least depth := %d total "
			                            "counts, so every sample was dropped (%s)",
			                            table_name, depth, e.what());
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("rarefy_feature_table: %s", e.what());
		}
	}();

	data->rows = FlattenToCoo(rarefied);

	names.emplace_back("sample_id");
	return_types.emplace_back(data->sample_id_type);
	names.emplace_back("feature_id");
	return_types.emplace_back(data->feature_id_type);
	names.emplace_back("value");
	return_types.emplace_back(LogicalType::DOUBLE);

	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> RarefyInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<RarefyData>();
	auto gstate = make_uniq<RarefyGlobalState>();
	gstate->rows = std::move(data.rows);
	gstate->sample_id_type = data.sample_id_type;
	gstate->feature_id_type = data.feature_id_type;
	return std::move(gstate);
}

void RarefyExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<RarefyGlobalState>();
	const idx_t total = gstate.rows.size();
	if (gstate.cursor >= total) {
		output.SetCardinality(0);
		return;
	}
	const idx_t n = std::min<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto &sample_id_vec = output.data[0];
	auto &feature_id_vec = output.data[1];
	auto value_data = FlatVector::GetData<double>(output.data[2]);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		// EmitIdCell mirrors the id type; its ""/"*"→NULL sentinel branch is
		// unreachable here (ReadFeatureTable drops NULL sample/feature ids).
		EmitIdCell(sample_id_vec, i, r.sample_id, gstate.sample_id_type);
		EmitIdCell(feature_id_vec, i, r.feature_id, gstate.feature_id_type);
		value_data[i] = r.count;
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterRarefyFeatureTable(ExtensionLoader &loader) {
	TableFunction fn("rarefy_feature_table", {LogicalType::VARCHAR}, RarefyExecute, RarefyBind, RarefyInitGlobal);
	fn.named_parameters["depth"] = LogicalType::INTEGER;
	fn.named_parameters["with_replacement"] = LogicalType::BOOLEAN;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	fn.named_parameters["threads"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
