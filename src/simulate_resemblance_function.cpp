#include "simulate_resemblance.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <cmath>
#include <stdexcept>

namespace duckdb {

namespace {

//! Parse + validate the positional `abundances` LIST(DOUBLE) argument shared by
//! both simulators: non-empty, no NULLs, every value finite and >= 0, not all zero.
std::vector<double> ParseAbundancesArg(const Value &list_val, const char *fn) {
	auto children = ListValue::GetChildren(list_val);
	if (children.empty()) {
		throw InvalidInputException("%s: abundances must be a non-empty list", fn);
	}
	std::vector<double> abundances;
	abundances.reserve(children.size());
	double abund_sum = 0.0;
	for (idx_t i = 0; i < children.size(); i++) {
		if (children[i].IsNull()) {
			throw InvalidInputException("%s: abundances must not contain NULL (index %d)", fn, static_cast<int>(i));
		}
		double v = children[i].GetValue<double>();
		if (!(v >= 0.0) || !std::isfinite(v)) {
			throw InvalidInputException("%s: abundances must be finite and non-negative (got %g at index %d)", fn, v,
			                            static_cast<int>(i));
		}
		abundances.push_back(v);
		abund_sum += v;
	}
	if (abund_sum <= 0.0) {
		throw InvalidInputException("%s: abundances must not be all zero", fn);
	}
	return abundances;
}

struct SimGradientBindData : public TableFunctionData {
	std::vector<double> abundances;
	int32_t num_samples = 0;
	int64_t seqs_per_sample = 0;
	double sp_width = 0.1;
	double noise = 0.0;
	std::string noise_type = "+species";
	double range_lo = 0.1;
	double range_hi = 0.9;
	int64_t seed = -1;
};

struct SimGlobalState : public GlobalTableFunctionState {
	miint::SimulationCOO coo;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1; // single-threaded generation + emission; cursor is not atomic
	}
};

unique_ptr<FunctionData> SimGradientBind(ClientContext &context, TableFunctionBindInput &input,
                                         vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<SimGradientBindData>();

	// Positional: abundances LIST(DOUBLE), num_samples INTEGER, seqs_per_sample BIGINT.
	data->abundances = ParseAbundancesArg(input.inputs[0], "simulate_gradient_otus");

	data->num_samples = input.inputs[1].GetValue<int32_t>();
	if (data->num_samples < 1) {
		throw InvalidInputException("simulate_gradient_otus: num_samples must be >= 1 (got %d)", data->num_samples);
	}
	data->seqs_per_sample = input.inputs[2].GetValue<int64_t>();
	if (data->seqs_per_sample < 1) {
		throw InvalidInputException("simulate_gradient_otus: seqs_per_sample must be >= 1 (got %lld)",
		                            static_cast<long long>(data->seqs_per_sample));
	}

	for (auto &kv : input.named_parameters) {
		auto key = StringUtil::Lower(kv.first);
		if (key == "sp_width") {
			data->sp_width = kv.second.GetValue<double>();
		} else if (key == "noise") {
			data->noise = kv.second.GetValue<double>();
		} else if (key == "noise_type") {
			data->noise_type = kv.second.GetValue<string>();
		} else if (key == "range_lo") {
			data->range_lo = kv.second.GetValue<double>();
		} else if (key == "range_hi") {
			data->range_hi = kv.second.GetValue<double>();
		} else if (key == "seed") {
			data->seed = kv.second.GetValue<int64_t>();
		}
	}

	if (!(data->sp_width > 0.0) || !std::isfinite(data->sp_width)) {
		throw InvalidInputException("simulate_gradient_otus: sp_width must be a finite positive number (got %g)",
		                            data->sp_width);
	}
	if (!(data->noise >= 0.0) || !std::isfinite(data->noise)) {
		throw InvalidInputException("simulate_gradient_otus: noise must be a finite non-negative number (got %g)",
		                            data->noise);
	}
	if (data->noise != 0.0 && data->noise_type != "+species" && data->noise_type != "*sample" &&
	    data->noise_type != "+sample") {
		throw InvalidInputException(
		    "simulate_gradient_otus: noise_type must be '+species', '*sample' or '+sample' (got '%s')",
		    data->noise_type);
	}
	// isfinite BEFORE the ordering check: NaN silently satisfies `range_lo > range_hi`
	// (NaN comparisons are false) and would then poison every position/weight.
	if (!std::isfinite(data->range_lo)) {
		throw InvalidInputException("simulate_gradient_otus: range_lo must be finite (got %g)", data->range_lo);
	}
	if (!std::isfinite(data->range_hi)) {
		throw InvalidInputException("simulate_gradient_otus: range_hi must be finite (got %g)", data->range_hi);
	}
	if (data->range_lo > data->range_hi) {
		throw InvalidInputException("simulate_gradient_otus: range_lo (%g) must be <= range_hi (%g)", data->range_lo,
		                            data->range_hi);
	}

	names = {"sample_id", "otu_id", "count", "gradient_position"};
	return_types = {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::BIGINT, LogicalType::DOUBLE};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> SimGradientInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<SimGradientBindData>();
	auto gstate = make_uniq<SimGlobalState>();
	// The pure core throws std::invalid_argument for data-dependent degeneracies
	// (e.g. a sample whose abundances all zero out under noise). Surface it as the
	// project's InvalidInputException (its message already carries the fn prefix).
	try {
		gstate->coo = miint::SimulateGradient(data.abundances, data.num_samples, data.seqs_per_sample, data.sp_width,
		                                      data.noise, data.noise_type, data.range_lo, data.range_hi, data.seed);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
	return std::move(gstate);
}

void SimGradientExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<SimGlobalState>();
	const auto &coo = gstate.coo;
	const idx_t total = coo.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto sample_id = FlatVector::GetData<int32_t>(output.data[0]);
	auto otu_id = FlatVector::GetData<int32_t>(output.data[1]);
	auto counts = FlatVector::GetData<int64_t>(output.data[2]);
	auto position = FlatVector::GetData<double>(output.data[3]);

	for (idx_t i = 0; i < count; i++) {
		const idx_t j = gstate.cursor + i;
		sample_id[i] = coo.sample_id[j];
		otu_id[i] = coo.otu_id[j];
		counts[i] = coo.count[j];
		position[i] = coo.ground_truth[j];
	}
	gstate.cursor += count;
	output.SetCardinality(count);
}

unique_ptr<LocalTableFunctionState> SimInitLocal(ExecutionContext &, TableFunctionInitInput &,
                                                 GlobalTableFunctionState *) {
	return nullptr;
}

// ---------------------------------------------------------------------------
// simulate_cluster_otus
// ---------------------------------------------------------------------------

struct SimClusterBindData : public TableFunctionData {
	std::vector<double> abundances;
	int64_t seqs_per_sample = 0;
	std::vector<int32_t> cluster_sizes = {30, 30, 30};
	double cluster_spacing = 1.0;
	double sample_spacing = 0.5;
	std::string noise_type = "*sample";
	std::string normalization = "clip";
	int64_t seed = -1;
};

unique_ptr<FunctionData> SimClusterBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<SimClusterBindData>();

	// Positional: abundances LIST(DOUBLE), seqs_per_sample BIGINT.
	data->abundances = ParseAbundancesArg(input.inputs[0], "simulate_cluster_otus");

	data->seqs_per_sample = input.inputs[1].GetValue<int64_t>();
	if (data->seqs_per_sample < 1) {
		throw InvalidInputException("simulate_cluster_otus: seqs_per_sample must be >= 1 (got %lld)",
		                            static_cast<long long>(data->seqs_per_sample));
	}

	for (auto &kv : input.named_parameters) {
		auto key = StringUtil::Lower(kv.first);
		if (key == "cluster_sizes") {
			auto sizes = ListValue::GetChildren(kv.second);
			if (sizes.empty()) {
				throw InvalidInputException("simulate_cluster_otus: cluster_sizes must be a non-empty list");
			}
			data->cluster_sizes.clear();
			data->cluster_sizes.reserve(sizes.size());
			for (idx_t i = 0; i < sizes.size(); i++) {
				if (sizes[i].IsNull()) {
					throw InvalidInputException("simulate_cluster_otus: cluster_sizes must not contain NULL");
				}
				int32_t sz = sizes[i].GetValue<int32_t>();
				if (sz < 1) {
					throw InvalidInputException("simulate_cluster_otus: cluster_sizes must all be >= 1 (got %d)", sz);
				}
				data->cluster_sizes.push_back(sz);
			}
		} else if (key == "cluster_spacing") {
			data->cluster_spacing = kv.second.GetValue<double>();
		} else if (key == "sample_spacing") {
			data->sample_spacing = kv.second.GetValue<double>();
		} else if (key == "noise_type") {
			data->noise_type = kv.second.GetValue<string>();
		} else if (key == "normalization") {
			data->normalization = kv.second.GetValue<string>();
		} else if (key == "seed") {
			data->seed = kv.second.GetValue<int64_t>();
		}
	}

	if (!(data->cluster_spacing >= 0.0) || !std::isfinite(data->cluster_spacing)) {
		throw InvalidInputException("simulate_cluster_otus: cluster_spacing must be a finite non-negative number "
		                            "(got %g)",
		                            data->cluster_spacing);
	}
	if (!(data->sample_spacing >= 0.0) || !std::isfinite(data->sample_spacing)) {
		throw InvalidInputException("simulate_cluster_otus: sample_spacing must be a finite non-negative number "
		                            "(got %g)",
		                            data->sample_spacing);
	}
	if (data->noise_type != "*sample" && data->noise_type != "+sample") {
		throw InvalidInputException("simulate_cluster_otus: noise_type must be '*sample' or '+sample' (got '%s')",
		                            data->noise_type);
	}
	if (data->normalization != "clip" && data->normalization != "rescale") {
		throw InvalidInputException("simulate_cluster_otus: normalization must be 'clip' or 'rescale' (got '%s')",
		                            data->normalization);
	}

	names = {"sample_id", "otu_id", "count", "cluster_id"};
	return_types = {LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::BIGINT, LogicalType::INTEGER};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> SimClusterInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<SimClusterBindData>();
	auto gstate = make_uniq<SimGlobalState>();
	// See SimGradientInitGlobal: surface core degeneracies as InvalidInputException.
	try {
		gstate->coo =
		    miint::SimulateCluster(data.abundances, data.cluster_sizes, data.seqs_per_sample, data.cluster_spacing,
		                           data.sample_spacing, data.noise_type, data.normalization, data.seed);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
	return std::move(gstate);
}

void SimClusterExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<SimGlobalState>();
	const auto &coo = gstate.coo;
	const idx_t total = coo.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto sample_id = FlatVector::GetData<int32_t>(output.data[0]);
	auto otu_id = FlatVector::GetData<int32_t>(output.data[1]);
	auto counts = FlatVector::GetData<int64_t>(output.data[2]);
	auto cluster_id = FlatVector::GetData<int32_t>(output.data[3]);

	for (idx_t i = 0; i < count; i++) {
		const idx_t j = gstate.cursor + i;
		sample_id[i] = coo.sample_id[j];
		otu_id[i] = coo.otu_id[j];
		counts[i] = coo.count[j];
		cluster_id[i] = static_cast<int32_t>(coo.ground_truth[j]);
	}
	gstate.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterSimulateResemblance(ExtensionLoader &loader) {
	TableFunction gradient("simulate_gradient_otus",
	                       {LogicalType::LIST(LogicalType::DOUBLE), LogicalType::INTEGER, LogicalType::BIGINT},
	                       SimGradientExecute, SimGradientBind, SimGradientInitGlobal, SimInitLocal);
	gradient.named_parameters["sp_width"] = LogicalType::DOUBLE;
	gradient.named_parameters["noise"] = LogicalType::DOUBLE;
	gradient.named_parameters["noise_type"] = LogicalType::VARCHAR;
	gradient.named_parameters["range_lo"] = LogicalType::DOUBLE;
	gradient.named_parameters["range_hi"] = LogicalType::DOUBLE;
	gradient.named_parameters["seed"] = LogicalType::BIGINT;
	gradient.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(gradient);

	TableFunction cluster("simulate_cluster_otus", {LogicalType::LIST(LogicalType::DOUBLE), LogicalType::BIGINT},
	                      SimClusterExecute, SimClusterBind, SimClusterInitGlobal, SimInitLocal);
	cluster.named_parameters["cluster_sizes"] = LogicalType::LIST(LogicalType::INTEGER);
	cluster.named_parameters["cluster_spacing"] = LogicalType::DOUBLE;
	cluster.named_parameters["sample_spacing"] = LogicalType::DOUBLE;
	cluster.named_parameters["noise_type"] = LogicalType::VARCHAR;
	cluster.named_parameters["normalization"] = LogicalType::VARCHAR;
	cluster.named_parameters["seed"] = LogicalType::BIGINT;
	cluster.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(cluster);
}

} // namespace duckdb
