#include "mmvec.hpp"
#include "mmvec_relation.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace duckdb {

namespace {

// ResolveSampleIdOutputType is named for the pattern (BIGINT/UUID pass through,
// everything else becomes VARCHAR), not for the sample_id column specifically --
// it is used here on feature_id.
using unifrac_internal::ResolveSampleIdOutputType;

// Everything the fit needs, resolved and validated at bind so a mistyped
// hyperparameter costs a binder error rather than a multi-minute fit.
struct MmvecFitBindData : public TableFunctionData {
	std::string x_table;
	std::string y_table;

	int32_t dimensions = 3;          //!< the latent dimension p
	std::string optimizer = "lbfgs"; //!< "lbfgs" or "adam"
	int64_t max_iter = 1000;         //!< L-BFGS iterations, or Adam EPOCHS
	miint::mmvec::Priors priors;     //!< scikit-bio defaults: means 0, scales 1
	miint::mmvec::AdamParams adam;   //!< scikit-bio defaults; ignored by lbfgs
	uint64_t seed = 0;               //!< the init draw; scikit-bio's None -> 0

	// Output types for the two id columns, each mirrored from its OWN relation's
	// feature_id column so a BIGINT-keyed X beside a VARCHAR-keyed Y round-trips
	// both exactly.
	LogicalType x_feature_id_type = LogicalType::VARCHAR;
	LogicalType y_feature_id_type = LogicalType::VARCHAR;
};

// Mirror one relation's feature_id type onto the output, and confirm the
// feature-table columns are present while we already hold the catalog entry.
//
// The type has to be resolved HERE rather than in ReadFeatureTable, which also
// captures id types but runs in InitGlobal -- far too late for return_types.
// community_distances_function.cpp does the same for sample_id.
LogicalType ResolveFeatureIdType(ClientContext &context, const std::string &table_name, const char *side) {
	auto cols = GetTableOrViewColumns(context, table_name, "feature-table");

	// Presence only; ReadFeatureTable's probe enforces castability when it reads.
	// Checking at bind means a mis-shaped relation fails before any I/O.
	for (const char *required : {"sample_id", "feature_id", "value"}) {
		if (!HasColumn(cols, required)) {
			throw BinderException(
			    "mmvec_fit: %s feature-table '%s' must expose (sample_id, feature_id, value); column '%s' is missing",
			    side, table_name, required);
		}
	}

	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == "feature_id") {
			return ResolveSampleIdOutputType(cols.types[i]);
		}
	}
	throw InternalException("mmvec_fit: feature_id vanished between the presence check and the type lookup");
}

unique_ptr<FunctionData> MmvecFitBind(ClientContext &context, TableFunctionBindInput &input,
                                      vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MmvecFitBindData>();
	data->x_table = input.inputs[0].GetValue<string>();
	data->y_table = input.inputs[1].GetValue<string>();

	std::string batch_norm = "unbiased";
	int64_t seed = 0;
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (kv.second.IsNull()) {
			throw BinderException("mmvec_fit: parameter '%s' cannot be NULL", key);
		}
		if (key == "dimensions") {
			data->dimensions = kv.second.GetValue<int32_t>();
		} else if (key == "optimizer") {
			data->optimizer = StringUtil::Lower(kv.second.GetValue<string>());
		} else if (key == "max_iter") {
			data->max_iter = kv.second.GetValue<int64_t>();
		} else if (key == "x_prior_mean") {
			data->priors.x_prior_mean = kv.second.GetValue<double>();
		} else if (key == "x_prior_scale") {
			data->priors.x_prior_scale = kv.second.GetValue<double>();
		} else if (key == "y_prior_mean") {
			data->priors.y_prior_mean = kv.second.GetValue<double>();
		} else if (key == "y_prior_scale") {
			data->priors.y_prior_scale = kv.second.GetValue<double>();
		} else if (key == "learning_rate") {
			data->adam.learning_rate = kv.second.GetValue<double>();
		} else if (key == "batch_size") {
			data->adam.batch_size = kv.second.GetValue<int64_t>();
		} else if (key == "beta_1") {
			data->adam.beta_1 = kv.second.GetValue<double>();
		} else if (key == "beta_2") {
			data->adam.beta_2 = kv.second.GetValue<double>();
		} else if (key == "clipnorm") {
			data->adam.clipnorm = kv.second.GetValue<double>();
		} else if (key == "batch_norm") {
			batch_norm = StringUtil::Lower(kv.second.GetValue<string>());
		} else if (key == "seed") {
			seed = kv.second.GetValue<int64_t>();
		}
		// No else: DuckDB rejects names absent from fn.named_parameters itself.
	}

	if (data->optimizer != "lbfgs" && data->optimizer != "adam") {
		throw BinderException("mmvec_fit: unknown optimizer '%s' (must be 'lbfgs' or 'adam')", data->optimizer);
	}
	if (batch_norm == "unbiased") {
		data->adam.batch_norm = miint::mmvec::BatchNorm::Unbiased;
	} else if (batch_norm == "legacy") {
		data->adam.batch_norm = miint::mmvec::BatchNorm::Legacy;
	} else {
		throw BinderException("mmvec_fit: unknown batch_norm '%s' (must be 'unbiased' or 'legacy')", batch_norm);
	}

	// These bounds duplicate the pure core's own contract checks on purpose: the
	// core validates defensively for any caller, while these give a SQL user a
	// binder error naming the parameter. Neither can drift without a test failing.
	if (data->dimensions < 1) {
		throw BinderException("mmvec_fit: dimensions must be >= 1 (got %d)", data->dimensions);
	}
	if (data->max_iter < 1) {
		throw BinderException("mmvec_fit: max_iter must be >= 1 (got %lld)", static_cast<long long>(data->max_iter));
	}
	if (!(data->priors.x_prior_scale > 0.0)) {
		throw BinderException("mmvec_fit: x_prior_scale must be > 0 (got %f)", data->priors.x_prior_scale);
	}
	if (!(data->priors.y_prior_scale > 0.0)) {
		throw BinderException("mmvec_fit: y_prior_scale must be > 0 (got %f)", data->priors.y_prior_scale);
	}
	if (!(data->adam.learning_rate > 0.0)) {
		throw BinderException("mmvec_fit: learning_rate must be > 0 (got %f)", data->adam.learning_rate);
	}
	if (data->adam.batch_size < 1) {
		throw BinderException("mmvec_fit: batch_size must be >= 1 (got %lld)",
		                      static_cast<long long>(data->adam.batch_size));
	}
	// Reported separately rather than as one combined message, so the error names
	// the offending parameter.
	if (!(data->adam.beta_1 >= 0.0 && data->adam.beta_1 < 1.0)) {
		throw BinderException("mmvec_fit: beta_1 must be in [0, 1) (got %f)", data->adam.beta_1);
	}
	if (!(data->adam.beta_2 >= 0.0 && data->adam.beta_2 < 1.0)) {
		throw BinderException("mmvec_fit: beta_2 must be in [0, 1) (got %f)", data->adam.beta_2);
	}
	if (!(data->adam.clipnorm > 0.0)) {
		throw BinderException("mmvec_fit: clipnorm must be > 0 (got %f)", data->adam.clipnorm);
	}
	if (seed < 0) {
		throw BinderException("mmvec_fit: seed must be >= 0 (got %lld)", static_cast<long long>(seed));
	}
	data->seed = static_cast<uint64_t>(seed);
	// Prior MEANS are deliberately unconstrained: a non-zero mean shifts the
	// posterior, which is a modeling choice rather than an error.

	data->x_feature_id_type = ResolveFeatureIdType(context, data->x_table, "X");
	data->y_feature_id_type = ResolveFeatureIdType(context, data->y_table, "Y");

	// One long-form relation carrying both modalities. `axis` is 1..dimensions for
	// an embedding coordinate and 0 for the bias term, which is exactly
	// scikit-bio's x_embeddings = hstack([x_main, x_bias]) column order; for
	// modality='loss' it is the 1-based evaluation ordinal.
	//
	// The five diagnostic columns are broadcast onto every row (the unifrac_pcoa
	// convention for per-result scalars) rather than encoded as extra rows, so
	// each keeps its natural type and needs no pivot to read. They are constant
	// per fit, so they cost nothing in a columnar store. n_evals is omitted as
	// derivable: count(*) WHERE modality='loss'.
	names = {"modality",  "x_feature_id", "y_feature_id", "axis",         "value",
	         "converged", "n_iter",       "final_loss",   "max_abs_grad", "message"};
	return_types = {LogicalType::VARCHAR, data->x_feature_id_type, data->y_feature_id_type, LogicalType::INTEGER,
	                LogicalType::DOUBLE,  LogicalType::BOOLEAN,    LogicalType::BIGINT,     LogicalType::DOUBLE,
	                LogicalType::DOUBLE,  LogicalType::VARCHAR};
	return std::move(data);
}

// The whole result, materialized in InitGlobal and streamed by a cursor. The fit
// is a single blocking computation that cannot be split across threads, and the
// row count is modest -- (d1 + d2) * (dimensions + 1) + one per loss evaluation,
// about 13k rows at cystic-fibrosis scale -- so there is nothing to gain from
// generating rows lazily (mirrors community_distances).
struct MmvecFitGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> x_feature_ids;
	std::vector<std::string> y_feature_ids;
	LogicalType x_type = LogicalType::VARCHAR;
	LogicalType y_type = LogicalType::VARCHAR;

	// Laid out by BuildModelRows, which is separately unit-tested by reconstructing
	// the logits from its output and comparing against the oracle -- theta's four
	// packed blocks are the one place an off-by-one stays silent.
	std::vector<miint::mmvec::ModelRow> rows;

	// Broadcast onto every row.
	bool converged = false;
	int64_t n_iter = 0;
	double final_loss = 0.0;
	double max_abs_grad = 0.0;
	std::string message;

	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

// Y as a dense row-major (n_samples x d2) matrix, reference column included,
// which is the shape the minibatch objective indexes by sample. Only the Adam
// path needs it; L-BFGS reads the summed-away sufficient statistics instead.
// A plain scatter is correct because duplicate cells were already rejected.
std::vector<double> DenseRowMajor(const miint::mmvec::SparseCounts &t) {
	std::vector<double> out(static_cast<size_t>(t.n_rows * t.n_cols), 0.0);
	for (size_t k = 0; k < t.val.size(); ++k) {
		out[static_cast<size_t>(t.row[k] * t.n_cols + t.col[k])] = t.val[k];
	}
	return out;
}

unique_ptr<GlobalTableFunctionState> MmvecFitInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<MmvecFitBindData>();
	auto gstate = make_uniq<MmvecFitGlobalState>();
	gstate->x_type = data.x_feature_id_type;
	gstate->y_type = data.y_feature_id_type;

	// ReadFeatureTable drops NULL/zero/NaN cells (the sparse feature-table
	// contract shared with community_distances and unifrac_distances) but passes
	// negatives and infinities straight through, which the core then rejects --
	// mmvec models counts, so a transformed table is a modeling error rather than
	// something to fit quietly.
	auto x_rows = unifrac_internal::ReadFeatureTable(context, data.x_table, "mmvec_fit");
	auto y_rows = unifrac_internal::ReadFeatureTable(context, data.y_table, "mmvec_fit");

	// The core reports data-dependent degeneracies (all-zero rows or columns,
	// negative or non-finite counts) and ingest reports unpairable tables, both as
	// std::invalid_argument with an already-prefixed message. Only that type is
	// caught: anything else is a bug here rather than bad user input, and should
	// surface as such instead of being relabelled.
	try {
		auto tables = miint::mmvec::IngestPairedTables(x_rows, y_rows);
		const miint::mmvec::ModelShape shape {tables.x.n_cols, tables.y.n_cols, data.dimensions};
		const auto stats = miint::mmvec::ComputeSufficientStats(tables.x, tables.y);
		const auto theta0 = miint::mmvec::InitTheta(shape, data.seed);

		miint::mmvec::Model model;
		if (data.optimizer == "lbfgs") {
			model = miint::mmvec::FitLbfgsFromInit(shape, data.priors, stats, theta0, data.max_iter);
		} else {
			miint::mmvec::MinibatchInputs mb;
			mb.n_samples = tables.x.n_rows;
			mb.x_row = tables.x.row;
			mb.x_col = tables.x.col;
			mb.x_val = tables.x.val;
			mb.y_dense = DenseRowMajor(tables.y);
			// max_iter is the EPOCH count on this path, matching scikit-bio, where
			// each epoch runs max(1, nnz / batch_size) parameter updates.
			model = miint::mmvec::FitAdam(shape, data.priors, mb, data.adam, theta0, data.max_iter, data.seed);
		}

		gstate->converged = model.converged;
		gstate->n_iter = model.n_iter;
		gstate->final_loss = model.final_loss;
		gstate->max_abs_grad = model.max_abs_grad;
		gstate->message = model.message;
		gstate->x_feature_ids = std::move(tables.x_feature_ids);
		gstate->y_feature_ids = std::move(tables.y_feature_ids);
		gstate->rows = miint::mmvec::BuildModelRows(model);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}

	return std::move(gstate);
}

void MmvecFitExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<MmvecFitGlobalState>();
	const idx_t total = g.rows.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &v_modality = output.data[0];
	auto &v_x = output.data[1];
	auto &v_y = output.data[2];
	auto axis = FlatVector::GetData<int32_t>(output.data[3]);
	auto value = FlatVector::GetData<double>(output.data[4]);
	auto converged = FlatVector::GetData<bool>(output.data[5]);
	auto n_iter = FlatVector::GetData<int64_t>(output.data[6]);
	auto final_loss = FlatVector::GetData<double>(output.data[7]);
	auto max_abs_grad = FlatVector::GetData<double>(output.data[8]);
	auto &v_message = output.data[9];

	auto modality_data = FlatVector::GetData<string_t>(v_modality);
	auto message_data = FlatVector::GetData<string_t>(v_message);

	for (idx_t r = 0; r < count; ++r) {
		const auto &row = g.rows[g.cursor + r];
		const auto id = static_cast<size_t>(row.id_index);

		// Validity is set on both branches of every id column: DuckDB reuses the
		// output vectors between calls, so a row left untouched would inherit the
		// previous chunk's flag.
		switch (row.kind) {
		case miint::mmvec::ModelRow::Kind::X:
			modality_data[r] = StringVector::AddString(v_modality, "x");
			EmitIdCell(v_x, r, g.x_feature_ids[id], g.x_type);
			FlatVector::Validity(v_y).SetInvalid(r);
			break;
		case miint::mmvec::ModelRow::Kind::Y:
			modality_data[r] = StringVector::AddString(v_modality, "y");
			FlatVector::Validity(v_x).SetInvalid(r);
			EmitIdCell(v_y, r, g.y_feature_ids[id], g.y_type);
			break;
		case miint::mmvec::ModelRow::Kind::Loss:
			modality_data[r] = StringVector::AddString(v_modality, "loss");
			FlatVector::Validity(v_x).SetInvalid(r);
			FlatVector::Validity(v_y).SetInvalid(r);
			break;
		}

		axis[r] = row.axis;
		value[r] = row.value;
		converged[r] = g.converged;
		n_iter[r] = g.n_iter;
		final_loss[r] = g.final_loss;
		max_abs_grad[r] = g.max_abs_grad;
		message_data[r] = StringVector::AddString(v_message, g.message);
	}

	g.cursor += count;
	output.SetCardinality(count);
}

} // namespace

void RegisterMmvecFit(ExtensionLoader &loader) {
	TableFunction fn("mmvec_fit", {LogicalType::VARCHAR, LogicalType::VARCHAR}, MmvecFitExecute, MmvecFitBind,
	                 MmvecFitInitGlobal);

	// scikit-bio's mmvec() signature, one named parameter each. `max_iter` serves
	// both optimizers -- L-BFGS iterations for 'lbfgs', EPOCHS for 'adam' -- which
	// is why there is no separate `epochs` parameter. There is deliberately no
	// `threads` parameter either: the core is single-threaded and the fit pins
	// Eigen to one thread, so a seeded fit is bit-reproducible.
	fn.named_parameters["dimensions"] = LogicalType::INTEGER;
	fn.named_parameters["optimizer"] = LogicalType::VARCHAR;
	fn.named_parameters["max_iter"] = LogicalType::BIGINT;
	fn.named_parameters["x_prior_mean"] = LogicalType::DOUBLE;
	fn.named_parameters["x_prior_scale"] = LogicalType::DOUBLE;
	fn.named_parameters["y_prior_mean"] = LogicalType::DOUBLE;
	fn.named_parameters["y_prior_scale"] = LogicalType::DOUBLE;
	fn.named_parameters["learning_rate"] = LogicalType::DOUBLE;
	fn.named_parameters["batch_size"] = LogicalType::BIGINT;
	fn.named_parameters["beta_1"] = LogicalType::DOUBLE;
	fn.named_parameters["beta_2"] = LogicalType::DOUBLE;
	fn.named_parameters["clipnorm"] = LogicalType::DOUBLE;
	fn.named_parameters["batch_norm"] = LogicalType::VARCHAR;
	fn.named_parameters["seed"] = LogicalType::BIGINT;

	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
