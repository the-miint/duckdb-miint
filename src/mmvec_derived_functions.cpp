#include "mmvec.hpp"
#include "mmvec_relation.hpp"

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "unifrac_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace duckdb {

namespace {

// Named for the pattern (BIGINT/UUID pass through, everything else becomes
// VARCHAR), not for the sample_id column specifically -- used here on feature_id.
using unifrac_internal::ResolveSampleIdOutputType;

// ---------------------------------------------------------------------------
// Shared by every function that consumes a model relation. All of them take a
// fitted model as input (locked decision D1: fit once, then read the model
// table), so resolving its id types at bind and reading its cells in InitGlobal
// is the one piece of plumbing they have in common.
//
// What is deliberately NOT shared, decided at the M4 review rather than left
// unstated:
//
//  * The three bind / InitGlobal / Execute skeletons. They look alike because
//    every DuckDB table function does, but they differ in arity, in what they
//    compute and in their cursor arithmetic -- ranks walks a (d1 x d2) grid,
//    predict an (n_samples x d2) one, score emits a single row. A template over
//    those would be longer than the three bodies and would hide the one line in
//    each that matters. The same shallow similarity is repeated across every
//    table function in this codebase; conforming to it beats forking style here.
//  * `ResolveFeatureTableIdType` is a near-duplicate of `ResolveFeatureIdType` in
//    mmvec_fit_function.cpp, differing in taking the column name and the caller
//    as parameters instead of hard-coding "feature_id"/"mmvec_fit". The M6
//    simplification pass looked at unifying the two and decided against it, for a
//    reason that was not obvious beforehand: they differ on a THIRD axis, the X/Y
//    side infix in mmvec_fit's message ("Y feature-table '...' must expose ..."),
//    which test/sql/mmvec_fit.test asserts precisely because it is what makes the
//    Y-side case more than a repeat of the X-side one. Preserving it means a
//    fifth parameter on a helper with two callers, in two different TUs, whose
//    only natural shared home is unifrac_function_common.hpp -- i.e. the
//    cross-feature consolidation that is deferred anyway. That consolidation
//    (this, plus the equivalents in community_distances, cluster_kmeans and the
//    unifrac functions, plus the three dictionary builders) remains worth doing
//    as its own PR; it restructures the shared reader layer and does not belong
//    in a feature commit.
// ---------------------------------------------------------------------------

// The two id types a model relation carries, mirrored onto whatever output the
// caller emits so a BIGINT- or UUID-keyed model joins back to typed metadata
// without a cast.
struct ModelIdTypes {
	LogicalType x = LogicalType::VARCHAR;
	LogicalType y = LogicalType::VARCHAR;
};

// Presence only; whether the columns actually convert is the reader's business.
// Checking at bind turns a mis-shaped relation into a binder error instead of a
// failure after the model has been materialized -- and every function that takes a
// model checks it the same way, including the ones with no id type to mirror.
TableOrViewColumns RequireModelColumns(ClientContext &context, const std::string &table_name, const char *caller) {
	auto cols = GetTableOrViewColumns(context, table_name, "model relation");
	for (const char *required : {"modality", "x_feature_id", "y_feature_id", "axis", "value"}) {
		if (!HasColumn(cols, required)) {
			throw BinderException("%s: model relation '%s' must expose (modality, x_feature_id, y_feature_id, axis, "
			                      "value) as produced by mmvec_fit; column '%s' is missing",
			                      caller, table_name, required);
		}
	}
	return cols;
}

TableOrViewColumns RequireFeatureTableColumns(ClientContext &context, const std::string &table_name,
                                              const char *caller) {
	auto cols = GetTableOrViewColumns(context, table_name, "feature-table");
	for (const char *required : {"sample_id", "feature_id", "value"}) {
		if (!HasColumn(cols, required)) {
			throw BinderException("%s: feature-table '%s' must expose (sample_id, feature_id, value); column '%s' is "
			                      "missing",
			                      caller, table_name, required);
		}
	}
	return cols;
}

// Resolved at bind, where return_types must be declared -- too early for the
// relation to be read, so this is a catalog metadata lookup rather than a query
// (the same split mmvec_fit uses for its two feature tables).
ModelIdTypes ResolveModelIdTypes(ClientContext &context, const std::string &table_name, const char *caller) {
	auto cols = RequireModelColumns(context, table_name, caller);

	ModelIdTypes types;
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		const auto name = StringUtil::Lower(cols.names[i]);
		if (name == "x_feature_id") {
			types.x = ResolveSampleIdOutputType(cols.types[i]);
		} else if (name == "y_feature_id") {
			types.y = ResolveSampleIdOutputType(cols.types[i]);
		}
	}
	return types;
}

// Read a model relation into cells. Ids are carried as canonical text, exactly as
// the feature-table readers do, so BIGINT and UUID keys survive the round trip.
//
// NULLs: a NULL `modality`, `axis` or `value` is rejected here, because none of
// them has a meaning -- an unclassifiable row, or a parameter that is not a
// number. A NULL feature id is passed through as nullopt so ParseModelCells
// reports it, keeping that message in one place.
//
// The id is selected BY MODALITY -- x_feature_id for 'x', y_feature_id for 'y' --
// and never "whichever is non-NULL", which on a mislabelled row would read a Y id
// as an X feature name and invent a feature rather than complain.
std::vector<miint::mmvec::ModelCell> ReadModelCells(ClientContext &context, const std::string &table_name,
                                                    const char *caller) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const std::string projection = "SELECT modality::VARCHAR, x_feature_id::VARCHAR, y_feature_id::VARCHAR, "
	                               "axis::INTEGER, value::DOUBLE FROM " +
	                               qname;

	auto probe = conn.Query(projection + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: model relation '%s' must expose (modality VARCHAR, x_feature_id, "
		                            "y_feature_id, axis INTEGER, value DOUBLE): %s",
		                            caller, table_name, probe->GetError());
	}
	auto result = conn.Query(projection);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read model relation '%s': %s", caller, table_name,
		                            result->GetError());
	}

	std::vector<miint::mmvec::ModelCell> cells;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat mod_u, x_u, y_u, axis_u, val_u;
		chunk->data[0].ToUnifiedFormat(n, mod_u);
		chunk->data[1].ToUnifiedFormat(n, x_u);
		chunk->data[2].ToUnifiedFormat(n, y_u);
		chunk->data[3].ToUnifiedFormat(n, axis_u);
		chunk->data[4].ToUnifiedFormat(n, val_u);
		auto mod_data = UnifiedVectorFormat::GetData<string_t>(mod_u);
		auto x_data = UnifiedVectorFormat::GetData<string_t>(x_u);
		auto y_data = UnifiedVectorFormat::GetData<string_t>(y_u);
		auto axis_data = UnifiedVectorFormat::GetData<int32_t>(axis_u);
		auto val_data = UnifiedVectorFormat::GetData<double>(val_u);

		for (idx_t i = 0; i < n; ++i) {
			const auto mi = mod_u.sel->get_index(i);
			const auto ai = axis_u.sel->get_index(i);
			const auto vi = val_u.sel->get_index(i);
			if (!mod_u.validity.RowIsValid(mi)) {
				throw InvalidInputException("%s: model relation '%s' has a row with a NULL modality", caller,
				                            table_name);
			}
			miint::mmvec::ModelCell cell;
			cell.modality = StringUtil::Lower(mod_data[mi].GetString());
			if (!axis_u.validity.RowIsValid(ai) || !val_u.validity.RowIsValid(vi)) {
				throw InvalidInputException("%s: model relation '%s' has a modality '%s' row with a NULL axis or value",
				                            caller, table_name, cell.modality);
			}
			cell.axis = axis_data[ai];
			cell.value = val_data[vi];

			if (cell.modality == "x") {
				const auto xi = x_u.sel->get_index(i);
				if (x_u.validity.RowIsValid(xi)) {
					cell.feature_id = x_data[xi].GetString();
				}
			} else if (cell.modality == "y") {
				const auto yi = y_u.sel->get_index(i);
				if (y_u.validity.RowIsValid(yi)) {
					cell.feature_id = y_data[yi].GetString();
				}
			}
			cells.push_back(std::move(cell));
		}
	}
	return cells;
}

// The core reports every structural problem as std::invalid_argument with an
// already-prefixed message. Only that type is caught: anything else is a bug here
// rather than bad user input, and should surface as such instead of being
// relabelled (mirrors mmvec_fit_function.cpp).
miint::mmvec::ParsedModel ParseOrThrow(const std::vector<miint::mmvec::ModelCell> &cells) {
	try {
		return miint::mmvec::ParseModelCells(cells);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
}

// ---------------------------------------------------------------------------
// mmvec_ranks
// ---------------------------------------------------------------------------

struct MmvecRanksBindData : public TableFunctionData {
	std::string model_table;
	ModelIdTypes types;
};

unique_ptr<FunctionData> MmvecRanksBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MmvecRanksBindData>();
	data->model_table = input.inputs[0].GetValue<string>();
	data->types = ResolveModelIdTypes(context, data->model_table, "mmvec_ranks");

	// `rank` and `prob` side by side rather than as two functions: Probs is
	// algebraically softmax(Ranks), so both come from one read of the model, and a
	// separate mmvec_probs would make the caller pay for that read twice.
	//
	// They do NOT share a logit computation, though an earlier version of this
	// comment claimed they did. Ranks() and Probs() each call ComputeLogits, so
	// this function runs the (d1 x p)(p x d2-1) product twice -- measured at 4-5 ms
	// of a 68 ms mmvec_ranks over the cystic-fibrosis model. Left alone
	// deliberately: collapsing it needs a core entry point returning both, with
	// exactly one caller, and the two must stay separately reachable because
	// softmax(logits) and softmax(ranks) differ in the last ulps and a test pins
	// their agreement.
	names = {"x_feature_id", "y_feature_id", "rank", "prob"};
	return_types = {data->types.x, data->types.y, LogicalType::DOUBLE, LogicalType::DOUBLE};
	return std::move(data);
}

// Both matrices materialized up front, following mmvec_fit. d1*d2 rows: 1.26M and
// ~20 MB of values at cystic-fibrosis scale (d1 2720, d2 462), which is why this is
// a plain cursor and not worth streaming per X feature.
struct MmvecRanksGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> x_feature_ids;
	std::vector<std::string> y_feature_ids;
	LogicalType x_type = LogicalType::VARCHAR;
	LogicalType y_type = LogicalType::VARCHAR;
	std::vector<double> ranks; //!< (d1 x d2) row-major
	std::vector<double> probs; //!< (d1 x d2) row-major, same indexing
	int64_t d2 = 0;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<GlobalTableFunctionState> MmvecRanksInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<MmvecRanksBindData>();
	auto gstate = make_uniq<MmvecRanksGlobalState>();
	gstate->x_type = data.types.x;
	gstate->y_type = data.types.y;

	const auto model = ParseOrThrow(ReadModelCells(context, data.model_table, "mmvec_ranks"));
	try {
		gstate->ranks = miint::mmvec::Ranks(model.shape, model.theta);
		gstate->probs = miint::mmvec::Probs(model.shape, model.theta);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
	gstate->d2 = model.shape.n_features_y;
	gstate->x_feature_ids = model.x_feature_ids;
	gstate->y_feature_ids = model.y_feature_ids;
	return std::move(gstate);
}

void MmvecRanksExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<MmvecRanksGlobalState>();
	const idx_t total = g.ranks.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &v_x = output.data[0];
	auto &v_y = output.data[1];
	auto rank = FlatVector::GetData<double>(output.data[2]);
	auto prob = FlatVector::GetData<double>(output.data[3]);

	for (idx_t r = 0; r < count; ++r) {
		const idx_t k = g.cursor + r;
		// Row-major over (X feature, Y feature), so the cell index carries both ids.
		const idx_t i = k / static_cast<idx_t>(g.d2);
		const idx_t j = k % static_cast<idx_t>(g.d2);
		EmitIdCell(v_x, r, g.x_feature_ids[i], g.x_type);
		EmitIdCell(v_y, r, g.y_feature_ids[j], g.y_type);
		rank[r] = g.ranks[k];
		prob[r] = g.probs[k];
	}
	g.cursor += count;
	output.SetCardinality(count);
}

// ---------------------------------------------------------------------------
// mmvec_predict
// ---------------------------------------------------------------------------

// Mirror one column of a feature-table onto the output. Same bind-time necessity
// as ResolveModelIdTypes: ReadFeatureTable also captures id types, but runs in
// InitGlobal, far too late for return_types.
LogicalType ResolveFeatureTableIdType(ClientContext &context, const std::string &table_name, const char *column,
                                      const char *caller) {
	auto cols = RequireFeatureTableColumns(context, table_name, caller);
	for (idx_t i = 0; i < cols.names.size(); ++i) {
		if (StringUtil::Lower(cols.names[i]) == column) {
			return ResolveSampleIdOutputType(cols.types[i]);
		}
	}
	throw InternalException("%s: column '%s' vanished between the presence check and the type lookup", caller, column);
}

struct MmvecPredictBindData : public TableFunctionData {
	std::string model_table;
	std::string x_table;
	//! sample_id comes from the NEW X table -- these are samples the model has never
	//! seen, and they are the point -- while y_feature_id comes from the model.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	ModelIdTypes types;
};

unique_ptr<FunctionData> MmvecPredictBind(ClientContext &context, TableFunctionBindInput &input,
                                          vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MmvecPredictBindData>();
	data->model_table = input.inputs[0].GetValue<string>();
	data->x_table = input.inputs[1].GetValue<string>();
	data->types = ResolveModelIdTypes(context, data->model_table, "mmvec_predict");
	data->sample_id_type = ResolveFeatureTableIdType(context, data->x_table, "sample_id", "mmvec_predict");

	names = {"sample_id", "y_feature_id", "proportion"};
	return_types = {data->sample_id_type, data->types.y, LogicalType::DOUBLE};
	return std::move(data);
}

struct MmvecPredictGlobalState : public GlobalTableFunctionState {
	std::vector<std::string> sample_ids;
	std::vector<std::string> y_feature_ids;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType y_type = LogicalType::VARCHAR;
	std::vector<double> proportions; //!< (n_samples x d2) row-major
	int64_t d2 = 0;
	idx_t cursor = 0;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<GlobalTableFunctionState> MmvecPredictInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<MmvecPredictBindData>();
	auto gstate = make_uniq<MmvecPredictGlobalState>();
	gstate->sample_id_type = data.sample_id_type;
	gstate->y_type = data.types.y;

	const auto model = ParseOrThrow(ReadModelCells(context, data.model_table, "mmvec_predict"));
	auto x_rows = unifrac_internal::ReadFeatureTable(context, data.x_table, "mmvec_predict");

	try {
		// Aligned against the MODEL's X dictionary, never a fresh one built from the
		// held-out table: a fresh dictionary would produce in-range indices naming
		// different microbes, which is a wrong answer rather than a failure.
		const auto aligned = miint::mmvec::AlignXToModel(x_rows, model.x_feature_ids);
		gstate->proportions = miint::mmvec::Predict(model.shape, model.theta, aligned.counts);
		gstate->sample_ids = aligned.sample_ids;
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
	gstate->d2 = model.shape.n_features_y;
	gstate->y_feature_ids = model.y_feature_ids;
	return std::move(gstate);
}

void MmvecPredictExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<MmvecPredictGlobalState>();
	const idx_t total = g.proportions.size();
	const idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - g.cursor);

	auto &v_sample = output.data[0];
	auto &v_y = output.data[1];
	auto proportion = FlatVector::GetData<double>(output.data[2]);

	for (idx_t r = 0; r < count; ++r) {
		const idx_t k = g.cursor + r;
		const idx_t n = k / static_cast<idx_t>(g.d2);
		const idx_t j = k % static_cast<idx_t>(g.d2);
		EmitIdCell(v_sample, r, g.sample_ids[n], g.sample_id_type);
		EmitIdCell(v_y, r, g.y_feature_ids[j], g.y_type);
		proportion[r] = g.proportions[k];
	}
	g.cursor += count;
	output.SetCardinality(count);
}

// ---------------------------------------------------------------------------
// mmvec_score
// ---------------------------------------------------------------------------

struct MmvecScoreBindData : public TableFunctionData {
	std::string model_table;
	std::string x_table;
	std::string y_table;
};

unique_ptr<FunctionData> MmvecScoreBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MmvecScoreBindData>();
	data->model_table = input.inputs[0].GetValue<string>();
	data->x_table = input.inputs[1].GetValue<string>();
	data->y_table = input.inputs[2].GetValue<string>();

	// Nothing here mirrors an id type -- the output is one number -- but all three
	// relations are still validated at bind, so a mis-shaped input is a binder error
	// exactly as it is for mmvec_ranks and mmvec_predict.
	RequireModelColumns(context, data->model_table, "mmvec_score");
	RequireFeatureTableColumns(context, data->x_table, "mmvec_score");
	RequireFeatureTableColumns(context, data->y_table, "mmvec_score");

	names = {"q_squared"};
	return_types = {LogicalType::DOUBLE};
	return std::move(data);
}

struct MmvecScoreGlobalState : public GlobalTableFunctionState {
	double q_squared = 0.0;
	bool emitted = false;
	idx_t MaxThreads() const override {
		return 1;
	}
};

unique_ptr<GlobalTableFunctionState> MmvecScoreInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<MmvecScoreBindData>();
	auto gstate = make_uniq<MmvecScoreGlobalState>();

	const auto model = ParseOrThrow(ReadModelCells(context, data.model_table, "mmvec_score"));
	auto x_rows = unifrac_internal::ReadFeatureTable(context, data.x_table, "mmvec_score");
	auto y_rows = unifrac_internal::ReadFeatureTable(context, data.y_table, "mmvec_score");

	try {
		// Both tables against the MODEL's dictionaries and one shared sample
		// dictionary, so row n of each is the same sample -- Score compares them cell
		// by cell.
		const auto aligned = miint::mmvec::AlignPairedToModel(x_rows, y_rows, model.x_feature_ids, model.y_feature_ids);
		gstate->q_squared = miint::mmvec::Score(model.shape, model.theta, aligned.x, aligned.y);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s", e.what());
	}
	return std::move(gstate);
}

void MmvecScoreExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &g = data_p.global_state->Cast<MmvecScoreGlobalState>();
	if (g.emitted) {
		output.SetCardinality(0);
		return;
	}
	FlatVector::GetData<double>(output.data[0])[0] = g.q_squared;
	g.emitted = true;
	output.SetCardinality(1);
}

} // namespace

void RegisterMmvecRanks(ExtensionLoader &loader) {
	TableFunction fn("mmvec_ranks", {LogicalType::VARCHAR}, MmvecRanksExecute, MmvecRanksBind, MmvecRanksInitGlobal);
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

void RegisterMmvecPredict(ExtensionLoader &loader) {
	TableFunction fn("mmvec_predict", {LogicalType::VARCHAR, LogicalType::VARCHAR}, MmvecPredictExecute,
	                 MmvecPredictBind, MmvecPredictInitGlobal);
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

void RegisterMmvecScore(ExtensionLoader &loader) {
	TableFunction fn("mmvec_score", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                 MmvecScoreExecute, MmvecScoreBind, MmvecScoreInitGlobal);
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
