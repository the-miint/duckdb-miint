#include "woltka_ogu_function.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <atomic>

namespace duckdb {

// ── woltka_ogu() table function ─────────────────────────────────────────────
// Usage:
//   SELECT * FROM woltka_ogu('my_alignments', 'read_id');
//   SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_col');
//
// Computes Woltka OGU counts over SAM-like alignment data. With sample_id,
// runs one aggregation per sample on a per-thread Connection, parallelised
// via an atomic claim counter; without, runs a single aggregation once at
// bind time and drains its materialised result from Execute. Read IDs are
// assumed unique across samples, so the per-sample subset preserves the
// window-function result.

struct WoltkaOguData : public TableFunctionData {
	string source;
	string seq_id_col;

	bool has_sample_id = false;
	string sample_id_col;
	LogicalType sample_id_type;
	vector<Value> sample_values;

	// Non-sample path: pre-run at Bind; ownership moved to GlobalState at InitGlobal.
	unique_ptr<MaterializedQueryResult> non_sample_result;
};

struct WoltkaOguGlobalState : public GlobalTableFunctionState {
	// Per-sample mode: atomic index of the next sample to claim.
	atomic<idx_t> next_sample_idx {0};
	idx_t max_threads = 1;

	// Non-sample path: pre-run result transferred from bind data at InitGlobal.
	// Single thread drains it from Execute; no synchronization needed.
	unique_ptr<MaterializedQueryResult> non_sample_result;

	idx_t MaxThreads() const override {
		return max_threads;
	}
};

struct WoltkaOguLocalState : public LocalTableFunctionState {
	unique_ptr<Connection> conn;                // per-sample mode only
	unique_ptr<MaterializedQueryResult> result; // current sample's result
	// Keeps the fetched chunk alive while output.Reference() points into its
	// buffers. Overwritten on the next Execute call, after the upstream
	// operator has consumed the previous output (DuckDB pull-based guarantee).
	unique_ptr<DataChunk> current_chunk;
};

static string BuildAggregationSql(const string &source, const string &seq_id_col) {
	auto q_src = KeywordHelper::WriteOptionallyQuoted(source);
	auto q_seq = KeywordHelper::WriteOptionallyQuoted(seq_id_col);
	return "WITH base AS ("
	       "  SELECT DISTINCT "
	       "    " +
	       q_seq +
	       " AS qid, "
	       "    reference AS feature_id, "
	       "    alignment_is_read1(flags::USMALLINT) AS is_fwd "
	       "  FROM " +
	       q_src +
	       "), "
	       "with_counts AS ("
	       "  SELECT feature_id, "
	       "         1.0 / COUNT(*) OVER (PARTITION BY qid, is_fwd) AS local_value "
	       "  FROM base) "
	       "SELECT feature_id, SUM(local_value) AS value "
	       "FROM with_counts "
	       "GROUP BY feature_id";
}

static unique_ptr<MaterializedQueryResult> RunGlobalAggregation(Connection &conn, const string &source,
                                                                const string &seq_id_col) {
	auto sql = BuildAggregationSql(source, seq_id_col);
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("woltka_ogu query failed: %s\nGenerated SQL: %s", result->GetError(), sql);
	}
	return result;
}

// Run the woltka_ogu pipeline for a single sample value. Creates a TEMP VIEW
// scoped to `conn`. Each thread has its own Connection, so names don't collide.
static unique_ptr<MaterializedQueryResult> RunSampleAggregation(Connection &conn, const string &source,
                                                                const string &seq_id_col, const string &sample_col,
                                                                const Value &sample_value) {
	auto q_src = KeywordHelper::WriteOptionallyQuoted(source);
	auto q_sample = KeywordHelper::WriteOptionallyQuoted(sample_col);
	// ToSQLString handles all Value types (integers, timestamps, strings) safely.
	auto sample_literal = sample_value.ToSQLString();

	auto view_sql = "CREATE OR REPLACE TEMP VIEW __woltka_per_sample AS SELECT * FROM " + q_src + " WHERE CAST(" +
	                q_sample + " AS VARCHAR) = CAST(" + sample_literal + " AS VARCHAR)";
	auto view_result = conn.Query(view_sql);
	if (view_result->HasError()) {
		throw InvalidInputException("woltka_ogu: failed to create per-sample view: %s", view_result->GetError());
	}

	auto inner_sql = BuildAggregationSql("__woltka_per_sample", seq_id_col);
	auto wrapped_sql = "SELECT " + sample_literal + " AS " + q_sample + ", __q.* FROM (" + inner_sql + ") __q";

	auto result = conn.Query(wrapped_sql);
	if (result->HasError()) {
		throw InvalidInputException("woltka_ogu query failed: %s\nGenerated SQL: %s", result->GetError(), wrapped_sql);
	}
	return result;
}

static unique_ptr<FunctionData> WoltkaOguBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<WoltkaOguData>();
	data->source = input.inputs[0].GetValue<string>();
	data->seq_id_col = input.inputs[1].GetValue<string>();

	if (data->source.empty()) {
		throw InvalidInputException("woltka_ogu: relation name must not be empty");
	}
	if (data->seq_id_col.empty()) {
		throw InvalidInputException("woltka_ogu: sequence_id_field must not be empty");
	}

	auto it = input.named_parameters.find("sample_id");
	if (it != input.named_parameters.end()) {
		data->sample_id_col = it->second.GetValue<string>();
		if (data->sample_id_col.empty()) {
			throw InvalidInputException("woltka_ogu: sample_id column name must not be empty");
		}
		// Reject column names that collide with the fixed output schema.
		// DuckDB identifiers are case-insensitive.
		auto col_lower = StringUtil::Lower(data->sample_id_col);
		if (col_lower == "feature_id" || col_lower == "value") {
			throw InvalidInputException("woltka_ogu: sample_id column name '%s' collides with an output column "
			                            "(feature_id, value); rename the column or alias it in a view",
			                            data->sample_id_col);
		}
		data->has_sample_id = true;
	}

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	auto q_src = KeywordHelper::WriteOptionallyQuoted(data->source);
	auto q_seq = KeywordHelper::WriteOptionallyQuoted(data->seq_id_col);

	// Validate source resolves AND required columns exist AND their types are compatible
	// with the aggregation (reference → VARCHAR, flags → USMALLINT). LIMIT 0 binds the
	// casts without scanning rows.
	auto probe = conn.Query("SELECT " + q_seq + ", reference::VARCHAR, flags::USMALLINT FROM " + q_src + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("woltka_ogu: source, required columns (reference, flags), or their types "
		                            "are not compatible: %s",
		                            probe->GetError());
	}

	if (data->has_sample_id) {
		auto q_sample = KeywordHelper::WriteOptionallyQuoted(data->sample_id_col);

		auto sample_probe = conn.Query("SELECT " + q_sample + " FROM " + q_src + " LIMIT 0");
		if (sample_probe->HasError()) {
			throw InvalidInputException("woltka_ogu: sample_id column '%s' not found", data->sample_id_col);
		}
		data->sample_id_type = sample_probe->types[0];

		// Single scan: distinct sample values, NULLs first so we can reject up-front.
		auto distinct_result =
		    conn.Query("SELECT DISTINCT " + q_sample + " FROM " + q_src + " ORDER BY " + q_sample + " NULLS FIRST");
		if (distinct_result->HasError()) {
			throw InvalidInputException("woltka_ogu: failed to query sample_id values: %s",
			                            distinct_result->GetError());
		}
		auto &materialized = distinct_result->Cast<MaterializedQueryResult>();
		while (auto chunk = materialized.Fetch()) {
			for (idx_t i = 0; i < chunk->size(); i++) {
				auto val = chunk->data[0].GetValue(i);
				if (val.IsNull()) {
					throw InvalidInputException("NULL values in sample_id column '%s'", data->sample_id_col);
				}
				data->sample_values.push_back(std::move(val));
			}
		}
		if (data->sample_values.empty()) {
			throw InvalidInputException("sample_id column '%s' has no non-NULL values", data->sample_id_col);
		}

		names.push_back(data->sample_id_col);
		return_types.push_back(data->sample_id_type);
	} else {
		// Non-sample path: pre-run the aggregation once here. Ownership of the result
		// transfers to GlobalState at InitGlobal, matching MassQLFunction's pattern.
		data->non_sample_result = RunGlobalAggregation(conn, data->source, data->seq_id_col);
	}

	names.emplace_back("feature_id");
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("value");
	return_types.emplace_back(LogicalType::DOUBLE);

	return data;
}

static unique_ptr<GlobalTableFunctionState> WoltkaOguInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	// CastNoConst: we move result ownership out of bind data into global state.
	// InitGlobal is called exactly once before any Execute thread starts, so this is safe.
	auto &data = input.bind_data->CastNoConst<WoltkaOguData>();
	auto gstate = make_uniq<WoltkaOguGlobalState>();
	if (data.has_sample_id) {
		idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		gstate->max_threads = std::max<idx_t>(1, std::min<idx_t>(db_threads, data.sample_values.size()));
	} else {
		gstate->max_threads = 1;
		gstate->non_sample_result = std::move(data.non_sample_result);
	}
	return gstate;
}

static unique_ptr<LocalTableFunctionState> WoltkaOguInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                              GlobalTableFunctionState * /*global_state*/) {
	auto &data = input.bind_data->Cast<WoltkaOguData>();
	auto lstate = make_uniq<WoltkaOguLocalState>();
	if (data.has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context.client);
		lstate->conn = make_uniq<Connection>(db);
	}
	return lstate;
}

static void WoltkaOguExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &data = input.bind_data->Cast<WoltkaOguData>();
	auto &gstate = input.global_state->Cast<WoltkaOguGlobalState>();
	auto &lstate = input.local_state->Cast<WoltkaOguLocalState>();

	if (!data.has_sample_id) {
		// Single thread drains the pre-run result from global state.
		lstate.current_chunk = gstate.non_sample_result->Fetch();
		if (lstate.current_chunk && lstate.current_chunk->size() > 0) {
			output.Reference(*lstate.current_chunk);
			return;
		}
		output.SetCardinality(0);
		return;
	}

	// Per-sample: each thread atomically claims sample indices.
	while (true) {
		if (lstate.result) {
			lstate.current_chunk = lstate.result->Fetch();
			if (lstate.current_chunk && lstate.current_chunk->size() > 0) {
				output.Reference(*lstate.current_chunk);
				return;
			}
			lstate.result.reset();
		}
		idx_t sample_idx = gstate.next_sample_idx.fetch_add(1, std::memory_order_relaxed);
		if (sample_idx >= data.sample_values.size()) {
			output.SetCardinality(0);
			return;
		}
		lstate.result = RunSampleAggregation(*lstate.conn, data.source, data.seq_id_col, data.sample_id_col,
		                                     data.sample_values[sample_idx]);
	}
}

void WoltkaOguFunction::Register(ExtensionLoader &loader) {
	TableFunction fn("woltka_ogu", {LogicalType::VARCHAR, LogicalType::VARCHAR}, WoltkaOguExecute, WoltkaOguBind,
	                 WoltkaOguInitGlobal, WoltkaOguInitLocal);
	fn.named_parameters["sample_id"] = LogicalType::VARCHAR;
	fn.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
