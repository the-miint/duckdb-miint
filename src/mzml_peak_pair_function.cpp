#include "mzml_peak_pair_function.hpp"
#include "catalog_utils.hpp"
#include "formula_parser.hpp"

#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <atomic>

namespace duckdb {

// ── mzml_peak_pair() table function ─────────────────────────────────────────
// Usage: SELECT * FROM mzml_peak_pair('spectra', 'Fe')
//
// Finds MS2 spectra containing a peak pair where one peak is at m/z = X and
// another is at m/z = 2*X - formula(formula_str), within 0.1 Da tolerance.
// Returns all peaks from matching spectra (composable with mzml_scaninfo).
//
// This is a C++ table function replacement for the former SQL macro. The macro
// suffered from severe performance issues because the recursive CTE's
// `SELECT DISTINCT mz FROM ms2` was re-evaluated at every recursion step (~3k
// iterations). This implementation materializes both the MS2 peaks and their
// distinct m/z values into temp tables, reducing the query from ~18s to ~1s.

static std::atomic<uint64_t> peak_pair_invocation_counter {0};

struct MzmlPeakPairData : public TableFunctionData {
	unique_ptr<MaterializedQueryResult> result;
};

struct MzmlPeakPairGlobalState : public GlobalTableFunctionState {
	unique_ptr<MaterializedQueryResult> result;
	idx_t MaxThreads() const override {
		return 1;
	}
};

struct MzmlPeakPairLocalState : public LocalTableFunctionState {
	unique_ptr<DataChunk> current_chunk;
};

static void ExtractSchema(MaterializedQueryResult &result, vector<LogicalType> &return_types, vector<string> &names) {
	for (idx_t i = 0; i < result.ColumnCount(); i++) {
		names.push_back(result.ColumnName(i));
		return_types.push_back(result.types[i]);
	}
}

static unique_ptr<FunctionData> PeakPairBind(ClientContext &context, TableFunctionBindInput &input,
                                             vector<LogicalType> &return_types, vector<string> &names) {
	auto data = make_uniq<MzmlPeakPairData>();

	auto relation = input.inputs[0].GetValue<string>();
	auto formula_str = input.inputs[1].GetValue<string>();

	auto quoted_relation = KeywordHelper::WriteOptionallyQuoted(relation);

	// Compute formula mass in C++ — avoids embedding user input into SQL.
	double formula_mass;
	try {
		formula_mass = miint::ParseFormula(formula_str);
	} catch (const std::exception &e) {
		throw InvalidInputException("mzml_peak_pair: invalid formula '%s': %s", formula_str, e.what());
	}
	auto mass_literal = StringUtil::Format("%.17g", formula_mass);

	// Unique suffix per invocation to avoid temp table name collisions.
	auto suffix = std::to_string(peak_pair_invocation_counter.fetch_add(1, std::memory_order_relaxed));
	auto ms2_table = "__peak_pair_ms2_" + suffix;
	auto mz_table = "__peak_pair_mz_" + suffix;

	// Inherits the caller's TEMP catalog so `relation` may itself be a TEMP table or
	// view (#193). The two tables created below are already uniquely named by the
	// counter above, which is the first half of what inheriting requires; the
	// HelperTempRelation guards are the second half — they now live in the caller's
	// session and no longer vanish when this connection is torn down at the end of
	// Bind. Both guards outlive the final query below, whose result is materialized
	// into bind data before they run.
	auto conn = MakeReadOnlyHelperConnection(context);

	// Materialize MS2 peaks into a temp table.
	auto mat_sql = "CREATE TEMP TABLE " + ms2_table +
	               " AS "
	               "SELECT * FROM mzml_peaks(" +
	               quoted_relation + ") WHERE ms_level = 2 AND intensity > 0";
	auto mat_result = conn.Query(mat_sql);
	if (mat_result->HasError()) {
		throw InvalidInputException("mzml_peak_pair: failed to materialize MS2 peaks: %s", mat_result->GetError());
	}
	HelperTempRelation ms2_guard(conn, ms2_table);

	// Materialize distinct m/z values. The recursive CTE re-evaluates
	// SELECT DISTINCT mz at every recursion step (~3k iterations); pre-materializing
	// this set is the key optimization (18s → <1s on real data).
	auto distinct_sql = "CREATE TEMP TABLE " + mz_table +
	                    " AS "
	                    "SELECT DISTINCT mz FROM " +
	                    ms2_table;
	auto distinct_result = conn.Query(distinct_sql);
	if (distinct_result->HasError()) {
		throw InvalidInputException("mzml_peak_pair: failed to materialize distinct mz: %s",
		                            distinct_result->GetError());
	}
	HelperTempRelation mz_guard(conn, mz_table);

	// Run the recursive query against the materialized tables.
	auto query_sql = "WITH RECURSIVE "
	                 "x_candidates(x_val, next_min) AS ( "
	                 "    (SELECT mz, mz + 0.05 FROM " +
	                 mz_table +
	                 " ORDER BY mz LIMIT 1) "
	                 "    UNION ALL "
	                 "    (SELECT s.mz, s.mz + 0.05 "
	                 "     FROM x_candidates g "
	                 "     JOIN " +
	                 mz_table +
	                 " s ON s.mz >= g.next_min "
	                 "     ORDER BY s.mz "
	                 "     LIMIT 1) "
	                 ") "
	                 "SELECT * FROM " +
	                 ms2_table +
	                 " WHERE spectrum_index IN ( "
	                 "    SELECT DISTINCT p1.spectrum_index "
	                 "    FROM x_candidates xc "
	                 "    JOIN " +
	                 ms2_table +
	                 " p1 ON p1.mz > xc.x_val - 0.1 AND p1.mz < xc.x_val + 0.1 "
	                 "    JOIN " +
	                 ms2_table +
	                 " p2 ON p2.spectrum_index = p1.spectrum_index "
	                 "        AND p2.mz > 2.0 * xc.x_val - " +
	                 mass_literal +
	                 " - 0.1 "
	                 "        AND p2.mz < 2.0 * xc.x_val - " +
	                 mass_literal +
	                 " + 0.1 "
	                 ")";
	auto query_result = conn.Query(query_sql);
	if (query_result->HasError()) {
		throw InvalidInputException("mzml_peak_pair: query failed: %s", query_result->GetError());
	}

	data->result = unique_ptr_cast<QueryResult, MaterializedQueryResult>(std::move(query_result));
	ExtractSchema(*data->result, return_types, names);

	return data;
}

static unique_ptr<GlobalTableFunctionState> PeakPairInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->CastNoConst<MzmlPeakPairData>();
	auto gstate = make_uniq<MzmlPeakPairGlobalState>();
	gstate->result = std::move(data.result);
	return gstate;
}

static unique_ptr<LocalTableFunctionState> PeakPairInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                             GlobalTableFunctionState *) {
	return make_uniq<MzmlPeakPairLocalState>();
}

static void PeakPairExecute(ClientContext &context, TableFunctionInput &input, DataChunk &output) {
	auto &gstate = input.global_state->Cast<MzmlPeakPairGlobalState>();
	auto &lstate = input.local_state->Cast<MzmlPeakPairLocalState>();

	lstate.current_chunk = gstate.result->Fetch();
	if (lstate.current_chunk && lstate.current_chunk->size() > 0) {
		output.Reference(*lstate.current_chunk);
		return;
	}
	output.SetCardinality(0);
}

void MzmlPeakPairFunction::Register(ExtensionLoader &loader) {
	TableFunction func("mzml_peak_pair", {LogicalType::VARCHAR, LogicalType::VARCHAR}, PeakPairExecute, PeakPairBind,
	                   PeakPairInitGlobal, PeakPairInitLocal);
	func.order_preservation_type = OrderPreservationType::NO_ORDER;
	loader.RegisterFunction(func);
}

} // namespace duckdb
