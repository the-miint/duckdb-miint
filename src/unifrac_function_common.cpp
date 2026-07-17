#include "unifrac_function_common.hpp"

#include <cmath>

#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb::unifrac_internal {

std::vector<miint::unifrac::CooRow> ReadFeatureTable(ClientContext &context, const std::string &table_name,
                                                     const std::string &caller_name, LogicalType *sample_id_type) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	// Schema probe via LIMIT 0 — surfaces missing columns or unsafe casts as a
	// binder-time error before we materialize the full table.
	auto probe = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException(
		    "%s: feature-table '%s' must expose (sample_id VARCHAR, feature_id VARCHAR, value DOUBLE): %s", caller_name,
		    table_name, probe->GetError());
	}

	// Capture the sample_id column's original SQL type (before the ::VARCHAR cast
	// above) for callers that mirror the id type onto their output. A catalog
	// metadata lookup — the same mechanism sequence_table_reader/tree_table_reader
	// use — rather than a second query; the probe above already proved sample_id
	// exists and is VARCHAR-castable, so the lookup here always finds it.
	if (sample_id_type != nullptr) {
		auto cols = GetTableOrViewColumns(context, table_name, "feature-table");
		for (idx_t i = 0; i < cols.names.size(); ++i) {
			if (StringUtil::Lower(cols.names[i]) == "sample_id") {
				*sample_id_type = cols.types[i];
				break;
			}
		}
	}

	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read feature-table '%s': %s", caller_name, table_name,
		                            result->GetError());
	}

	std::vector<miint::unifrac::CooRow> rows;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat sid_u, fid_u, val_u;
		chunk->data[0].ToUnifiedFormat(n, sid_u);
		chunk->data[1].ToUnifiedFormat(n, fid_u);
		chunk->data[2].ToUnifiedFormat(n, val_u);
		auto sid_data = UnifiedVectorFormat::GetData<string_t>(sid_u);
		auto fid_data = UnifiedVectorFormat::GetData<string_t>(fid_u);
		auto val_data = UnifiedVectorFormat::GetData<double>(val_u);
		for (idx_t i = 0; i < n; ++i) {
			const auto si = sid_u.sel->get_index(i);
			const auto fi = fid_u.sel->get_index(i);
			const auto vi = val_u.sel->get_index(i);
			if (!sid_u.validity.RowIsValid(si) || !fid_u.validity.RowIsValid(fi) || !val_u.validity.RowIsValid(vi)) {
				continue;
			}
			const double v = val_data[vi];
			if (v == 0.0 || std::isnan(v)) {
				continue; // sparse-storage invariant; UnifracSupportBiomView would drop these anyway
			}
			rows.push_back({sid_data[si].GetString(), fid_data[fi].GetString(), v});
		}
	}
	return rows;
}

int ResolveThreadsParameter(ClientContext &context, int32_t user_value, const std::string &caller_name) {
	if (user_value < 0) {
		throw BinderException("%s: threads must be >= 0 (got %d; use 0 to follow DuckDB's thread count)", caller_name,
		                      user_value);
	}
	if (user_value > 0) {
		return static_cast<int>(user_value);
	}
	// user_value == 0 → follow DuckDB. NumberOfThreads() is always >= 1.
	const auto db_threads = TaskScheduler::GetScheduler(context).NumberOfThreads();
	return NumericCast<int>(db_threads);
}

} // namespace duckdb::unifrac_internal
