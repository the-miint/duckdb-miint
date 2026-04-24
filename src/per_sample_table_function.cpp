#include "per_sample_table_function.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

namespace duckdb {

void DiscoverSamples(Connection &conn, const string &source_relation, const string &sample_id_col,
                     const std::vector<std::string> &reserved_lowercase_output_names, const string &fn_label,
                     PerSampleBindInfo &out) {
	if (sample_id_col.empty()) {
		throw InvalidInputException("%s: sample_id column name must not be empty", fn_label);
	}
	// Callers that don't want a collision check pass {} and skip this loop entirely.
	auto col_lower = StringUtil::Lower(sample_id_col);
	for (const auto &reserved : reserved_lowercase_output_names) {
		if (col_lower == reserved) {
			// Construct a duckdb::vector<std::string> for StringUtil::Join which expects it.
			duckdb::vector<std::string> joined(reserved_lowercase_output_names.begin(),
			                                   reserved_lowercase_output_names.end());
			throw InvalidInputException("%s: sample_id column name '%s' collides with an output column (%s); "
			                            "rename the column or alias it in a view",
			                            fn_label, sample_id_col, StringUtil::Join(joined, ", "));
		}
	}

	auto q_src = KeywordHelper::WriteOptionallyQuoted(source_relation);
	auto q_col = KeywordHelper::WriteOptionallyQuoted(sample_id_col);

	auto probe = conn.Query("SELECT " + q_col + " FROM " + q_src + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: sample_id column '%s' not found", fn_label, sample_id_col);
	}
	out.sample_id_type = probe->types[0];

	// Single scan: DISTINCT with NULLs first so we can reject up-front.
	auto distinct_result =
	    conn.Query("SELECT DISTINCT " + q_col + " FROM " + q_src + " ORDER BY " + q_col + " NULLS FIRST");
	if (distinct_result->HasError()) {
		throw InvalidInputException("%s: failed to query sample_id values: %s", fn_label, distinct_result->GetError());
	}
	auto &materialized = distinct_result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto val = chunk->data[0].GetValue(i);
			if (val.IsNull()) {
				throw InvalidInputException("%s: NULL values in sample_id column '%s'", fn_label, sample_id_col);
			}
			out.sample_values.push_back(std::move(val));
		}
	}
	if (out.sample_values.empty()) {
		throw InvalidInputException("%s: sample_id column '%s' has no non-NULL values", fn_label, sample_id_col);
	}
	out.sample_id_col = sample_id_col;
}

void InitPerSampleGlobal(ClientContext &context, PerSampleGlobalState &gstate, idx_t num_samples,
                         idx_t max_threads_hint) {
	if (num_samples == 0) {
		gstate.max_threads = 1;
		return;
	}
	idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	idx_t cap = (max_threads_hint == 0) ? db_threads : max_threads_hint;
	gstate.max_threads = std::max<idx_t>(1, std::min<idx_t>(cap, num_samples));
}

bool ClaimNextSample(PerSampleGlobalState &gstate, idx_t num_samples, idx_t &out_idx) {
	auto claimed = gstate.next_sample_idx.fetch_add(1, std::memory_order_relaxed);
	if (claimed >= num_samples) {
		return false;
	}
	out_idx = claimed;
	return true;
}

} // namespace duckdb
