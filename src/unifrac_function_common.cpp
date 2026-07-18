#include "unifrac_function_common.hpp"

#include <algorithm>
#include <cmath>
#include <unordered_map>

#include "catalog_utils.hpp"
#include "unifrac_dense_distance.hpp"
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

DenseDistanceMatrix ReadDistanceTable(ClientContext &context, const std::string &table_name,
                                      const std::string &caller_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	// Schema probe via LIMIT 0 — surfaces missing columns or unsafe casts as a
	// binder-time error before we materialize the full relation (mirrors
	// ReadFeatureTable).
	auto probe = conn.Query("SELECT sample_a::VARCHAR, sample_b::VARCHAR, distance::DOUBLE FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: distance-table '%s' must expose (sample_a, sample_b, distance DOUBLE): %s",
		                            caller_name, table_name, probe->GetError());
	}

	// Capture sample_a's and sample_b's original SQL types (before the ::VARCHAR
	// cast) so the caller can mirror BIGINT/UUID ids onto its output. The probe
	// above already proved both columns exist and are VARCHAR-castable.
	LogicalType sample_a_type = LogicalType::VARCHAR;
	LogicalType sample_b_type = LogicalType::VARCHAR;
	{
		auto cols = GetTableOrViewColumns(context, table_name, "distance-table");
		for (idx_t i = 0; i < cols.names.size(); ++i) {
			const auto lname = StringUtil::Lower(cols.names[i]);
			if (lname == "sample_a") {
				sample_a_type = cols.types[i];
			} else if (lname == "sample_b") {
				sample_b_type = cols.types[i];
			}
		}
	}
	// sample_a and sample_b are merged into ONE id dictionary emitted under a
	// single output type. If they resolve to different output types (e.g. a
	// BIGINT sample_a but a VARCHAR sample_b), the mirror would emit sample_b's
	// ids under sample_a's type — a silent wrong-namespace merge if they happen
	// to parse, or a confusing deferred EmitIdCell failure if they don't. Reject
	// at bind instead (fail loud).
	if (ResolveSampleIdOutputType(sample_a_type) != ResolveSampleIdOutputType(sample_b_type)) {
		throw BinderException("%s: distance-table '%s' has mismatched sample_a (%s) and sample_b (%s) id types; "
		                      "both must map to the same output type (BIGINT, UUID, or VARCHAR)",
		                      caller_name, table_name, sample_a_type.ToString(), sample_b_type.ToString());
	}

	auto result = conn.Query("SELECT sample_a::VARCHAR, sample_b::VARCHAR, distance::DOUBLE FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read distance-table '%s': %s", caller_name, table_name,
		                            result->GetError());
	}

	// Collect surviving (a, b, distance) triples and the raw id list. NULL
	// sample ids or NULL/NaN distances are dropped ("not provided"); an unfilled
	// pair is caught downstream by the completeness check with a named message.
	struct RawTriple {
		std::string a;
		std::string b;
		double distance;
	};
	std::vector<RawTriple> triples;
	std::vector<std::string> ids_raw;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat a_u, b_u, d_u;
		chunk->data[0].ToUnifiedFormat(rn, a_u);
		chunk->data[1].ToUnifiedFormat(rn, b_u);
		chunk->data[2].ToUnifiedFormat(rn, d_u);
		auto a_data = UnifiedVectorFormat::GetData<string_t>(a_u);
		auto b_data = UnifiedVectorFormat::GetData<string_t>(b_u);
		auto d_data = UnifiedVectorFormat::GetData<double>(d_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ai = a_u.sel->get_index(i);
			const auto bi = b_u.sel->get_index(i);
			const auto di = d_u.sel->get_index(i);
			// A NULL sample id has no identity — skip the row entirely.
			if (!a_u.validity.RowIsValid(ai) || !b_u.validity.RowIsValid(bi)) {
				continue;
			}
			std::string a = a_data[ai].GetString();
			std::string b = b_data[bi].GetString();
			// Record both ids BEFORE gating on the distance: a sample that appears
			// in the table must enter the dictionary even if this particular
			// distance is missing, so a sample whose every distance is NULL/NaN
			// (e.g. a Bray-Curtis all-zero sample, NaN against everything)
			// surfaces as a completeness error rather than silently vanishing
			// from the result.
			ids_raw.push_back(a);
			ids_raw.push_back(b);
			if (!d_u.validity.RowIsValid(di)) {
				continue; // NULL distance → "not provided"
			}
			const double dv = d_data[di];
			if (std::isnan(dv)) {
				continue; // NaN distance → "not provided"
			}
			triples.push_back({std::move(a), std::move(b), dv});
		}
	}

	// Stable dictionary: distinct ids from both columns, sorted lexicographically
	// (matches build_dictionary / UnifracSupportBiomView::FromCoo). ids_raw is
	// dead after this move.
	std::vector<std::string> sample_ids = std::move(ids_raw);
	std::sort(sample_ids.begin(), sample_ids.end());
	sample_ids.erase(std::unique(sample_ids.begin(), sample_ids.end()), sample_ids.end());
	const auto n = static_cast<uint32_t>(sample_ids.size());
	if (n < 2) {
		throw InvalidInputException(
		    "%s: distance-table '%s' has %u distinct sample(s) after dropping NULL rows; at least 2 are required",
		    caller_name, table_name, n);
	}
	std::unordered_map<std::string, uint32_t> index;
	index.reserve(sample_ids.size());
	for (uint32_t i = 0; i < n; ++i) {
		index.emplace(sample_ids[i], i);
	}

	std::vector<miint::unifrac::DistanceEntry> entries;
	entries.reserve(triples.size());
	for (const auto &t : triples) {
		entries.push_back({index.at(t.a), index.at(t.b), t.distance});
	}

	DenseDistanceMatrix out;
	try {
		out.matrix = miint::unifrac::BuildDenseDistanceMatrix(entries, n, sample_ids);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", caller_name, e.what());
	}
	out.sample_ids = std::move(sample_ids);
	out.n_samples = n;
	out.sample_id_type = ResolveSampleIdOutputType(sample_a_type);
	return out;
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
