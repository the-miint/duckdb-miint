#include "unifrac_function_common.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
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
#include "duckdb/storage/buffer_manager.hpp"

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

DistanceRelationIds EnumerateDistanceIds(ClientContext &context, const std::string &table_name,
                                         const std::string &caller_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	// Schema probe via LIMIT 0 — surfaces missing columns or unsafe casts as a
	// binder-time error before any scan (mirrors ReadFeatureTable).
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

	auto res = conn.Query("SELECT id FROM (SELECT sample_a::VARCHAR AS id FROM " + qname +
	                      " UNION SELECT sample_b::VARCHAR FROM " + qname + ") WHERE id IS NOT NULL ORDER BY id");
	if (res->HasError()) {
		throw InvalidInputException("%s: failed to enumerate ids of distance-table '%s': %s", caller_name, table_name,
		                            res->GetError());
	}
	DistanceRelationIds out;
	out.sample_id_type = ResolveSampleIdOutputType(sample_a_type);
	auto &mat = res->Cast<MaterializedQueryResult>();
	while (auto chunk = mat.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat id_u;
		chunk->data[0].ToUnifiedFormat(rn, id_u);
		auto id_data = UnifiedVectorFormat::GetData<string_t>(id_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ii = id_u.sel->get_index(i);
			if (id_u.validity.RowIsValid(ii)) {
				out.sorted_ids.emplace_back(id_data[ii].GetString());
			}
		}
	}
	return out;
}

DenseDistanceMatrix ReadDistanceTable(ClientContext &context, const std::string &table_name,
                                      const std::string &caller_name) {
	// ── Pass 1: the id dictionary (schema probe + id-type check + sorted ids) ──
	// Learning N before reading any distance is what keeps this bounded. The
	// previous single-pass design had to park every row somewhere first (N was
	// unknown until the scan ended, so the matrix could not be allocated yet):
	// a fully materialized query result (~40 B/row) PLUS an entries vector
	// (16 B/row) — together ~28·N² bytes of pure intermediate, on top of the
	// matrix. That is what drove the 25k-sample EMP matrix into an OOM SIGKILL.
	// DuckDB does the DISTINCT + ORDER BY here, which it can spill; we hold only
	// the N-bounded dictionary.
	auto ids = EnumerateDistanceIds(context, table_name, caller_name);
	auto sample_ids = std::move(ids.sorted_ids);
	const auto n = static_cast<uint32_t>(sample_ids.size());
	if (n < 2) {
		throw InvalidInputException(
		    "%s: distance-table '%s' has %u distinct sample(s) after dropping NULL rows; at least 2 are required",
		    caller_name, table_name, n);
	}

	// ── Fail-loud size guard, before a single N² byte is allocated ──
	// A full dense analysis needs an N×N fp32 matrix plus a fill bitmap (~5·N²
	// bytes); with the in-place fsvd that is also the process peak. Refuse rather
	// than let a large N drive the process into an OOM SIGKILL. The budget is what
	// DuckDB's own memory_limit leaves unused: the matrix and skbb's workspace are
	// plain extension heap and therefore NOT buffer-manager tracked, so
	// memory_limit cannot stop them by itself — this check is what brings them
	// under the user's declared limit. Placing it after pass 1 (which holds only
	// the dictionary) means the estimate is compared against a budget that the
	// reader has not already eaten into.
	{
		const auto est_bytes = static_cast<idx_t>(5.0 * static_cast<double>(n) * static_cast<double>(n));
		auto &buffer_manager = BufferManager::GetBufferManager(context);
		const auto max_memory = buffer_manager.GetMaxMemory();
		const auto used_memory = buffer_manager.GetUsedMemory();
		const idx_t budget_bytes = max_memory > used_memory ? max_memory - used_memory : 0;
		if (est_bytes > budget_bytes) {
			throw InvalidInputException(
			    "%s: distance-table '%s' has %u distinct samples; a full dense analysis needs ~%s, exceeding the %s "
			    "left of the %s memory_limit. Raise memory_limit, or use progressive_pcoa_from_distances (or "
			    "progressive_pcoa_from_unifrac) for large sample counts.",
			    caller_name, table_name, n, StringUtil::BytesToHumanReadableString(est_bytes),
			    StringUtil::BytesToHumanReadableString(budget_bytes),
			    StringUtil::BytesToHumanReadableString(max_memory));
		}
	}

	// id → sorted index. The dictionary is already in the sorted order that is the
	// contract shared with build_dictionary / FromCoo, so cells can be written at
	// their final positions during the scan — no insertion-order remap step.
	std::unordered_map<std::string, uint32_t> index;
	index.reserve(n * 2);
	for (uint32_t k = 0; k < n; ++k) {
		index.emplace(sample_ids[k], k);
	}

	// ── Pass 2: stream the distances straight into the matrix ──
	// SendQuery (not Query) so chunks arrive lazily instead of the whole relation
	// being materialized: nothing per-row is retained.
	//
	// NULL sample ids or NULL/NaN distances are dropped ("not provided"); an
	// unfilled pair is caught by the builder's completeness check with a named
	// message. Pass 1 enumerated ids independently of the distance value, so a
	// sample whose every distance is NULL/NaN is in the dictionary and surfaces as
	// a completeness error rather than silently vanishing.
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	auto result = conn.SendQuery("SELECT sample_a::VARCHAR, sample_b::VARCHAR, distance::DOUBLE FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read distance-table '%s': %s", caller_name, table_name,
		                            result->GetError());
	}

	DenseDistanceMatrix out;
	try {
		miint::unifrac::DenseDistanceMatrixBuilder builder(n, sample_ids);
		while (auto chunk = result->Fetch()) {
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
				if (!d_u.validity.RowIsValid(di)) {
					continue; // NULL distance → "not provided"
				}
				const double dv = d_data[di];
				if (std::isnan(dv)) {
					continue; // NaN distance → "not provided"
				}
				// Both ids came from this same relation in pass 1, so a miss means
				// the relation did not return the same rows twice — a volatile or
				// concurrently-changing source (e.g. a view over random()/nextval,
				// or a file rewritten between passes). Fail loud: silently dropping
				// the row would report an incomplete matrix and blame the data.
				const auto ia = index.find(a_data[ai].GetString());
				const auto ib = index.find(b_data[bi].GetString());
				if (ia == index.end() || ib == index.end()) {
					throw InvalidInputException(
					    "%s: distance-table '%s' returned different rows on a second scan (sample id '%s' was not in "
					    "the enumerated id set); it must be stable across scans — materialize it into a table first",
					    caller_name, table_name, ia == index.end() ? a_data[ai].GetString() : b_data[bi].GetString());
				}
				builder.Add(ia->second, ib->second, dv);
			}
		}
		out.matrix = builder.Finish();
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("%s: %s", caller_name, e.what());
	}
	// A streaming result surfaces execution errors during Fetch, not at SendQuery.
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read distance-table '%s': %s", caller_name, table_name,
		                            result->GetError());
	}
	out.sample_ids = std::move(sample_ids);
	out.n_samples = n;
	out.sample_id_type = ids.sample_id_type;
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
