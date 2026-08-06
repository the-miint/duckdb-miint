#include "absquant_readers.hpp"

#include "catalog_utils.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/query_result.hpp"

#include <cmath>

namespace duckdb::absquant_internal {

namespace {

// "(feature_id VARCHAR, syndna_indiv_ng_ul DOUBLE)" -- the column contract as it
// appears in the probe's error message. Rendered from the same strings the
// projection is built from, so the two can never drift.
std::string DescribeContract(const char *key_column, const std::vector<const char *> &value_columns) {
	std::string out = "(" + std::string(key_column) + " VARCHAR";
	for (const auto *column : value_columns) {
		out += ", " + std::string(column) + " DOUBLE";
	}
	return out + ")";
}

} // namespace

KeyedColumns ReadKeyedColumns(ClientContext &context, const std::string &table_name, const char *key_column,
                              const std::vector<const char *> &value_columns, const char *entity,
                              const std::string &caller_name) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);

	std::string projection = std::string(key_column) + "::VARCHAR";
	for (const auto *column : value_columns) {
		projection += ", " + std::string(column) + "::DOUBLE";
	}

	// Schema probe first, so a missing column or an impossible cast surfaces
	// before the whole relation is materialized (mirrors ReadFeatureTable).
	auto probe = conn.Query("SELECT " + projection + " FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: %s '%s' must expose %s: %s", caller_name, entity, table_name,
		                            DescribeContract(key_column, value_columns), probe->GetError());
	}
	auto result = conn.Query("SELECT " + projection + " FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read %s '%s': %s", caller_name, entity, table_name,
		                            result->GetError());
	}

	KeyedColumns out;
	std::vector<std::vector<double> *> column_data;
	column_data.reserve(value_columns.size());
	for (const auto *column : value_columns) {
		// A repeated column name would collapse two entries into one and leave
		// the caller silently reading the wrong column's values. The names are
		// compile-time literals, so this can only ever be a caller bug.
		const auto inserted = out.values.emplace(column, std::vector<double> {});
		if (!inserted.second) {
			throw InternalException("ReadKeyedColumns: duplicate value column '%s'", column);
		}
		column_data.push_back(&inserted.first->second);
	}

	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t n = chunk->size();
		if (n == 0) {
			break;
		}
		UnifiedVectorFormat key_u;
		chunk->data[0].ToUnifiedFormat(n, key_u);
		auto key_data = UnifiedVectorFormat::GetData<string_t>(key_u);

		std::vector<UnifiedVectorFormat> val_u(value_columns.size());
		std::vector<const double *> val_data(value_columns.size());
		for (size_t c = 0; c < value_columns.size(); ++c) {
			chunk->data[c + 1].ToUnifiedFormat(n, val_u[c]);
			val_data[c] = UnifiedVectorFormat::GetData<double>(val_u[c]);
		}

		for (idx_t i = 0; i < n; ++i) {
			const auto ki = key_u.sel->get_index(i);
			if (!key_u.validity.RowIsValid(ki)) {
				throw InvalidInputException("%s: %s '%s' has a row with a NULL %s", caller_name, entity, table_name,
				                            key_column);
			}
			out.keys.push_back(key_data[ki].GetString());
			for (size_t c = 0; c < value_columns.size(); ++c) {
				const auto vi = val_u[c].sel->get_index(i);
				column_data[c]->push_back(val_u[c].validity.RowIsValid(vi) ? val_data[c][vi] : std::nan(""));
			}
		}
	}
	return out;
}

std::vector<LongFormRow> ReadLongFormValues(ClientContext &context, const std::string &table_name,
                                            const char *value_column, const char *entity,
                                            const std::string &caller_name) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const std::string projection = "sample_id::VARCHAR, feature_id::VARCHAR, " + std::string(value_column) + "::DOUBLE";

	auto probe = conn.Query("SELECT " + projection + " FROM " + qname + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: %s '%s' must expose (sample_id VARCHAR, feature_id VARCHAR, %s DOUBLE): %s",
		                            caller_name, entity, table_name, value_column, probe->GetError());
	}
	auto result = conn.Query("SELECT " + projection + " FROM " + qname);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read %s '%s': %s", caller_name, entity, table_name,
		                            result->GetError());
	}

	std::vector<LongFormRow> rows;
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
			if (!sid_u.validity.RowIsValid(si) || !fid_u.validity.RowIsValid(fi)) {
				continue; // an unidentified cell, dropped as ReadFeatureTable does
			}
			// No zero-drop, and no NaN-drop either: see the header. Both would
			// turn a cell the core must reject into one it never hears about.
			rows.push_back({sid_data[si].GetString(), fid_data[fi].GetString(),
			                val_u.validity.RowIsValid(vi) ? val_data[vi] : std::nan("")});
		}
	}
	return rows;
}

} // namespace duckdb::absquant_internal
