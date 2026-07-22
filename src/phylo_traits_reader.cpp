#include "phylo_traits_reader.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"
#include <map>
#include <string>
#include <unordered_map>

namespace duckdb {

// Validate that the traits table/view has the required long-form columns.
void ValidateTraitsTableSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Traits table");
	for (const char *required : {"name", "trait", "value"}) {
		if (!HasColumn(info, required)) {
			throw BinderException("Traits table '%s' missing required column '%s'", table_name, required);
		}
	}
}

// Read long-form traits into trait -> (tip name -> value). Uses a separate
// Connection (see docs/internals/reading-tables-views.md). `name` is cast to
// VARCHAR so UUID/BIGINT/VARCHAR tip keys all match tip labels by canonical
// text. A std::map keeps trait order stable (deterministic output). NULLs and
// duplicate (name, trait) pairs are hard errors.
std::map<std::string, std::unordered_map<std::string, double>> ReadContinuousTraits(ClientContext &context,
                                                                                    const std::string &table_name) {
	std::map<std::string, std::unordered_map<std::string, double>> traits;

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	std::string query =
	    "SELECT name::VARCHAR, trait::VARCHAR, value::DOUBLE FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from traits table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		UnifiedVectorFormat name_data, trait_data, value_data;
		chunk->data[0].ToUnifiedFormat(chunk->size(), name_data);
		chunk->data[1].ToUnifiedFormat(chunk->size(), trait_data);
		chunk->data[2].ToUnifiedFormat(chunk->size(), value_data);
		auto names = UnifiedVectorFormat::GetData<string_t>(name_data);
		auto trait_names = UnifiedVectorFormat::GetData<string_t>(trait_data);
		auto values = UnifiedVectorFormat::GetData<double>(value_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			auto ni = name_data.sel->get_index(i);
			auto ti = trait_data.sel->get_index(i);
			auto vi = value_data.sel->get_index(i);

			if (!name_data.validity.RowIsValid(ni)) {
				throw InvalidInputException("NULL name in traits table '%s'", table_name);
			}
			if (!trait_data.validity.RowIsValid(ti)) {
				throw InvalidInputException("NULL trait in traits table '%s'", table_name);
			}
			std::string nm = names[ni].GetString();
			std::string tr = trait_names[ti].GetString();
			if (!value_data.validity.RowIsValid(vi)) {
				throw InvalidInputException("NULL value for tip '%s' trait '%s' in traits table '%s'", nm, tr,
				                            table_name);
			}

			auto [it, inserted] = traits[tr].emplace(std::move(nm), values[vi]);
			if (!inserted) {
				throw InvalidInputException("Duplicate trait value for tip '%s' trait '%s' in traits table '%s'",
				                            it->first, tr, table_name);
			}
		}
	}

	if (traits.empty()) {
		throw InvalidInputException("Traits table '%s' is empty", table_name);
	}

	return traits;
}

// Read long-form traits into trait -> (tip name -> state string). Mirrors
// ReadContinuousTraits, but `value` is cast to VARCHAR so integer- or text-coded
// states both work; the caller builds the state alphabet. NULLs and duplicate
// (name, trait) pairs are hard errors.
std::map<std::string, std::unordered_map<std::string, std::string>> ReadDiscreteTraits(ClientContext &context,
                                                                                       const std::string &table_name) {
	std::map<std::string, std::unordered_map<std::string, std::string>> traits;

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	std::string query =
	    "SELECT name::VARCHAR, trait::VARCHAR, value::VARCHAR FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from traits table '%s': %s", table_name, query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		UnifiedVectorFormat name_data, trait_data, value_data;
		chunk->data[0].ToUnifiedFormat(chunk->size(), name_data);
		chunk->data[1].ToUnifiedFormat(chunk->size(), trait_data);
		chunk->data[2].ToUnifiedFormat(chunk->size(), value_data);
		auto names = UnifiedVectorFormat::GetData<string_t>(name_data);
		auto trait_names = UnifiedVectorFormat::GetData<string_t>(trait_data);
		auto values = UnifiedVectorFormat::GetData<string_t>(value_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			auto ni = name_data.sel->get_index(i);
			auto ti = trait_data.sel->get_index(i);
			auto vi = value_data.sel->get_index(i);

			if (!name_data.validity.RowIsValid(ni)) {
				throw InvalidInputException("NULL name in traits table '%s'", table_name);
			}
			if (!trait_data.validity.RowIsValid(ti)) {
				throw InvalidInputException("NULL trait in traits table '%s'", table_name);
			}
			std::string nm = names[ni].GetString();
			std::string tr = trait_names[ti].GetString();
			if (!value_data.validity.RowIsValid(vi)) {
				throw InvalidInputException("NULL value for tip '%s' trait '%s' in traits table '%s'", nm, tr,
				                            table_name);
			}

			auto [it, inserted] = traits[tr].emplace(std::move(nm), values[vi].GetString());
			if (!inserted) {
				throw InvalidInputException("Duplicate trait value for tip '%s' trait '%s' in traits table '%s'",
				                            it->first, tr, table_name);
			}
		}
	}

	if (traits.empty()) {
		throw InvalidInputException("Traits table '%s' is empty", table_name);
	}

	return traits;
}

} // namespace duckdb
