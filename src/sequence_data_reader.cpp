#include "sequence_data_reader.hpp"
#include "catalog_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"

namespace duckdb {

SequenceDataSchema ValidateSequenceDataSchema(ClientContext &context, const std::string &table_name) {
	auto info = GetTableOrViewColumns(context, table_name, "Sequence data table");
	auto &col_names = info.names;
	auto &col_types = info.types;

	// Build case-insensitive name-to-index map
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	auto require_varchar = [&](const string &col_name) {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			throw BinderException("Sequence data table '%s' missing required column '%s'", table_name, col_name);
		}
		if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
			throw BinderException("Column '%s' in sequence data table '%s' must be VARCHAR", col_name, table_name);
		}
	};

	auto require_list_utinyint = [&](const string &col_name) {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			throw BinderException("Sequence data table '%s' missing required column '%s'", table_name, col_name);
		}
		auto &t = col_types[it->second];
		if (t.id() != LogicalTypeId::LIST || ListType::GetChildType(t).id() != LogicalTypeId::UTINYINT) {
			throw BinderException("Column '%s' in sequence data table '%s' must be LIST(UTINYINT)", col_name,
			                      table_name);
		}
	};

	auto has_varchar = [&](const string &col_name) -> bool {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			return false;
		}
		if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
			throw BinderException("Column '%s' in sequence data table '%s' must be VARCHAR", col_name, table_name);
		}
		return true;
	};

	auto has_list_utinyint = [&](const string &col_name) -> bool {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			return false;
		}
		auto &t = col_types[it->second];
		if (t.id() != LogicalTypeId::LIST || ListType::GetChildType(t).id() != LogicalTypeId::UTINYINT) {
			throw BinderException("Column '%s' in sequence data table '%s' must be LIST(UTINYINT)", col_name,
			                      table_name);
		}
		return true;
	};

	// Required columns
	require_varchar("read_id");
	require_varchar("sequence1");
	require_list_utinyint("qual1");

	// Optional columns (for paired-end)
	SequenceDataSchema schema;
	schema.has_sequence2 = has_varchar("sequence2");
	schema.has_qual2 = has_list_utinyint("qual2");

	return schema;
}

// Extract a LIST(UTINYINT) value as vector<uint8_t>.
// Returns empty vector for NULL values.
static std::vector<uint8_t> ExtractQualVector(DataChunk &chunk, idx_t col_idx, UnifiedVectorFormat &qual_data,
                                              idx_t row) {
	auto qual_row = qual_data.sel->get_index(row);
	if (!qual_data.validity.RowIsValid(qual_row)) {
		return {};
	}

	auto &qual_vec = chunk.data[col_idx];
	auto &qual_list = ListVector::GetEntry(qual_vec);
	auto qual_list_data = FlatVector::GetData<uint8_t>(qual_list);
	auto qual_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual_data);

	idx_t length = qual_entries[qual_row].length;
	idx_t offset = qual_entries[qual_row].offset;

	return std::vector<uint8_t>(qual_list_data + offset, qual_list_data + offset + length);
}

SequenceDataMap ReadSequenceDataTable(ClientContext &context, const std::string &table_name) {
	SequenceDataMap result;

	auto schema = ValidateSequenceDataSchema(context, table_name);

	// Separate connection to avoid deadlocking the current context
	auto conn = MakeReadOnlyHelperConnection(context);

	// Build query based on which columns exist
	std::string columns = "read_id, sequence1, qual1";
	if (schema.has_sequence2) {
		columns += ", sequence2";
	}
	if (schema.has_qual2) {
		columns += ", qual2";
	}
	std::string query = "SELECT " + columns + " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from sequence data table '%s': %s", table_name,
		                            query_result->GetError());
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	idx_t row_number = 0;

	// Column indices in our SELECT — dynamic based on schema
	constexpr idx_t COL_READ_ID = 0;
	constexpr idx_t COL_SEQ1 = 1;
	constexpr idx_t COL_QUAL1 = 2;
	idx_t col_seq2 = DConstants::INVALID_INDEX;
	idx_t col_qual2 = DConstants::INVALID_INDEX;
	idx_t next_col = 3;
	if (schema.has_sequence2) {
		col_seq2 = next_col++;
	}
	if (schema.has_qual2) {
		col_qual2 = next_col++;
	}

	// Temp vectors for two-pass extraction (strings first, then lists).
	// ListVector::GetEntry() can corrupt string pointers, so all strings
	// must be copied to std::string before any list extraction.
	std::vector<std::string> temp_read_ids;
	std::vector<std::string> temp_seq1;
	std::vector<std::string> temp_seq2;
	std::vector<bool> temp_seq2_valid;

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		UnifiedVectorFormat read_id_data, seq1_data, qual1_data;
		chunk->data[COL_READ_ID].ToUnifiedFormat(chunk->size(), read_id_data);
		chunk->data[COL_SEQ1].ToUnifiedFormat(chunk->size(), seq1_data);
		chunk->data[COL_QUAL1].ToUnifiedFormat(chunk->size(), qual1_data);

		UnifiedVectorFormat seq2_data, qual2_data;
		const string_t *seq2_ptr = nullptr;
		if (schema.has_sequence2) {
			chunk->data[col_seq2].ToUnifiedFormat(chunk->size(), seq2_data);
			seq2_ptr = UnifiedVectorFormat::GetData<string_t>(seq2_data);
		}
		if (schema.has_qual2) {
			chunk->data[col_qual2].ToUnifiedFormat(chunk->size(), qual2_data);
		}

		auto read_id_ptr = UnifiedVectorFormat::GetData<string_t>(read_id_data);
		auto seq1_ptr = UnifiedVectorFormat::GetData<string_t>(seq1_data);

		// Pass 1: extract all strings
		temp_read_ids.clear();
		temp_seq1.clear();
		temp_seq2.clear();
		temp_seq2_valid.clear();

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;
			auto rid_idx = read_id_data.sel->get_index(i);
			auto s1_idx = seq1_data.sel->get_index(i);

			if (!read_id_data.validity.RowIsValid(rid_idx)) {
				throw InvalidInputException("NULL read_id at row %llu in sequence data table '%s'", row_number,
				                            table_name);
			}
			if (!seq1_data.validity.RowIsValid(s1_idx)) {
				throw InvalidInputException("NULL sequence1 at row %llu in sequence data table '%s'", row_number,
				                            table_name);
			}

			temp_read_ids.push_back(read_id_ptr[rid_idx].GetString());
			temp_seq1.push_back(seq1_ptr[s1_idx].GetString());

			if (schema.has_sequence2) {
				auto s2_idx = seq2_data.sel->get_index(i);
				if (seq2_data.validity.RowIsValid(s2_idx)) {
					temp_seq2.push_back(seq2_ptr[s2_idx].GetString());
					temp_seq2_valid.push_back(true);
				} else {
					temp_seq2.emplace_back();
					temp_seq2_valid.push_back(false);
				}
			} else {
				temp_seq2.emplace_back();
				temp_seq2_valid.push_back(false);
			}
		}

		// Pass 2: extract list columns (safe now that strings are copied)
		for (idx_t i = 0; i < chunk->size(); i++) {
			SequenceDataEntry entry;
			entry.sequence1 = std::move(temp_seq1[i]);
			entry.qual1 = ExtractQualVector(*chunk, COL_QUAL1, qual1_data, i);
			entry.sequence2 = std::move(temp_seq2[i]);
			if (temp_seq2_valid[i] && schema.has_qual2) {
				entry.qual2 = ExtractQualVector(*chunk, col_qual2, qual2_data, i);
			}

			auto &read_id = temp_read_ids[i];
			auto [it, inserted] = result.emplace(read_id, std::move(entry));
			if (!inserted) {
				throw InvalidInputException("Duplicate read_id '%s' in sequence data table '%s'", read_id, table_name);
			}
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Sequence data table '%s' is empty", table_name);
	}

	return result;
}

} // namespace duckdb
