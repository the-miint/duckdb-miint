#include "sequence_table_reader.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"

namespace duckdb {

// Helper to get column names and types from either a table or view
// Returns true if it's a physical table (has rowid), false if it's a view
static bool GetTableOrViewColumns(ClientContext &context, const std::string &table_name,
                                  vector<string> &out_names, vector<LogicalType> &out_types) {
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);

	if (!entry) {
		throw BinderException("Table or view '%s' does not exist", table_name);
	}

	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &col = columns.GetColumn(LogicalIndex(i));
			out_names.push_back(col.Name());
			out_types.push_back(col.Type());
		}
		return true; // Physical table
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		out_names = view.names;
		out_types = view.types;
		return false; // View
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}
}

SequenceTableSchema ValidateSequenceTableSchema(ClientContext &context, const std::string &table_name,
                                                bool allow_paired) {
	vector<string> col_names;
	vector<LogicalType> col_types;
	bool is_physical_table = GetTableOrViewColumns(context, table_name, col_names, col_types);

	// Build name-to-index map (case-insensitive)
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	SequenceTableSchema schema;
	schema.is_physical_table = is_physical_table;

	// Check required columns
	auto check_column = [&](const string &col_name, const vector<LogicalTypeId> &allowed_types,
	                        const string &type_desc, bool required) -> bool {
		auto it = name_to_idx.find(col_name);
		if (it == name_to_idx.end()) {
			if (required) {
				throw BinderException("Sequence table '%s' missing required column '%s'", table_name, col_name);
			}
			return false;
		}
		auto &col_type = col_types[it->second];
		bool valid = false;
		for (auto &allowed : allowed_types) {
			if (col_type.id() == allowed) {
				valid = true;
				break;
			}
		}
		if (!valid) {
			throw BinderException("Column '%s' in table '%s' must be %s", col_name, table_name, type_desc);
		}
		return true;
	};

	// Required columns: read_id, sequence1
	check_column("read_id", {LogicalTypeId::VARCHAR}, "VARCHAR", true);
	check_column("sequence1", {LogicalTypeId::VARCHAR}, "VARCHAR", true);

	// Optional columns: sequence2, qual1, qual2
	schema.has_sequence2 = check_column("sequence2", {LogicalTypeId::VARCHAR}, "VARCHAR", false);
	schema.has_qual1 = check_column("qual1", {LogicalTypeId::LIST}, "LIST", false);
	schema.has_qual2 = check_column("qual2", {LogicalTypeId::LIST}, "LIST", false);

	if (schema.has_sequence2 && !allow_paired) {
		throw BinderException("Subject table '%s' has sequence2 column but subjects cannot be paired-end", table_name);
	}

	return schema;
}

std::vector<miint::AlignmentSubject> ReadSubjectTable(ClientContext &context, const std::string &table_name) {
	std::vector<miint::AlignmentSubject> result;

	// Create a new connection to avoid deadlocking
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Query only required columns - try with sequence2 first to detect paired data
	std::string query = "SELECT read_id, sequence1, sequence2 FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		// If sequence2 doesn't exist, try without it
		query = "SELECT read_id, sequence1, NULL as sequence2 FROM " +
		        KeywordHelper::WriteOptionallyQuoted(table_name);
		query_result = conn.Query(query);
		if (query_result->HasError()) {
			throw InvalidInputException("Failed to read from subject table '%s': %s", table_name,
			                            query_result->GetError());
		}
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	idx_t row_number = 0;

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		auto &read_id_vec = chunk->data[0];
		auto &seq1_vec = chunk->data[1];
		auto &seq2_vec = chunk->data[2];

		UnifiedVectorFormat read_id_data, seq1_data, seq2_data;
		read_id_vec.ToUnifiedFormat(chunk->size(), read_id_data);
		seq1_vec.ToUnifiedFormat(chunk->size(), seq1_data);
		seq2_vec.ToUnifiedFormat(chunk->size(), seq2_data);

		auto read_ids = UnifiedVectorFormat::GetData<string_t>(read_id_data);
		auto sequences1 = UnifiedVectorFormat::GetData<string_t>(seq1_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto rid_idx = read_id_data.sel->get_index(i);
			auto seq1_idx = seq1_data.sel->get_index(i);
			auto seq2_idx = seq2_data.sel->get_index(i);

			// Skip rows with NULL read_id or sequence1
			if (!read_id_data.validity.RowIsValid(rid_idx)) {
				continue;
			}
			if (!seq1_data.validity.RowIsValid(seq1_idx)) {
				continue;
			}

			// Check sequence2 is NULL (subjects can't be paired)
			if (seq2_data.validity.RowIsValid(seq2_idx)) {
				throw InvalidInputException(
				    "Subject table '%s' has non-NULL sequence2 at row %llu. Subjects cannot be paired-end.",
				    table_name, row_number);
			}

			miint::AlignmentSubject subject;
			subject.read_id = read_ids[rid_idx].GetString();
			subject.sequence = sequences1[seq1_idx].GetString();

			result.push_back(std::move(subject));
		}
	}

	if (result.empty()) {
		throw InvalidInputException("Subject table '%s' is empty", table_name);
	}

	return result;
}

// Helper to extract quality scores from a LIST<UTINYINT> vector
static miint::QualScore ExtractQualScore(DataChunk &chunk, idx_t qual_col_idx, UnifiedVectorFormat &qual_data, idx_t row) {
	auto qual_row = qual_data.sel->get_index(row);

	if (!qual_data.validity.RowIsValid(qual_row)) {
		// Return empty quality scores if NULL
		return miint::QualScore("");
	}

	auto &qual_vec = chunk.data[qual_col_idx];
	auto &qual_list = ListVector::GetEntry(qual_vec);
	auto qual_list_data = FlatVector::GetData<uint8_t>(qual_list);
	auto qual_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual_data);

	idx_t qual_length = qual_entries[qual_row].length;
	idx_t qual_offset = qual_entries[qual_row].offset;

	// Build vector of quality scores
	std::vector<uint8_t> qual_vec_data(qual_list_data + qual_offset, qual_list_data + qual_offset + qual_length);

	// Construct QualScore from uint8 vector (uses offset 33 by default)
	return miint::QualScore(qual_vec_data);
}

// Helper to process query result chunks into SequenceRecordBatch.
// Handles two-pass extraction (strings first, then quals) to avoid pointer corruption.
// Returns total number of rows processed.
static idx_t ProcessQueryResultChunks(MaterializedQueryResult &materialized, const SequenceTableSchema &schema,
                                       miint::SequenceRecordBatch &output) {
	idx_t total_rows = 0;

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		total_rows += chunk->size();

		// Column indices based on schema
		idx_t read_id_col = 0;
		idx_t seq1_col = 1;
		idx_t seq2_col = schema.has_sequence2 ? 2 : DConstants::INVALID_INDEX;
		idx_t qual1_col = DConstants::INVALID_INDEX;
		idx_t qual2_col = DConstants::INVALID_INDEX;

		idx_t next_col = 2;
		if (schema.has_sequence2) {
			next_col = 3;
		}
		if (schema.has_qual1) {
			qual1_col = next_col++;
		}
		if (schema.has_qual2) {
			qual2_col = next_col++;
		}

		// Prepare unified formats
		UnifiedVectorFormat read_id_data, seq1_data;
		chunk->data[read_id_col].ToUnifiedFormat(chunk->size(), read_id_data);
		chunk->data[seq1_col].ToUnifiedFormat(chunk->size(), seq1_data);

		auto read_ids = UnifiedVectorFormat::GetData<string_t>(read_id_data);
		auto sequences1 = UnifiedVectorFormat::GetData<string_t>(seq1_data);

		UnifiedVectorFormat seq2_data, qual1_data, qual2_data;
		const string_t *sequences2 = nullptr;
		if (schema.has_sequence2) {
			chunk->data[seq2_col].ToUnifiedFormat(chunk->size(), seq2_data);
			sequences2 = UnifiedVectorFormat::GetData<string_t>(seq2_data);
		}
		if (schema.has_qual1) {
			chunk->data[qual1_col].ToUnifiedFormat(chunk->size(), qual1_data);
		}
		if (schema.has_qual2) {
			chunk->data[qual2_col].ToUnifiedFormat(chunk->size(), qual2_data);
		}

		// IMPORTANT: Extract ALL string data FIRST before calling ExtractQualScore
		// ExtractQualScore calls ListVector::GetEntry() which may corrupt string pointers
		std::vector<std::string> batch_read_ids;
		std::vector<std::string> batch_seq1;
		std::vector<std::string> batch_seq2;

		for (idx_t i = 0; i < chunk->size(); i++) {
			auto rid_idx = read_id_data.sel->get_index(i);
			auto seq1_idx = seq1_data.sel->get_index(i);

			// Skip rows with NULL read_id or sequence1
			if (!read_id_data.validity.RowIsValid(rid_idx)) {
				continue;
			}
			if (!seq1_data.validity.RowIsValid(seq1_idx)) {
				continue;
			}

			// Extract strings NOW before any LIST operations
			batch_read_ids.push_back(read_ids[rid_idx].GetString());
			batch_seq1.push_back(sequences1[seq1_idx].GetString());

			if (output.is_paired && schema.has_sequence2 && sequences2) {
				auto seq2_idx = seq2_data.sel->get_index(i);
				if (seq2_data.validity.RowIsValid(seq2_idx)) {
					batch_seq2.push_back(sequences2[seq2_idx].GetString());
				} else {
					batch_seq2.push_back("");
				}
			} else if (output.is_paired) {
				batch_seq2.push_back("");
			}
		}

		// Now process the extracted strings and quality scores
		idx_t batch_idx = 0;
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto rid_idx = read_id_data.sel->get_index(i);
			auto seq1_idx = seq1_data.sel->get_index(i);

			// Skip rows with NULL read_id or sequence1 (same logic as above)
			if (!read_id_data.validity.RowIsValid(rid_idx)) {
				continue;
			}
			if (!seq1_data.validity.RowIsValid(seq1_idx)) {
				continue;
			}

			output.read_ids.push_back(std::move(batch_read_ids[batch_idx]));
			output.comments.push_back("");
			output.sequences1.push_back(std::move(batch_seq1[batch_idx]));

			// Extract qual1 - NOW it's safe to call ExtractQualScore
			if (schema.has_qual1) {
				output.quals1.push_back(ExtractQualScore(*chunk, qual1_col, qual1_data, i));
			} else {
				output.quals1.push_back(miint::QualScore(""));
			}

			// Handle sequence2 and qual2 for paired reads
			if (output.is_paired) {
				output.sequences2.push_back(std::move(batch_seq2[batch_idx]));

				if (schema.has_qual2) {
					output.quals2.push_back(ExtractQualScore(*chunk, qual2_col, qual2_data, i));
				} else {
					output.quals2.push_back(miint::QualScore(""));
				}
			}

			batch_idx++;
		}
	}

	return total_rows;
}

// Helper to build column list for sequence queries based on schema
static std::string BuildSequenceColumnList(const SequenceTableSchema &schema, const std::string &prefix = "") {
	std::string columns = prefix + "read_id, " + prefix + "sequence1";
	if (schema.has_sequence2) {
		columns += ", " + prefix + "sequence2";
	}
	if (schema.has_qual1) {
		columns += ", " + prefix + "qual1";
	}
	if (schema.has_qual2) {
		columns += ", " + prefix + "qual2";
	}
	return columns;
}

bool ReadQueryBatch(ClientContext &context, const std::string &table_name, const SequenceTableSchema &schema,
                    idx_t batch_size, idx_t &offset, miint::SequenceRecordBatch &output) {
	// Create a new connection to avoid deadlocking
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Build query with ORDER BY for deterministic pagination
	// Use rowid for physical tables (fast), read_id for views
	std::string order_col = schema.is_physical_table ? "rowid" : "read_id";
	std::string query = "SELECT " + BuildSequenceColumnList(schema) + " FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(table_name) +
	                    " ORDER BY " + order_col +
	                    " LIMIT " + std::to_string(batch_size) +
	                    " OFFSET " + std::to_string(offset);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from query table '%s': %s", table_name,
		                            query_result->GetError());
	}

	// Clear output and set paired flag
	output.clear();
	output.is_paired = schema.has_sequence2;

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	idx_t total_rows = ProcessQueryResultChunks(materialized, schema, output);

	// Update offset for next batch
	offset += total_rows;

	// Return true if we got a full batch (more rows may exist)
	return total_rows == batch_size;
}

bool ReadShardQueryBatch(ClientContext &context, const std::string &query_table,
                         const std::string &read_to_shard_table, const std::string &shard_name,
                         const SequenceTableSchema &schema, idx_t batch_size, idx_t &offset,
                         miint::SequenceRecordBatch &output) {
	// Create a new connection to avoid deadlocking
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Build query with JOIN to read_to_shard table, filtering by shard_name
	// Use "q." prefix for query table columns
	// ORDER BY is required for deterministic LIMIT/OFFSET pagination
	// Use rowid for physical tables (fast), read_id for views
	std::string order_col = schema.is_physical_table ? "q.rowid" : "q.read_id";
	std::string query = "SELECT " + BuildSequenceColumnList(schema, "q.") + " FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(query_table) + " q JOIN " +
	                    KeywordHelper::WriteOptionallyQuoted(read_to_shard_table) + " r " +
	                    "ON q.read_id = r.read_id WHERE r.shard_name = " +
	                    KeywordHelper::WriteQuoted(shard_name, '\'') +
	                    " ORDER BY " + order_col +
	                    " LIMIT " + std::to_string(batch_size) +
	                    " OFFSET " + std::to_string(offset);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read queries for shard '%s': %s", shard_name,
		                            query_result->GetError());
	}

	// Clear output and set paired flag
	output.clear();
	output.is_paired = schema.has_sequence2;

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	idx_t total_rows = ProcessQueryResultChunks(materialized, schema, output);

	// Update offset for next batch
	offset += total_rows;

	// Return true if we got a full batch (more rows may exist)
	return total_rows == batch_size;
}

} // namespace duckdb
