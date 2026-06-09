#include "sequence_table_reader.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/appender.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"

namespace duckdb {

SequenceTableSchema ValidateSequenceTableSchema(ClientContext &context, const std::string &table_name,
                                                bool allow_bigint) {
	auto info = GetTableOrViewColumns(context, table_name, "Sequence table");
	auto &col_names = info.names;
	auto &col_types = info.types;
	bool is_physical_table = info.is_physical_table;

	// Build name-to-index map (case-insensitive)
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	SequenceTableSchema schema;
	schema.is_physical_table = is_physical_table;

	// Check required columns
	auto check_column = [&](const string &col_name, const vector<LogicalTypeId> &allowed_types, const string &type_desc,
	                        bool required) -> bool {
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

	// Required columns: read_id (VARCHAR; or BIGINT/UUID if allow_bigint), sequence1 (VARCHAR)
	if (allow_bigint) {
		check_column("read_id", {LogicalTypeId::VARCHAR, LogicalTypeId::BIGINT, LogicalTypeId::UUID},
		             AllowedIdTypeList(), true);
	} else {
		check_column("read_id", {LogicalTypeId::VARCHAR}, "VARCHAR", true);
	}
	check_column("sequence1", {LogicalTypeId::VARCHAR}, "VARCHAR", true);

	// Record the read_id column's storage type so downstream readers can
	// stringify on ingress and the bind layer can mirror the type on output.
	{
		auto it = name_to_idx.find("read_id");
		schema.id_type = col_types[it->second];
	}

	// Optional columns: sequence2, qual1, qual2
	schema.has_sequence2 = check_column("sequence2", {LogicalTypeId::VARCHAR}, "VARCHAR", false);
	schema.has_qual1 = check_column("qual1", {LogicalTypeId::LIST}, "LIST", false);
	schema.has_qual2 = check_column("qual2", {LogicalTypeId::LIST}, "LIST", false);

	return schema;
}

std::vector<miint::AlignmentSubject> ReadSubjectTable(ClientContext &context, const std::string &table_name,
                                                      const SequenceTableSchema &schema) {
	std::vector<miint::AlignmentSubject> result;

	// Create a new connection to avoid deadlocking
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Query only required columns - try with sequence2 first to detect paired data
	std::string query = "SELECT read_id, sequence1, sequence2 FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		// If sequence2 doesn't exist, try without it
		query = "SELECT read_id, sequence1, NULL as sequence2 FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
		query_result = conn.Query(query);
		if (query_result->HasError()) {
			throw InvalidInputException("Failed to read from subject table '%s': %s", table_name,
			                            query_result->GetError());
		}
	}

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	idx_t row_number = 0;

	// Reusable buffers across chunks — populated by ExtractIdColumnAsStrings.
	std::vector<std::string> id_strings;
	std::vector<bool> id_nulls;

	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}

		ExtractIdColumnAsStrings(*chunk, /*col_idx=*/0, schema.id_type, id_strings, id_nulls);

		auto &seq1_vec = chunk->data[1];
		auto &seq2_vec = chunk->data[2];

		UnifiedVectorFormat seq1_data, seq2_data;
		seq1_vec.ToUnifiedFormat(chunk->size(), seq1_data);
		seq2_vec.ToUnifiedFormat(chunk->size(), seq2_data);

		auto sequences1 = UnifiedVectorFormat::GetData<string_t>(seq1_data);

		for (idx_t i = 0; i < chunk->size(); i++) {
			row_number++;

			auto seq1_idx = seq1_data.sel->get_index(i);
			auto seq2_idx = seq2_data.sel->get_index(i);

			// Skip rows with NULL read_id or sequence1
			if (id_nulls[i]) {
				continue;
			}
			if (!seq1_data.validity.RowIsValid(seq1_idx)) {
				continue;
			}

			// Check sequence2 is NULL (subjects can't be paired)
			if (seq2_data.validity.RowIsValid(seq2_idx)) {
				throw InvalidInputException(
				    "Subject table '%s' has non-NULL sequence2 at row %llu. Subjects cannot be paired-end.", table_name,
				    row_number);
			}

			miint::AlignmentSubject subject;
			subject.read_id = std::move(id_strings[i]);
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
static miint::QualScore ExtractQualScore(DataChunk &chunk, idx_t qual_col_idx, UnifiedVectorFormat &qual_data,
                                         idx_t row) {
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

// Helper to create a new sub-batch with reserved capacity.
static miint::SequenceRecordBatch MakeSubBatch(bool is_paired, idx_t capacity) {
	miint::SequenceRecordBatch batch(is_paired);
	if (capacity > 0) {
		batch.reserve(capacity);
	}
	return batch;
}

// Process a single DataChunk into a SequenceRecordBatch.
// Handles two-pass extraction (strings first, then quals) to avoid pointer corruption.
// temp_* vectors are caller-owned and reused across calls to avoid re-allocation.
// When sub_batches and sub_batch_size are provided, flushes full sub-batches as rows
// are appended — `output` acts as the current partial sub-batch.
static void ProcessSingleChunk(DataChunk &chunk, const SequenceTableSchema &schema, miint::SequenceRecordBatch &output,
                               std::vector<std::string> &temp_read_ids, std::vector<std::string> &temp_seq1,
                               std::vector<std::string> &temp_seq2,
                               std::vector<miint::SequenceRecordBatch> *sub_batches = nullptr,
                               idx_t sub_batch_size = 0) {
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

	// Pre-extract column 0 (read_id) as strings, dispatching on the captured
	// id_type. This is the single point where VARCHAR-vs-BIGINT ingress
	// branching lives — every downstream consumer sees a uniform string id.
	std::vector<std::string> chunk_id_strings;
	std::vector<bool> chunk_id_nulls;
	ExtractIdColumnAsStrings(chunk, read_id_col, schema.id_type, chunk_id_strings, chunk_id_nulls);

	// Prepare unified formats for the remaining columns
	UnifiedVectorFormat seq1_data;
	chunk.data[seq1_col].ToUnifiedFormat(chunk.size(), seq1_data);

	auto sequences1 = UnifiedVectorFormat::GetData<string_t>(seq1_data);

	UnifiedVectorFormat seq2_data, qual1_data, qual2_data;
	const string_t *sequences2 = nullptr;
	if (schema.has_sequence2) {
		chunk.data[seq2_col].ToUnifiedFormat(chunk.size(), seq2_data);
		sequences2 = UnifiedVectorFormat::GetData<string_t>(seq2_data);
	}
	if (schema.has_qual1) {
		chunk.data[qual1_col].ToUnifiedFormat(chunk.size(), qual1_data);
	}
	if (schema.has_qual2) {
		chunk.data[qual2_col].ToUnifiedFormat(chunk.size(), qual2_data);
	}

	// IMPORTANT: Extract ALL string data FIRST before calling ExtractQualScore.
	// ExtractQualScore calls ListVector::GetEntry() which may corrupt string pointers.
	// clear() preserves heap capacity — no re-allocation after first chunk.
	temp_read_ids.clear();
	temp_seq1.clear();
	temp_seq2.clear();

	for (idx_t i = 0; i < chunk.size(); i++) {
		auto seq1_idx = seq1_data.sel->get_index(i);

		if (chunk_id_nulls[i]) {
			continue;
		}
		if (!seq1_data.validity.RowIsValid(seq1_idx)) {
			continue;
		}

		temp_read_ids.push_back(std::move(chunk_id_strings[i]));
		temp_seq1.push_back(sequences1[seq1_idx].GetString());

		if (output.is_paired && schema.has_sequence2 && sequences2) {
			auto seq2_idx = seq2_data.sel->get_index(i);
			if (seq2_data.validity.RowIsValid(seq2_idx)) {
				temp_seq2.push_back(sequences2[seq2_idx].GetString());
			} else {
				temp_seq2.push_back("");
			}
		} else if (output.is_paired) {
			temp_seq2.push_back("");
		}
	}

	// Now process the extracted strings and quality scores
	idx_t batch_idx = 0;
	for (idx_t i = 0; i < chunk.size(); i++) {
		auto seq1_idx = seq1_data.sel->get_index(i);

		if (chunk_id_nulls[i]) {
			continue;
		}
		if (!seq1_data.validity.RowIsValid(seq1_idx)) {
			continue;
		}

		output.read_ids.push_back(std::move(temp_read_ids[batch_idx]));
		output.comments.push_back("");
		output.sequences1.push_back(std::move(temp_seq1[batch_idx]));

		if (schema.has_qual1) {
			output.quals1.push_back(ExtractQualScore(chunk, qual1_col, qual1_data, i));
		} else {
			output.quals1.push_back(miint::QualScore(""));
		}

		if (output.is_paired) {
			output.sequences2.push_back(std::move(temp_seq2[batch_idx]));

			if (schema.has_qual2) {
				output.quals2.push_back(ExtractQualScore(chunk, qual2_col, qual2_data, i));
			} else {
				output.quals2.push_back(miint::QualScore(""));
			}
		}

		batch_idx++;

		if (sub_batches && sub_batch_size > 0 && output.size() >= sub_batch_size) {
			sub_batches->push_back(std::move(output));
			output = MakeSubBatch(schema.has_sequence2, sub_batch_size);
		}
	}
}

// Helper to process all chunks from a query result into SequenceRecordBatch.
// Returns total number of rows processed.
static idx_t ProcessQueryResultChunks(QueryResult &result, const SequenceTableSchema &schema,
                                      miint::SequenceRecordBatch &output,
                                      std::vector<miint::SequenceRecordBatch> *sub_batches = nullptr,
                                      idx_t sub_batch_size = 0) {
	idx_t total_rows = 0;
	// Reusable temp vectors — allocated once, cleared per chunk
	std::vector<std::string> temp_read_ids, temp_seq1, temp_seq2;

	while (true) {
		auto chunk = result.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		total_rows += chunk->size();
		ProcessSingleChunk(*chunk, schema, output, temp_read_ids, temp_seq1, temp_seq2, sub_batches, sub_batch_size);
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
	                    KeywordHelper::WriteOptionallyQuoted(table_name) + " ORDER BY " + order_col + " LIMIT " +
	                    std::to_string(batch_size) + " OFFSET " + std::to_string(offset);

	auto query_result = conn.Query(query);

	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read from query table '%s': %s", table_name, query_result->GetError());
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

std::vector<std::string> ReadShardIds(ClientContext &context, const std::string &read_to_shard_table,
                                      const std::string &shard_name, const LogicalType &id_type) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	std::string query = "SELECT read_id FROM " + KeywordHelper::WriteOptionallyQuoted(read_to_shard_table) +
	                    " WHERE shard_name = " + KeywordHelper::WriteQuoted(shard_name, '\'') + " ORDER BY read_id";

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read IDs for shard '%s': %s", shard_name, query_result->GetError());
	}

	std::vector<std::string> ids;
	auto &materialized = query_result->Cast<MaterializedQueryResult>();

	if (!IsAllowedIdType(id_type)) {
		throw InternalException("ReadShardIds: id_type must be %s, got '%s'", AllowedIdTypeList(), id_type.ToString());
	}

	// Stringify each chunk through the shared id ingress dispatcher (VARCHAR /
	// BIGINT / UUID) and keep only the non-NULL ids, in order.
	std::vector<std::string> chunk_vals;
	std::vector<bool> chunk_nulls;
	while (true) {
		auto chunk = materialized.Fetch();
		if (!chunk || chunk->size() == 0) {
			break;
		}
		ExtractIdColumnAsStrings(*chunk, 0, id_type, chunk_vals, chunk_nulls);
		for (idx_t i = 0; i < chunk_vals.size(); i++) {
			if (!chunk_nulls[i]) {
				ids.push_back(std::move(chunk_vals[i]));
			}
		}
	}

	return ids;
}

void ReadBatchByIds(ClientContext &context, const std::string &query_table, const SequenceTableSchema &schema,
                    const std::vector<std::string> &ids, idx_t offset, idx_t count,
                    miint::SequenceRecordBatch &output) {
	// Clamp count to available IDs
	if (offset >= ids.size()) {
		return;
	}
	count = std::min(count, static_cast<idx_t>(ids.size()) - offset);
	if (count == 0) {
		return;
	}

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);

	// Declare the temp table with the same id_type as the query table so the
	// downstream JOIN type-checks naturally. The `ids` vector holds stringified
	// ids regardless of source type; for BIGINT we parse back through the codec
	// before appending. INVALID is rejected here — the project convention is
	// that any SequenceTableSchema reaching this layer has been through
	// ValidateSequenceTableSchema, which always resolves to VARCHAR or BIGINT.
	if (!IsAllowedIdType(schema.id_type)) {
		throw InternalException("ReadBatchByIds: schema.id_type must be %s, got '%s'", AllowedIdTypeList(),
		                        schema.id_type.ToString());
	}
	const LogicalType &id_type = schema.id_type;
	const std::string create_sql = "CREATE TEMPORARY TABLE _batch_ids (read_id " + id_type.ToString() + ")";
	auto create_result = conn.Query(create_sql);
	if (create_result->HasError()) {
		throw InvalidInputException("Failed to create temp table for batch IDs: %s", create_result->GetError());
	}

	{
		Appender appender(conn, "_batch_ids");
		for (idx_t i = offset; i < offset + count; i++) {
			if (id_type.id() == LogicalTypeId::BIGINT) {
				auto parsed = miint::ParseIdAsInt64(ids[i]);
				if (parsed.has_value()) {
					appender.AppendRow(Value::BIGINT(*parsed));
				} else {
					appender.AppendRow(Value(LogicalType::BIGINT));
				}
			} else if (id_type.id() == LogicalTypeId::UUID) {
				// ids[i] is a canonical UUID string produced by ReadShardIds. Parse
				// via the bool-checked helper and pass the INT128 to the hugeint
				// overload — Value::UUID(string) silently swallows parse failures,
				// so this fails loud on a malformed id, mirroring the BIGINT branch.
				appender.AppendRow(Value::UUID(ParseUuidOrThrow(ids[i])));
			} else {
				appender.AppendRow(Value(ids[i]));
			}
		}
		appender.Close();
	}

	// Join against query table using the temp table of exact IDs
	// No ORDER BY needed — alignment doesn't depend on order
	std::string query = "SELECT " + BuildSequenceColumnList(schema, "q.") + " FROM " +
	                    KeywordHelper::WriteOptionallyQuoted(query_table) +
	                    " q JOIN _batch_ids b ON q.read_id = b.read_id";

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read batch sequences: %s", query_result->GetError());
	}

	// Clear output and set paired flag
	output.clear();
	output.is_paired = schema.has_sequence2;

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	ProcessQueryResultChunks(materialized, schema, output);
}

QuerySequenceStream::QuerySequenceStream(ClientContext &context, const std::string &table_name,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size)
    : owned_conn_(make_uniq<Connection>(DatabaseInstance::GetDatabase(context))), conn_ptr_(owned_conn_.get()),
      schema_(schema), sub_batch_size_(sub_batch_size), partial_(schema.has_sequence2) {
	InitStream(table_name);
}

QuerySequenceStream::QuerySequenceStream(Connection &conn, const std::string &table_name,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size)
    : owned_conn_(nullptr), conn_ptr_(&conn), schema_(schema), sub_batch_size_(sub_batch_size),
      partial_(schema.has_sequence2) {
	InitStream(table_name);
}

void QuerySequenceStream::InitStream(const std::string &table_name) {
	partial_.reserve(sub_batch_size_);

	std::string query =
	    "SELECT " + BuildSequenceColumnList(schema_) + " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);

	stream_ = conn_ptr_->SendQuery(query);
	if (stream_->HasError()) {
		throw InvalidInputException("Failed to read from query table '%s': %s", table_name, stream_->GetError());
	}
}

miint::SequenceRecordBatch QuerySequenceStream::FetchSubBatch() {
	std::lock_guard<std::mutex> lock(mutex_);

	if (exhausted_ && partial_.empty()) {
		return miint::SequenceRecordBatch(schema_.has_sequence2);
	}

	// Fetch chunks from stream until we have a full sub-batch or stream is exhausted
	while (!exhausted_ && partial_.size() < sub_batch_size_) {
		auto chunk = stream_->Fetch();
		if (!chunk || chunk->size() == 0) {
			exhausted_ = true;
			break;
		}
		ProcessSingleChunk(*chunk, schema_, partial_, temp_read_ids_, temp_seq1_, temp_seq2_);
	}

	if (partial_.size() >= sub_batch_size_) {
		if (partial_.size() == sub_batch_size_) {
			// Exact fit — return as-is
			miint::SequenceRecordBatch result = std::move(partial_);
			partial_ = MakeSubBatch(schema_.has_sequence2, sub_batch_size_);
			return result;
		}
		// Chunk pushed partial_ past sub_batch_size_ — split: return exactly
		// sub_batch_size_ rows and keep the remainder for the next call.
		miint::SequenceRecordBatch result = partial_.SubRange(0, sub_batch_size_);
		miint::SequenceRecordBatch remainder = partial_.SubRange(sub_batch_size_, partial_.size() - sub_batch_size_);
		partial_ = std::move(remainder);
		return result;
	}

	// Stream exhausted — return whatever remains
	miint::SequenceRecordBatch result = std::move(partial_);
	partial_ = miint::SequenceRecordBatch(schema_.has_sequence2);
	return result;
}

LoadedSingleEndSequences LoadSingleEndSequences(Connection &conn, const std::string &table_name,
                                                const std::string &function_name, bool strict,
                                                const std::string &where_sql) {
	auto sql = "SELECT read_id, sequence1 FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	if (!where_sql.empty()) {
		sql += " WHERE " + where_sql;
	}
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("Failed to read table '%s': %s", table_name, result->GetError());
	}

	LoadedSingleEndSequences loaded;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	auto row_count = materialized.RowCount();
	loaded.labels.reserve(row_count);
	loaded.sequences.reserve(row_count);

	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto read_id_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			if (strict) {
				if (read_id_val.IsNull()) {
					throw InvalidInputException("%s: NULL read_id found in table '%s'. "
					                            "All rows must have a non-NULL read_id.",
					                            function_name, table_name);
				}
				if (seq_val.IsNull()) {
					throw InvalidInputException("%s: NULL sequence1 found for read_id '%s' in table '%s'",
					                            function_name, read_id_val.GetValue<std::string>(), table_name);
				}
				auto seq_str = seq_val.GetValue<std::string>();
				if (seq_str.empty()) {
					throw InvalidInputException("%s: empty sequence1 found for read_id '%s' in table '%s'",
					                            function_name, read_id_val.GetValue<std::string>(), table_name);
				}
				loaded.labels.push_back(read_id_val.GetValue<std::string>());
				loaded.sequences.push_back(std::move(seq_str));
			} else {
				if (read_id_val.IsNull() || seq_val.IsNull()) {
					continue;
				}
				auto seq_str = seq_val.GetValue<std::string>();
				if (seq_str.empty()) {
					continue;
				}
				loaded.labels.push_back(read_id_val.GetValue<std::string>());
				loaded.sequences.push_back(std::move(seq_str));
			}
		}
	}

	if (loaded.labels.empty()) {
		throw InvalidInputException("Table '%s' is empty (or contains only NULL/empty sequences)", table_name);
	}

	return loaded;
}

LoadedSingleEndSequences LoadSingleEndSequences(ClientContext &context, const std::string &table_name,
                                                const std::string &function_name, bool strict) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	return LoadSingleEndSequences(conn, table_name, function_name, strict, /*where_sql=*/"");
}

} // namespace duckdb
