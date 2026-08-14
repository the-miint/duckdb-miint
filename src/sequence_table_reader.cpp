#include "sequence_table_reader.hpp"
#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/main/appender.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/types/uuid.hpp"

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
	auto conn = MakeReadOnlyHelperConnection(context);

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

std::string BuildShardedQueryReadsSelect(const std::string &query_table, const std::string &read_to_shard_table,
                                         const SequenceTableSchema &schema) {
	// No ORDER BY: alignment does not depend on read order. Clustering by
	// shard_name would let a snapshot's zonemaps prune each shard's scan, but it
	// costs a payload-carrying sort of the whole corpus on the blocking startup
	// path — a losing trade at the shard counts seen in practice (single digits to
	// low tens), and only worth revisiting in the hundreds.
	//
	// The join is on native types: ValidateReadToShardSchema enforces that both
	// read_id columns share a type, so VARCHAR/BIGINT/UUID all compare directly.
	return "SELECT rts.shard_name, " + BuildSequenceColumnList(schema, "q.") + " FROM " +
	       KeywordHelper::WriteOptionallyQuoted(query_table) + " q JOIN " +
	       KeywordHelper::WriteOptionallyQuoted(read_to_shard_table) + " rts ON q.read_id = rts.read_id";
}

std::string MaterializeShardedQueryReads(Connection &conn, const std::string &query_table,
                                         const std::string &read_to_shard_table, const SequenceTableSchema &schema) {
	// Uniquified per call: this lands in the *caller's* TEMP catalog (the
	// connection inherits it, which is what lets worker connections see it), so a
	// fixed name would collide across concurrent queries in one session. Name
	// shape follows MaterializeRypeInputTempTable.
	const std::string tmp_name =
	    "_miint_shard_reads_" + StringUtil::Replace(UUID::ToString(UUID::GenerateRandomUUID()), "-", "");
	const std::string tmp_quoted = KeywordHelper::WriteOptionallyQuoted(tmp_name);

	auto create_result = conn.Query("CREATE TEMP TABLE " + tmp_quoted + " AS " +
	                                BuildShardedQueryReadsSelect(query_table, read_to_shard_table, schema));
	if (create_result->HasError()) {
		throw InvalidInputException("Failed to materialize shard-assigned reads from query table '%s': %s", query_table,
		                            create_result->GetError());
	}
	return tmp_name;
}

void ReadShardReadsFrom(ClientContext &context, const std::string &source_sql, const SequenceTableSchema &schema,
                        const std::string &shard_name, miint::SequenceRecordBatch &output) {
	auto conn = MakeReadOnlyHelperConnection(context);

	// shard_name is a user-supplied row value from read_to_shard, so it must be
	// quoted rather than concatenated — same reasoning as ReadShardNameCounts.
	const std::string query = "SELECT " + BuildSequenceColumnList(schema) + " FROM " + source_sql +
	                          " src WHERE src.shard_name = " + KeywordHelper::WriteQuoted(shard_name, '\'');

	auto query_result = conn.Query(query);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read reads for shard '%s': %s", shard_name, query_result->GetError());
	}

	output.clear();
	output.is_paired = schema.has_sequence2;

	auto &materialized = query_result->Cast<MaterializedQueryResult>();
	ProcessQueryResultChunks(materialized, schema, output);
}

QuerySequenceStream::QuerySequenceStream(ClientContext &context, const std::string &table_name,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size)
    : owned_conn_(make_uniq<Connection>(DatabaseInstance::GetDatabase(context))), conn_ptr_(owned_conn_.get()),
      schema_(schema), sub_batch_size_(sub_batch_size), partial_(schema.has_sequence2) {
	// Before InitStream: the stream's own SELECT must be able to resolve a TEMP
	// relation. This constructor owns its connection and creates nothing on it, so
	// inheriting is safe — the Connection& overload below deliberately does not,
	// because those callers pass a connection they created TEMP objects on.
	InheritTempObjects(context, *owned_conn_);
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
	auto conn = MakeReadOnlyHelperConnection(context);
	return LoadSingleEndSequences(conn, table_name, function_name, strict, /*where_sql=*/"");
}

} // namespace duckdb
