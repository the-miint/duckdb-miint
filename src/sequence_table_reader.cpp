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
static void ProcessSingleChunk(DataChunk &chunk, const SequenceTableSchema &schema, miint::SequenceRecordBatch &output,
                               std::vector<std::string> &temp_read_ids, std::vector<std::string> &temp_seq1,
                               std::vector<std::string> &temp_seq2) {
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
	}
}

// Helper to build column list for sequence queries based on schema
static std::string BuildSequenceColumnList(const SequenceTableSchema &schema) {
	std::string columns = "read_id, sequence1";
	if (schema.has_sequence2) {
		columns += ", sequence2";
	}
	if (schema.has_qual1) {
		columns += ", qual1";
	}
	if (schema.has_qual2) {
		columns += ", qual2";
	}
	return columns;
}

std::string BuildShardReadsSelect(const std::string &snapshot_table, const std::string &routing_snapshot,
                                  const SequenceTableSchema &schema, idx_t shard_id) {
	// See the header for why membership is a projected IN (MARK join) and not a
	// JOIN. No ORDER BY: alignment does not depend on read order.
	//
	// The IN compares native types: ValidateReadToShardSchema enforces that both
	// read_id columns share a type, so VARCHAR/BIGINT/UUID all compare directly.
	// The shard is selected by its integer id (MaterializeReadToShard), so no
	// user-supplied value is concatenated into this SQL at all.
	const std::string columns = BuildSequenceColumnList(schema);
	return "SELECT " + columns + " FROM (SELECT " + columns + ", read_id IN (SELECT read_id FROM " +
	       KeywordHelper::WriteOptionallyQuoted(routing_snapshot) + " WHERE shard_id = " + std::to_string(shard_id) +
	       ") AS _miint_in_shard FROM " + KeywordHelper::WriteOptionallyQuoted(snapshot_table) +
	       ") WHERE _miint_in_shard";
}

std::string BuildQueryReadsSelect(const std::string &query_table, const SequenceTableSchema &schema) {
	return "SELECT " + BuildSequenceColumnList(schema) + " FROM " + KeywordHelper::WriteOptionallyQuoted(query_table);
}

// Uniquified per call: these TEMP tables land in the *caller's* catalog (the
// connection inherits it, which is what lets worker connections see them), so a
// fixed name would collide across concurrent queries in one session. Name shape
// follows MaterializeRypeInputTempTable.
static std::string UniqueTempRelationName(const std::string &prefix) {
	return prefix + StringUtil::Replace(UUID::ToString(UUID::GenerateRandomUUID()), "-", "");
}

// One streaming pass of `select_sql` into a fresh per-call TEMP table on `conn`.
// Shared by every #229 snapshot taken here, so they cannot drift apart in how
// they bound their own memory, type the destination, or clean up on failure.
static std::string MaterializeSelectIntoTemp(Connection &conn, const std::string &select_sql,
                                             const std::string &name_prefix, const std::string &error_context,
                                             idx_t &out_row_count) {
	const std::string tmp_name = UniqueTempRelationName(name_prefix);
	const std::string tmp_quoted = KeywordHelper::WriteOptionallyQuoted(tmp_name);

	// Stream the source and append each chunk into the snapshot as it arrives,
	// rather than one CREATE TABLE AS SELECT that pulls the whole query relation
	// through the pipeline before this call returns. A single-pass streaming
	// query (SendQuery, not Query) still reads the source exactly once — the
	// #229 guarantee is about pass count, not chunk size — but bounds this
	// materialization's own working set to O(one chunk) + O(the Appender's
	// internal flush buffer) instead of O(relation size). See the multi-part
	// memory bug this fixes: a full-corpus-sized query relation (tens of
	// millions of reads) made this the dominant unmanaged, memory_limit-
	// invisible cost regardless of index-part size or thread count.
	//
	// A dedicated connection drives the read: a Connection supports only one
	// active pending query at a time, and the Appender below issues its own
	// statements against `conn` as it flushes, which would otherwise collide
	// with `conn`'s still-open SendQuery stream mid-loop.
	Connection stream_conn = MakeReadOnlyHelperConnection(*conn.context);
	auto stream = stream_conn.SendQuery(select_sql);
	if (stream->HasError()) {
		throw InvalidInputException("%s: %s", error_context, stream->GetError());
	}

	// The destination's column types come from the stream's own output schema
	// (already resolved by SendQuery's bind, before any row is fetched) rather
	// than a second "CREATE TABLE AS SELECT ... WHERE FALSE" probe query against
	// the source. The source can be a view over something with bind-time work of
	// its own (e.g. read_fastx opening/sniffing the underlying file) — one query
	// against it here means that work happens once, not twice.
	std::string create_sql = "CREATE TEMP TABLE " + tmp_quoted + " (";
	for (idx_t i = 0; i < stream->types.size(); i++) {
		if (i > 0) {
			create_sql += ", ";
		}
		create_sql += KeywordHelper::WriteOptionallyQuoted(stream->names[i]) + " " + stream->types[i].ToString();
		// LogicalType::ToString() renders a collated VARCHAR as plain "VARCHAR",
		// so without this the snapshot silently loses the collation and every
		// later comparison against it turns case-sensitive. That is a silent
		// wrong-answer bug, not an error: with `read_id VARCHAR COLLATE NOCASE`,
		// rows whose ids differ only in case stop matching and are never
		// aligned, and nothing reports it.
		if (stream->types[i].id() == LogicalTypeId::VARCHAR) {
			const auto collation = StringType::GetCollation(stream->types[i]);
			if (!collation.empty()) {
				create_sql += " COLLATE " + KeywordHelper::WriteOptionallyQuoted(collation);
			}
		}
	}
	create_sql += ")";
	auto create_result = conn.Query(create_sql);
	if (create_result->HasError()) {
		throw InvalidInputException("%s: %s", error_context, create_result->GetError());
	}

	// From here on, the empty snapshot table above is committed in conn's TEMP
	// catalog. If the fill below throws partway (a mid-stream query error, or
	// OOM), drop it before propagating — otherwise the caller never receives
	// query_snapshot's name to clean up later and the empty table leaks in the
	// user's session for the lifetime of the connection.
	idx_t row_count = 0;
	try {
		Appender appender(conn, tmp_name);
		while (true) {
			auto chunk = stream->Fetch();
			if (!chunk || chunk->size() == 0) {
				// Fetch() returns null for both a clean end-of-stream AND a
				// mid-stream query error (e.g. a malformed row deep in the
				// source) — HasError() is what tells them apart. Missing
				// this check would silently truncate the snapshot to whatever
				// was read before the failure instead of surfacing it.
				if (stream->HasError()) {
					throw InvalidInputException("%s: %s", error_context, stream->GetError());
				}
				break;
			}
			row_count += chunk->size();
			appender.AppendDataChunk(*chunk);
		}
		appender.Close();
	} catch (...) {
		DropHelperTempRelation(conn, tmp_quoted);
		throw;
	}

	out_row_count = row_count;
	return tmp_name;
}

std::string MaterializeQueryReads(Connection &conn, const std::string &query_table, const SequenceTableSchema &schema,
                                  idx_t &out_row_count) {
	return MaterializeSelectIntoTemp(conn, BuildQueryReadsSelect(query_table, schema), "_miint_query_reads_",
	                                 "Failed to materialize query table '" + query_table + "'", out_row_count);
}

std::string MaterializeReadToShard(Connection &conn, const std::string &read_to_shard_table,
                                   const std::vector<std::string> &shard_names, idx_t &out_row_count) {
	// Routing is encoded to an integer shard id rather than carrying the shard
	// name, following DuckDB's dimension-table pattern
	// (https://duckdb.org/2026/10/02/dimension-tables): work on the narrow key,
	// resolve the string once at the edge. `shard_id` is the index into the
	// caller's shard list, which ActiveShard already carries as shard_idx, so
	// nothing has to be resolved back.
	//
	// Three things this buys, in order of importance:
	//
	// 1. Zonemap pruning that actually prunes. Every shard re-reads this
	//    snapshot once per index part with an equality filter, and pruning is
	//    decided by physical row order (docs/internals/duckdb-engine-notes.md
	//    § Indexes) — but DuckDB keeps only the first 8 bytes of a string in
	//    min/max statistics (StringStatsData::MAX_STRING_MINMAX_SIZE). Shard
	//    names that agree in their first 8 bytes — "mem_shard_0..59", or
	//    accessions like "GCF_000001405.40" — give every row group an identical
	//    min and max, so sorting by the name prunes NOTHING. Measured on 6M
	//    routing rows over 60 shards: filtering the name-sorted snapshot scanned
	//    all 6,000,000 rows, the id-sorted one 182,272. Integer statistics are
	//    exact.
	// 2. The user's collation decides membership. The join below compares
	//    shard names under the routing relation's own collation, before anything
	//    is copied; afterwards only integers are compared. Encoding is therefore
	//    collation-correct even where the snapshot's own DDL is not.
	// 3. A shard name the caller's list does not contain becomes a NULL id
	//    instead of a row that quietly matches no shard. Checked below.
	//
	// The snapshot is also narrower — 4 bytes per routing row instead of a
	// string — which shrinks the sort, its spill, and the temp_directory budget.
	const std::string dim_name = UniqueTempRelationName("_miint_shards_");
	const std::string dim_quoted = KeywordHelper::WriteOptionallyQuoted(dim_name);
	const std::string error_context = "Failed to materialize read_to_shard table '" + read_to_shard_table + "'";

	auto dim_result = conn.Query("CREATE TEMP TABLE " + dim_quoted + " (shard_name VARCHAR, shard_id UINTEGER)");
	if (dim_result->HasError()) {
		throw InvalidInputException("%s: %s", error_context, dim_result->GetError());
	}

	std::string snapshot;
	try {
		Appender dim_appender(conn, dim_name);
		for (idx_t i = 0; i < shard_names.size(); i++) {
			dim_appender.BeginRow();
			dim_appender.Append(Value(shard_names[i]));
			dim_appender.Append(Value::UINTEGER(NumericCast<uint32_t>(i)));
			dim_appender.EndRow();
		}
		dim_appender.Close();

		// LEFT JOIN, not INNER: an unmatched routing row must survive as a NULL
		// so the check below can see it. An INNER JOIN would drop exactly the
		// rows we want to complain about.
		const std::string select_sql = "SELECT r.read_id, s.shard_id FROM " +
		                               KeywordHelper::WriteOptionallyQuoted(read_to_shard_table) + " r LEFT JOIN " +
		                               dim_quoted + " s ON r.shard_name = s.shard_name ORDER BY s.shard_id";
		snapshot = MaterializeSelectIntoTemp(conn, select_sql, "_miint_read_to_shard_", error_context, out_row_count);
	} catch (...) {
		DropHelperTempRelation(conn, dim_quoted);
		throw;
	}
	// The dimension table has done its job once the ids are encoded; only the
	// snapshot is read from here on.
	DropHelperTempRelation(conn, dim_quoted);

	// A NULL id means the routing relation named a shard that shard discovery
	// did not see. Discovery runs at bind time against the user's relation and
	// this snapshot is taken at execution time, so an unstable relation can
	// disagree between them. We cannot make the two agree here, but we can
	// refuse to align a partial result silently.
	//
	// ORDER BY above puts the NULLs last and integer statistics are exact, so
	// this scan prunes to the tail row groups rather than reading the snapshot.
	const std::string snapshot_quoted = KeywordHelper::WriteOptionallyQuoted(snapshot);
	auto orphan_result = conn.Query("SELECT count(*) FROM " + snapshot_quoted + " WHERE shard_id IS NULL");
	if (orphan_result->HasError()) {
		DropHelperTempRelation(conn, snapshot_quoted);
		throw InvalidInputException("%s: %s", error_context, orphan_result->GetError());
	}
	const auto orphans = orphan_result->GetValue(0, 0).GetValue<idx_t>();
	if (orphans > 0) {
		DropHelperTempRelation(conn, snapshot_quoted);
		throw InvalidInputException(
		    "read_to_shard table '%s' changed while align_minimap2_sharded was starting: %llu routing row(s) name a "
		    "shard that was not present when the shard list was read. Materialize '%s' into a table before aligning.",
		    read_to_shard_table, static_cast<unsigned long long>(orphans), read_to_shard_table);
	}
	return snapshot;
}

QuerySequenceStream::QuerySequenceStream(ClientContext &context, const std::string &select_sql,
                                         const std::string &source, const SequenceTableSchema &schema,
                                         idx_t sub_batch_size, std::shared_ptr<void> source_keepalive)
    : owned_conn_(make_uniq<Connection>(DatabaseInstance::GetDatabase(context))), conn_ptr_(owned_conn_.get()),
      source_keepalive_(std::move(source_keepalive)), schema_(schema), sub_batch_size_(sub_batch_size),
      partial_(schema.has_sequence2) {
	// Before InitStream: the stream's own SELECT must be able to resolve a TEMP
	// relation. This constructor owns its connection and creates nothing on it, so
	// inheriting is safe — the Connection& overload below deliberately does not,
	// because those callers pass a connection they created TEMP objects on.
	InheritTempObjects(context, *owned_conn_);
	InitStream(select_sql, source);
}

QuerySequenceStream::QuerySequenceStream(ClientContext &context, const std::string &table_name,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size)
    // Same projection the snapshot was built with, so replaying a snapshot binds
    // against exactly the columns it holds.
    : QuerySequenceStream(context, BuildQueryReadsSelect(table_name, schema), "query table '" + table_name + "'",
                          schema, sub_batch_size) {
}

QuerySequenceStream::QuerySequenceStream(Connection &conn, const std::string &table_name,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size)
    : owned_conn_(nullptr), conn_ptr_(&conn), schema_(schema), sub_batch_size_(sub_batch_size),
      partial_(schema.has_sequence2) {
	InitStream(BuildQueryReadsSelect(table_name, schema_), "query table '" + table_name + "'");
}

QuerySequenceStream::QuerySequenceStream(ClientContext &context, const SelectSql &select,
                                         const SequenceTableSchema &schema, idx_t sub_batch_size,
                                         std::shared_ptr<void> source_keepalive)
    : QuerySequenceStream(context, select.sql, select.label.empty() ? "query '" + select.sql + "'" : select.label,
                          schema, sub_batch_size, std::move(source_keepalive)) {
}

QuerySequenceStream::~QuerySequenceStream() {
	stream_.reset();
}

void QuerySequenceStream::InitStream(const std::string &select_sql, const std::string &source) {
	partial_.reserve(sub_batch_size_);

	stream_ = conn_ptr_->SendQuery(select_sql);
	if (stream_->HasError()) {
		throw InvalidInputException("Failed to read from %s: %s", source, stream_->GetError());
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
			// Fetch() returns null for both a clean end-of-stream AND a mid-stream
			// query error (e.g. out of memory partway through) — HasError() is what
			// tells them apart. Same check as MaterializeQueryReads.
			if (stream_->HasError()) {
				throw InvalidInputException("Failed to read sequences mid-stream: %s", stream_->GetError());
			}
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
			rows_delivered_ += result.size();
			return result;
		}
		// Chunk pushed partial_ past sub_batch_size_ — split: return exactly
		// sub_batch_size_ rows and keep the remainder for the next call.
		miint::SequenceRecordBatch result = partial_.SubRange(0, sub_batch_size_);
		miint::SequenceRecordBatch remainder = partial_.SubRange(sub_batch_size_, partial_.size() - sub_batch_size_);
		partial_ = std::move(remainder);
		rows_delivered_ += result.size();
		return result;
	}

	// Stream exhausted — return whatever remains
	miint::SequenceRecordBatch result = std::move(partial_);
	partial_ = miint::SequenceRecordBatch(schema_.has_sequence2);
	rows_delivered_ += result.size();
	return result;
}

idx_t QuerySequenceStream::RowsDelivered() const {
	std::lock_guard<std::mutex> lock(mutex_);
	return rows_delivered_;
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
