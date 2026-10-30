#include "rype_input_stream.hpp"
#include "rype_common.hpp"

#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/common/optional_idx.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/parser/keyword_helper.hpp"

namespace duckdb {

namespace {

// Column order of the Arrow batch, and of `transformed`, which feeds it. Where
// the columns sit in the source scan is `source` (RypeSourceColumns).
constexpr idx_t COL_ID = 0;
constexpr idx_t COL_SEQUENCE = 1;
constexpr idx_t COL_PAIR_SEQUENCE = 2;

//! Stream the relation's `select_cols`. SubmitStream, not Query: the result must
//! stream so memory stays O(batch) instead of materializing the whole corpus.
unique_ptr<StreamingQuery> SubmitScan(Connection &conn, const std::string &select_cols,
                                      const RypeInputStreamOptions &options) {
	auto query_result = SubmitStream(conn, "SELECT " + select_cols + " FROM " + options.relation_quoted);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read sequences from '%s': %s", options.source_name,
		                            query_result->GetError());
	}
	return query_result;
}

} // namespace

idx_t GetNonNegativeOption(ClientContext &context, const char *name, idx_t if_null, const char *zero_means) {
	Value v;
	// Registered at extension load (LoadInternal in src/miint_extension.cpp).
	// Missing means the option was never registered, which is a load-order
	// programmer error rather than a user mistake.
	if (!context.TryGetCurrentSetting(name, v)) {
		throw InternalException("miint extension option '%s' is not registered — was LoadInternal completed?", name);
	}
	if (v.IsNull()) {
		return if_null;
	}
	const auto value = v.GetValue<int64_t>();
	if (value < 0) {
		throw InvalidInputException("%s must be >= 0 (%s)", name, zero_means);
	}
	return NumericCast<idx_t>(value);
}

idx_t GetRypeArrowBatchBytes(ClientContext &context) {
	return GetNonNegativeOption(context, "miint_rype_arrow_batch_bytes", RYPE_ARROW_BATCH_BYTES,
	                            "0 disables the Arrow batch byte ceiling");
}

// ============================================================================
// Construction
// ============================================================================

RypeInputStream::RypeInputStream(unique_ptr<StreamingQuery> result_p, optional_ptr<RypeIdMap> id_map_p,
                                 RypeInputStreamOptions options_p, RypeSourceColumns source_p)
    : result(std::move(result_p)), id_map(id_map_p), options(std::move(options_p)), source(source_p) {
	D_ASSERT(result);
	appender_capacity =
	    options.batch_bytes > 0 ? MinValue<idx_t>(options.batch_rows, STANDARD_VECTOR_SIZE) : options.batch_rows;
	client_properties = result->GetClientProperties();

	// The Arrow schema RYpe validates against (ext/rype/src/arrow/schema.rs):
	// id Int64, sequence + optional pair_sequence as a binary/string type. BLOB
	// exports as LargeBinary because ConfigureRypeArrowExport pinned the
	// connection to 64-bit offsets (#222).
	arrow_types = {LogicalType::BIGINT, LogicalType::BLOB};
	arrow_names = {"id", "sequence"};
	if (options.include_pair_column) {
		arrow_types.push_back(LogicalType::BLOB);
		arrow_names.push_back("pair_sequence");
	}

	D_ASSERT(client_properties.client_context);
	extension_types = ArrowTypeExtensionData::GetExtensionTypes(*client_properties.client_context, arrow_types);
	transformed.Initialize(*client_properties.client_context, arrow_types);

	ArmStream();
}

void RypeInputStream::ArmStream() {
	stream.private_data = this;
	stream.get_schema = StreamGetSchema;
	stream.get_next = StreamGetNext;
	stream.release = StreamRelease;
	stream.get_last_error = StreamGetLastError;
}

unique_ptr<RypeInputStream> BuildRypeInputStream(Connection &conn, RypeIdMap &id_map, RypeInputStreamOptions options) {
	// One scan, carrying the identifier and the sequence together. No ORDER BY:
	// RYpe echoes each row's surrogate id back as query_id and never relies on
	// input row order, and the id <-> sequence pairing travels within the row.
	//
	// sequence2 is projected as a typed NULL when the caller's relation has no
	// such column, matching what the old TEMP table did, so the batch shape RYpe
	// sees does not depend on the source schema.
	//
	// The sequences are read as they are, not cast to BLOB in SQL: AppendSlice
	// reinterprets a VARCHAR in place, which copies nothing, where a cast would
	// also read a backslash escape in it ('\x43') as the byte it names.
	std::string select_cols = options.id_column_quoted + " AS read_id, sequence1 AS sequence";
	if (options.include_pair_column) {
		select_cols += options.has_sequence2 ? ", sequence2 AS pair_sequence" : ", NULL::BLOB AS pair_sequence";
	}

	auto query_result = SubmitScan(conn, select_cols, options);
	// Bind refuses these types already; this check is what AppendSlice's
	// reinterpret rests on, for a relation redefined since bind.
	const auto &types = query_result->GetTypes();
	const RypeSourceColumns source;
	RequireRypeSequenceType(types[source.sequence], "sequence1", options.source_name);
	if (options.include_pair_column) {
		RequireRypeSequenceType(types[source.pair_sequence], "sequence2", options.source_name);
	}
	return make_uniq<RypeInputStream>(std::move(query_result), id_map, std::move(options));
}

unique_ptr<RypeInputStream> BuildRypeFilterStream(Connection &conn, const vector<string> &columns,
                                                  RypeInputStreamOptions options) {
	// The caller's columns as they are, so a group's rows go back out in their own
	// types; the sequences are reinterpreted, as BuildRypeInputStream's are.
	optional_idx sequence1;
	optional_idx sequence2;
	std::string select_cols;
	for (idx_t i = 0; i < columns.size(); i++) {
		if (StringUtil::CIEquals(columns[i], "sequence1")) {
			sequence1 = i;
		} else if (StringUtil::CIEquals(columns[i], "sequence2")) {
			sequence2 = i;
		}
		select_cols += (i == 0 ? "" : ", ") + SQLIdentifier::ToString(columns[i]);
	}
	if (!sequence1.IsValid()) {
		throw InvalidInputException("'%s' has no sequence1 column", options.source_name);
	}
	RypeSourceColumns source;
	source.sequence = sequence1.GetIndex();
	options.has_sequence2 = sequence2.IsValid();
	options.include_pair_column = options.has_sequence2;
	if (options.has_sequence2) {
		source.pair_sequence = sequence2.GetIndex();
	}

	auto query_result = SubmitScan(conn, select_cols, options);
	// As in BuildRypeInputStream: bind refuses these, and AppendSlice rests on it.
	const auto &types = query_result->GetTypes();
	RequireRypeSequenceType(types[source.sequence], "sequence1", options.source_name);
	if (options.has_sequence2) {
		RequireRypeSequenceType(types[source.pair_sequence], "sequence2", options.source_name);
	}
	return make_uniq<RypeInputStream>(std::move(query_result), nullptr, std::move(options), source);
}

// ============================================================================
// C-ABI callbacks
// ============================================================================

int RypeInputStream::StreamGetSchema(ArrowArrayStream *stream, ArrowSchema *out) {
	if (!stream->release) {
		return -1;
	}
	auto &self = *reinterpret_cast<RypeInputStream *>(stream->private_data);
	out->release = nullptr;
	try {
		ArrowConverter::ToArrowSchema(out, self.arrow_types, self.arrow_names, self.client_properties);
	} catch (std::exception &e) {
		self.last_error = ErrorData(e);
		return -1;
	}
	return 0;
}

int RypeInputStream::StreamGetNext(ArrowArrayStream *stream, ArrowArray *out) {
	if (!stream->release) {
		return -1;
	}
	auto &self = *reinterpret_cast<RypeInputStream *>(stream->private_data);
	try {
		self.FetchBatch(out);
	} catch (std::exception &e) {
		self.last_error = ErrorData(e);
		return -1;
	}
	return 0;
}

void RypeInputStream::StreamRelease(ArrowArrayStream *stream) {
	if (!stream || !stream->release) {
		return;
	}
	stream->release = nullptr;
	// Deliberately does not delete the RypeInputStream — see the ownership note
	// on RypeInputStream::stream. Free what the scan holds, which is the part
	// worth releasing promptly; the shell is destroyed with the GlobalState.
	//
	// RYpe calls this on a by-value copy of the struct (arrow-rs's from_raw
	// moves the contents out and leaves an empty stream behind), so the copy's
	// private_data is the only route back to the object. Nulling release on the
	// copy is what stops a second call.
	auto *self = reinterpret_cast<RypeInputStream *>(stream->private_data);
	// A filter stream's release ends the group, not the scan: the next group
	// carries on from where this one stopped, and the rows it sent stay held for
	// TakeGroup.
	if (self && self->id_map) {
		self->ReleaseScan();
	}
}

void RypeInputStream::ReleaseScan() {
	// transformed's sequence columns share the last source chunk's vectors, so
	// they keep that chunk's buffers (and its string heap) alive. Dropping the scan
	// without this would pin one DataChunk of sequence data for the rest of the
	// query — negligible for short reads, hundreds of MB for long ones. Reset
	// restores the vectors from their own caches, releasing the references; the
	// !result guard in FetchBatch means nothing reads transformed afterwards.
	transformed.Reset();
	current_chunk.reset();
	held.clear();
	result.reset();
}

const char *RypeInputStream::StreamGetLastError(ArrowArrayStream *stream) {
	if (!stream->release) {
		return "stream was released";
	}
	D_ASSERT(stream->private_data);
	return reinterpret_cast<RypeInputStream *>(stream->private_data)->last_error.Message().c_str();
}

// ============================================================================
// Batch construction
// ============================================================================

bool RypeInputStream::EnsureChunk() {
	if (current_chunk && chunk_offset < current_chunk->size()) {
		return true;
	}
	current_chunk = NextChunk();
	chunk_offset = 0;
	return current_chunk != nullptr;
}

shared_ptr<DataChunk> RypeInputStream::NextChunk() {
	std::unique_lock<std::mutex> guard(ahead_lock);
	while (ahead.empty() && !scan_ended) {
		if (aborted) {
			throw InterruptException();
		}
		if (fetching) {
			// A ReadAhead on another thread is fetching the next chunk.
			ahead_progress.wait(guard);
			continue;
		}
		fetching = true;
		FetchLocked(guard);
	}
	if (ahead.empty()) {
		// Rethrown with its original type.
		if (scan_error.HasError()) {
			scan_error.Throw();
		}
		return nullptr;
	}
	auto chunk = std::move(ahead.front().chunk);
	ahead_bytes -= ahead.front().bytes;
	ahead.pop_front();
	// Half of what may be read ahead is left: fetch more before the pulls run out.
	const bool start = start_read_ahead && !reading_ahead && !scan_ended && ahead_bytes <= read_ahead_bytes / 2;
	if (start) {
		reading_ahead = true;
	}
	guard.unlock();
	if (start) {
		start_read_ahead();
	}
	return chunk;
}

void RypeInputStream::FetchLocked(std::unique_lock<std::mutex> &guard) {
	D_ASSERT(fetching);
	guard.unlock();
	shared_ptr<DataChunk> chunk;
	ErrorData error;
	try {
		chunk = result->Fetch();
		if (chunk && chunk->size() == 0) {
			chunk.reset();
		}
		// A stream reports an execution error by ending early with HasError() set.
		if (!chunk && result->HasError()) {
			error = result->GetErrorObject();
		}
	} catch (std::exception &e) {
		error = ErrorData(e);
	}
	const idx_t bytes = chunk ? chunk->GetAllocationSize() : 0;
	guard.lock();
	fetching = false;
	// On a ReadAhead's thread nothing would see this throw, and the pulls would
	// carry on past the chunk it lost: end the scan with the error instead.
	try {
		if (chunk) {
			ahead.push_back(AheadChunk {chunk, bytes});
		}
	} catch (std::exception &e) {
		chunk.reset();
		error = ErrorData(e);
	}
	if (chunk) {
		ahead_bytes += bytes;
	} else {
		scan_ended = true;
		scan_error = std::move(error);
	}
	ahead_progress.notify_all();
}

void RypeInputStream::SetBatchRows(idx_t rows) {
	D_ASSERT(rows > 0);
	options.batch_rows = rows;
	// Mirrors the constructor: with a byte ceiling a batch closes well short of
	// the row ceiling, so reserving batch_rows worth of offsets would waste far
	// more than any batch uses.
	appender_capacity = options.batch_bytes > 0 ? MinValue<idx_t>(rows, STANDARD_VECTOR_SIZE) : rows;
}

void RypeInputStream::SetGroupRows(idx_t rows) {
	D_ASSERT(!id_map);
	D_ASSERT(rows > 0);
	group_rows = rows;
}

bool RypeInputStream::BeginGroup() {
	// Without a group size every group would be empty, and a caller looping until
	// BeginGroup returns false would never stop.
	if (id_map || group_rows == 0) {
		throw InternalException("RypeInputStream::BeginGroup needs a filter stream whose group size is set");
	}
	D_ASSERT(held.empty());
	if (!result || !EnsureChunk()) {
		ReleaseScan();
		return false;
	}
	group_remaining = group_rows;
	ArmStream();
	return true;
}

void RypeInputStream::SetReadAhead(idx_t bytes, std::function<void()> start) {
	// A stream with an id map drops its scan when RYpe releases it, which a
	// ReadAhead on another thread could be fetching from.
	if (id_map) {
		throw InternalException("RypeInputStream::SetReadAhead needs a filter stream");
	}
	std::lock_guard<std::mutex> guard(ahead_lock);
	read_ahead_bytes = bytes;
	start_read_ahead = std::move(start);
}

void RypeInputStream::ReadAhead() {
	std::unique_lock<std::mutex> guard(ahead_lock);
	while (!aborted && !scan_ended && ahead_bytes < read_ahead_bytes) {
		if (fetching) {
			// A pull is fetching for itself; carry on after it.
			ahead_progress.wait(guard);
			continue;
		}
		fetching = true;
		FetchLocked(guard);
	}
	reading_ahead = false;
}

idx_t RypeInputStream::ReadAheadBytes() {
	std::lock_guard<std::mutex> guard(ahead_lock);
	return ahead_bytes;
}

void RypeInputStream::Abort() {
	{
		std::lock_guard<std::mutex> guard(ahead_lock);
		aborted = true;
	}
	ahead_progress.notify_all();
}

vector<RypeHeldRows> RypeInputStream::TakeGroup() {
	D_ASSERT(!id_map);
	vector<RypeHeldRows> taken;
	std::swap(taken, held);
	return taken;
}

RypeInputStream::Sizing RypeInputStream::SampleSizing(size_t fallback_read_length) {
	Sizing sizing {fallback_read_length, false, false};
	if (!result || !EnsureChunk()) {
		return sizing;
	}

	// EnsureChunk leaves the chunk current with its offset untouched, so the rows
	// measured here are the same rows the first get_next will append. Nothing is
	// consumed and the relation is read once.
	auto &chunk = *current_chunk;
	const idx_t chunk_size = chunk.size();
	const idx_t from = chunk_offset;

	UnifiedVectorFormat sequence_format;
	chunk.data[source.sequence].ToUnifiedFormat(sequence_format);
	auto sequence_data = UnifiedVectorFormat::GetData<string_t>(sequence_format);

	idx_t total_bytes = 0;
	idx_t counted = 0;
	for (idx_t i = from; i < chunk_size; i++) {
		const idx_t idx = sequence_format.sel->get_index(i);
		if (sequence_format.validity.RowIsValid(idx)) {
			total_bytes += sequence_data[idx].GetSize();
			counted++;
		}
	}
	if (counted > 0) {
		sizing.sampled = true;
		// rype_calculate_batch_config requires avg_read_length > 0; a chunk of
		// empty sequences would otherwise round to zero.
		sizing.avg_read_length = MaxValue<size_t>(total_bytes / counted, 1);
		sizing.held_bytes_per_row = chunk.GetAllocationSize() / chunk_size;
	}

	// A relation with no sequence2 column is projected as NULL::BLOB by
	// BuildRypeInputStream, so its validity is uniformly false and the loop below
	// would reach the same answer. Skipping it makes "known, not estimated"
	// explicit.
	if (options.include_pair_column && options.has_sequence2) {
		UnifiedVectorFormat pair_format;
		chunk.data[source.pair_sequence].ToUnifiedFormat(pair_format);
		for (idx_t i = from; i < chunk_size; i++) {
			if (pair_format.validity.RowIsValid(pair_format.sel->get_index(i))) {
				sizing.is_paired = true;
				break;
			}
		}
	}
	return sizing;
}

void RypeInputStream::FetchBatch(ArrowArray *out) {
	out->release = nullptr;
	if (aborted) {
		throw InterruptException();
	}

	// The released-stream guard in StreamGetNext already covers a well-behaved
	// consumer; this keeps a pull after ReleaseScan a clean end-of-stream rather
	// than a null dereference.
	if (!result) {
		return;
	}

	// A filter stream's group ends after group_rows rows, wherever that falls in a
	// chunk, and the end-of-stream comes before the scan is touched again.
	const idx_t batch_rows = id_map ? options.batch_rows : MinValue<idx_t>(options.batch_rows, group_remaining);
	if (batch_rows == 0) {
		return;
	}

	ArrowAppender appender(arrow_types, appender_capacity, client_properties, extension_types);
	idx_t appended = 0;
	idx_t batch_bytes = 0;

	// Control flow mirrors DuckDB 1.5's ArrowUtil::TryFetchChunk: resume the partially
	// consumed chunk first, then pull further chunks until the batch is full or
	// the scan is drained. "Full" is a row count and, when configured, a byte
	// ceiling — a source chunk can be left partly consumed by either.
	while (appended < batch_rows) {
		if (!EnsureChunk()) {
			break;
		}
		auto &chunk = *current_chunk;
		const idx_t from = chunk_offset;
		const idx_t row_limit = MinValue<idx_t>(chunk.size(), from + (batch_rows - appended));
		const idx_t first_id = next_id;
		const idx_t to = AppendSlice(appender, chunk, from, row_limit, appended == 0, batch_bytes);
		if (!id_map && to > from) {
			// One run per chunk per group: a batch that carries on where the last
			// one stopped in the same chunk extends its run.
			if (!held.empty() && held.back().chunk == current_chunk && held.back().to == from) {
				held.back().to = to;
			} else {
				held.push_back(RypeHeldRows {current_chunk, from, to, first_id});
			}
		}
		chunk_offset += to - from;
		appended += to - from;
		if (to < row_limit) {
			// The byte ceiling closed the batch mid-chunk; the rest of this chunk
			// is picked up by the next get_next call.
			break;
		}
	}

	if (!id_map) {
		group_remaining -= appended;
	}
	if (appended > 0) {
		// Batches under a byte ceiling all hold a similar number of rows, so the
		// one just built is the best available hint for the next.
		appender_capacity = MaxValue<idx_t>(appender_capacity, appended);
		*out = appender.Finalize();
	}
}

idx_t RypeInputStream::AppendSlice(ArrowAppender &appender, DataChunk &chunk, idx_t from, idx_t row_limit,
                                   bool require_progress, idx_t &batch_bytes) {
	const idx_t chunk_size = chunk.size();

	UnifiedVectorFormat sequence_format;
	chunk.data[source.sequence].ToUnifiedFormat(sequence_format);
	UnifiedVectorFormat pair_format;
	if (options.include_pair_column) {
		chunk.data[source.pair_sequence].ToUnifiedFormat(pair_format);
	}

	// Decide how much of [from, row_limit) fits under the byte ceiling before
	// touching the id map, so the map and the Arrow batch stay the same length.
	idx_t to = row_limit;
	if (options.batch_bytes > 0) {
		auto sequence_data = UnifiedVectorFormat::GetData<string_t>(sequence_format);
		auto pair_data = options.include_pair_column ? UnifiedVectorFormat::GetData<string_t>(pair_format) : nullptr;
		to = from;
		for (idx_t i = from; i < row_limit; i++) {
			idx_t row_bytes = 0;
			const idx_t sequence_idx = sequence_format.sel->get_index(i);
			if (sequence_format.validity.RowIsValid(sequence_idx)) {
				row_bytes += sequence_data[sequence_idx].GetSize();
			}
			if (pair_data) {
				const idx_t pair_idx = pair_format.sel->get_index(i);
				if (pair_format.validity.RowIsValid(pair_idx)) {
					row_bytes += pair_data[pair_idx].GetSize();
				}
			}
			const bool forced = require_progress && i == from;
			if (!forced && batch_bytes + row_bytes > options.batch_bytes) {
				break;
			}
			batch_bytes += row_bytes;
			to = i + 1;
		}
		if (to == from) {
			return from;
		}
	}

	// The surrogate for logical row i is base + (i - from): AppendRange appends in
	// logical row order, so this stays in lockstep with the id map.
	const idx_t base = next_id;
	if (id_map) {
		D_ASSERT(base == id_map->size());
		UnifiedVectorFormat id_format;
		chunk.data[source.id].ToUnifiedFormat(id_format);
		id_map->AppendRange(id_format, from, to);
	}
	next_id += to - from;

	// Reject a NULL sequence at the boundary miint owns. RYpe rejects it too, but
	// only by position within an internal record batch the caller cannot see
	// (#243 reported "Unexpected null in column 'sequence' at row 36504" against a
	// relation the reporter had already verified was null-free). Naming the
	// identifier and the absolute row makes the report actionable.
	for (idx_t i = from; i < to; i++) {
		if (!sequence_format.validity.RowIsValid(sequence_format.sel->get_index(i))) {
			const idx_t surrogate = base + (i - from);
			// A filter stream keeps no id map, so it has only the row to name.
			const std::string id =
			    id_map ? StringUtil::Format(" (%s = %s)", options.id_column_name, id_map->ToString(surrogate)) : "";
			throw InvalidInputException(
			    "sequence1 is NULL in '%s' at row %llu%s. RYpe cannot classify a NULL sequence — filter "
			    "these rows out (WHERE sequence1 IS NOT NULL) or supply a value.",
			    options.source_name, static_cast<unsigned long long>(surrogate), id);
		}
	}

	// Column 0 is synthesized; the sequence columns reference the source vectors,
	// so nothing is copied here — ArrowAppender copies straight from the scan's
	// buffers into the Arrow batch. transformed is deliberately not Reset(): every
	// index in [from, to) of the id vector is written below, and the sequence
	// columns are re-referenced each call. Reinterpret, not Reference: a sequence
	// may be VARCHAR, whose physical type BLOB shares, and both builders refuse
	// any other type (RequireRypeSequenceType).
	auto ids = FlatVector::GetDataMutable<int64_t>(transformed.data[COL_ID]);
	for (idx_t i = from; i < to; i++) {
		ids[i] = static_cast<int64_t>(base + (i - from));
	}
	transformed.data[COL_SEQUENCE].Reinterpret(chunk.data[source.sequence]);
	if (options.include_pair_column) {
		transformed.data[COL_PAIR_SEQUENCE].Reinterpret(chunk.data[source.pair_sequence]);
	}
	transformed.SetChildCardinality(chunk_size);

	appender.Append(transformed, from, to, chunk_size);
	return to;
}

} // namespace duckdb
