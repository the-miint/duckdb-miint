#include "rype_input_stream.hpp"

#include "duckdb/common/arrow/arrow_converter.hpp"
#include "duckdb/main/chunk_scan_state/query_result.hpp"

#include <cstdlib>

namespace duckdb {

namespace {

// Column order of both the source scan and the Arrow batch. The two agree by
// construction so AppendSlice can reference source vectors positionally.
constexpr idx_t COL_ID = 0;
constexpr idx_t COL_SEQUENCE = 1;
constexpr idx_t COL_PAIR_SEQUENCE = 2;

//! Apply the MIINT_RYPE_ARROW_BATCH_BYTES override to a configured byte ceiling.
//!
//! Two reasons this override exists, both about being able to observe the
//! ceiling's effect rather than take it on faith:
//!
//!   * Lowering it. At the 256 MiB default the multi-batch path — surrogate ids
//!     and id-map entries staying in lockstep across Arrow batch boundaries —
//!     would only ever run on corpora far too large to keep as test fixtures.
//!     test/sql/rype_input_stream_batching.test sets a few hundred bytes.
//!   * Disabling it, with 0. That restores one unbounded Arrow batch, which is
//!     the A/B baseline for the memory claim and an escape hatch if the ceiling
//!     ever turns out to hurt a workload.
//!
//! The override never *raises* a ceiling and never introduces one where the
//! caller did not ask for one.
idx_t ResolveBatchBytes(idx_t configured) {
	if (configured == 0) {
		return 0;
	}
	const char *env = std::getenv("MIINT_RYPE_ARROW_BATCH_BYTES");
	if (!env) {
		return configured;
	}
	try {
		const auto parsed = std::stoull(env);
		if (parsed < configured) {
			return NumericCast<idx_t>(parsed);
		}
	} catch (const std::exception &) {
		// Unparseable: keep the configured ceiling rather than silently
		// disabling the thing that bounds memory.
	}
	return configured;
}

} // namespace

// ============================================================================
// Construction
// ============================================================================

RypeInputStream::RypeInputStream(unique_ptr<QueryResult> result_p, RypeIdMap &id_map_p,
                                 RypeInputStreamOptions options_p)
    : result(std::move(result_p)), id_map(id_map_p), options(std::move(options_p)) {
	D_ASSERT(result);
	options.batch_bytes = ResolveBatchBytes(options.batch_bytes);
	appender_capacity =
	    options.batch_bytes > 0 ? MinValue<idx_t>(options.batch_rows, STANDARD_VECTOR_SIZE) : options.batch_rows;
	client_properties = result->client_properties;
	scan_state = make_uniq<QueryResultChunkScanState>(*result);

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
	std::string select_cols = options.id_column_quoted + " AS read_id, sequence1::BLOB AS sequence";
	if (options.include_pair_column) {
		select_cols += options.has_sequence2 ? ", sequence2::BLOB AS pair_sequence" : ", NULL::BLOB AS pair_sequence";
	}

	// SendQuery, not Query: the result must stream so memory stays O(batch)
	// instead of materializing the whole corpus. A prepared statement would
	// force-materialize even with allow_stream_result.
	auto query_result = conn.SendQuery("SELECT " + select_cols + " FROM " + options.relation_quoted);
	if (query_result->HasError()) {
		throw InvalidInputException("Failed to read sequences from '%s': %s", options.source_name,
		                            query_result->GetError());
	}
	return make_uniq<RypeInputStream>(std::move(query_result), id_map, std::move(options));
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
	if (self) {
		self->ReleaseScan();
	}
}

void RypeInputStream::ReleaseScan() {
	// Order matters: the scan state borrows the QueryResult.
	scan_state.reset();
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

void RypeInputStream::FetchBatch(ArrowArray *out) {
	out->release = nullptr;

	// The released-stream guard in StreamGetNext already covers a well-behaved
	// consumer; this keeps a pull after ReleaseScan a clean end-of-stream rather
	// than a null dereference.
	if (!result) {
		return;
	}

	if (result->HasError()) {
		result->GetErrorObject().Throw();
	}

	ArrowAppender appender(arrow_types, appender_capacity, client_properties, extension_types);
	idx_t appended = 0;
	idx_t batch_bytes = 0;

	// Control flow mirrors ArrowUtil::TryFetchChunk: resume the partially
	// consumed chunk first, then pull further chunks until the batch is full or
	// the scan is drained. "Full" is a row count and, when configured, a byte
	// ceiling — a source chunk can be left partly consumed by either.
	while (appended < options.batch_rows) {
		if (scan_state->RemainingInChunk() == 0) {
			ErrorData error;
			if (!scan_state->LoadNextChunk(error)) {
				if (scan_state->HasError()) {
					error = scan_state->GetError();
				}
				error.Throw();
			}
			if (scan_state->ChunkIsEmpty() || scan_state->Finished()) {
				break;
			}
		}
		auto &chunk = scan_state->CurrentChunk();
		const idx_t from = scan_state->CurrentOffset();
		const idx_t row_limit = MinValue<idx_t>(chunk.size(), from + (options.batch_rows - appended));
		const idx_t to = AppendSlice(appender, chunk, from, row_limit, appended == 0, batch_bytes);
		scan_state->IncreaseOffset(to - from);
		appended += to - from;
		if (to < row_limit) {
			// The byte ceiling closed the batch mid-chunk; the rest of this chunk
			// is picked up by the next get_next call.
			break;
		}
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

	UnifiedVectorFormat id_format;
	chunk.data[COL_ID].ToUnifiedFormat(chunk_size, id_format);
	UnifiedVectorFormat sequence_format;
	chunk.data[COL_SEQUENCE].ToUnifiedFormat(chunk_size, sequence_format);
	UnifiedVectorFormat pair_format;
	if (options.include_pair_column) {
		chunk.data[COL_PAIR_SEQUENCE].ToUnifiedFormat(chunk_size, pair_format);
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
	const idx_t base = id_map.size();
	id_map.AppendRange(id_format, from, to);

	// Reject a NULL sequence at the boundary miint owns. RYpe rejects it too, but
	// only by position within an internal record batch the caller cannot see
	// (#243 reported "Unexpected null in column 'sequence' at row 36504" against a
	// relation the reporter had already verified was null-free). Naming the
	// identifier and the absolute row makes the report actionable.
	for (idx_t i = from; i < to; i++) {
		if (!sequence_format.validity.RowIsValid(sequence_format.sel->get_index(i))) {
			const idx_t surrogate = base + (i - from);
			throw InvalidInputException(
			    "sequence1 is NULL in '%s' at row %llu (%s = %s). RYpe cannot classify a NULL sequence — filter "
			    "these rows out (WHERE sequence1 IS NOT NULL) or supply a value.",
			    options.source_name, static_cast<unsigned long long>(surrogate), options.id_column_name,
			    id_map.ToString(surrogate));
		}
	}

	// Column 0 is synthesized; the sequence columns reference the source vectors,
	// so nothing is copied here — ArrowAppender copies straight from the scan's
	// buffers into the Arrow batch. transformed is deliberately not Reset(): every
	// index in [from, to) of the id vector is written below, and the sequence
	// columns are re-referenced each call.
	auto ids = FlatVector::GetData<int64_t>(transformed.data[COL_ID]);
	for (idx_t i = from; i < to; i++) {
		ids[i] = static_cast<int64_t>(base + (i - from));
	}
	transformed.data[COL_SEQUENCE].Reference(chunk.data[COL_SEQUENCE]);
	if (options.include_pair_column) {
		transformed.data[COL_PAIR_SEQUENCE].Reference(chunk.data[COL_PAIR_SEQUENCE]);
	}
	transformed.SetCardinality(chunk_size);

	appender.Append(transformed, from, to, chunk_size);
	return to;
}

} // namespace duckdb
