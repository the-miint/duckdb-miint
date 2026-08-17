#pragma once
//
// The Arrow input stream the RYpe table functions feed to
// rype_classify_arrow / rype_log_ratio_arrow / rype_extract_*.
//
// Replaces the previous two-step pipeline (materialize the caller's relation into
// a per-call TEMP table, then export that table through
// ResultArrowArrayStreamWrapper). That copied every sequence byte a second time
// before RYpe saw any of them: measured peak RSS on a 2.50 GB corpus was 11.86 GB,
// of which the TEMP table was ~1.8 GB, and it scaled linearly with the corpus
// (the-miint/Qiita#459).
//
// This reads the caller's relation exactly once. Per source DataChunk it
//   - records the row's identifier in a RypeIdMap, and
//   - appends (surrogate id, sequence[, pair_sequence]) to an Arrow batch,
//     referencing the source sequence vectors rather than copying them into an
//     intermediate.
// The Arrow batch is the single copy of the sequence bytes; RYpe borrows slices
// out of it (ext/rype/src/arrow/input.rs documents the zero-copy contract).
//
// Correctness note. The TEMP table existed to anchor the read_id <-> sequence
// correspondence: two independent scans of the caller's relation may return rows
// in different orders, so pairing a `SELECT read_id` scan with a separate
// sequence scan by row index is unsound (commit 8d62915,
// test/sql/rype_row_id_correspondence.test). A single scan that carries the
// identifier and the sequence together makes the correspondence intrinsic to the
// row rather than an emergent property of two scans, which is a stronger
// guarantee than the TEMP table gave, not a weaker one.
//
// Threading. get_next runs on whichever thread pulls RYpe's output stream. That
// is the table function's Execute, and every RYpe table function pins
// MaxThreads() to 1, so the RypeIdMap is written and read from one thread. RYpe's
// internal rayon parallelism happens inside classification, after the pull
// returns.

#include "rype_id_map.hpp"

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_appender.hpp"
#include "duckdb/common/error_data.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/function/table/arrow/arrow_duck_schema.hpp"
#include "duckdb/main/chunk_scan_state.hpp"
#include "duckdb/main/client_properties.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"

#include <string>

namespace duckdb {

//! Byte ceiling for one Arrow record batch's sequence payload.
//!
//! Why a power of two: ArrowBuffer::reserve rounds every growth to
//! NextPowerOfTwo and grows by realloc
//! (duckdb/src/include/duckdb/common/arrow/arrow_buffer.hpp). A batch that
//! stops just short of 2^k therefore has capacity exactly 2^k, so the
//! over-allocation is at most one read length rather than up to a full
//! doubling. Sizing by rows instead — which is all rype_recommend_batch_size
//! can give us, from a 1000-row sampled mean — leaves the byte size of a batch
//! unbounded, and that byte size is what actually drives peak memory
//! (the-miint/Qiita#459).
//!
//! 256 MiB is large enough that per-batch fixed costs stay negligible and
//! small enough that the ceiling is a constant rather than a fraction of the
//! corpus.
//!
//! MIINT_RYPE_ARROW_BATCH_BYTES overrides this downwards (0 removes the ceiling
//! entirely) — see ResolveBatchBytes in rype_input_stream.cpp.
constexpr idx_t RYPE_ARROW_BATCH_BYTES = 256ULL * 1024 * 1024;

//! Everything the stream needs that is not derivable from the connection.
struct RypeInputStreamOptions {
	//! The caller's relation, already SQL-quoted.
	std::string relation_quoted;
	//! The identifier column, already SQL-quoted.
	std::string id_column_quoted;
	//! Unquoted names used only in error messages.
	std::string source_name;
	std::string id_column_name;
	//! Storage type of the identifier column; drives RypeIdMap.
	LogicalType id_type;
	//! Emit a `pair_sequence` column (classify / log_ratio) or not (extract).
	bool include_pair_column = false;
	//! True when the caller's relation actually has a `sequence2` column. When
	//! false and include_pair_column is true, the column is projected as NULL.
	bool has_sequence2 = false;
	//! Row ceiling for one Arrow record batch.
	idx_t batch_rows = STANDARD_VECTOR_SIZE;
	//! Byte ceiling for one batch's sequence payload; 0 disables it. Whichever
	//! ceiling is reached first closes the batch. A batch always carries at
	//! least one row, so a single read larger than the ceiling still flows.
	idx_t batch_bytes = 0;
};

class RypeInputStream {
public:
	//! The C-ABI handle handed to RYpe. private_data points back at this object.
	//!
	//! This object is NOT owned by the stream. The release callback tears down
	//! the scan and marks the stream released, but never deletes; the owner is
	//! whoever holds the unique_ptr, in practice the table function's
	//! GlobalState.
	//!
	//! Why not transfer ownership to RYpe, which is the usual C Data Interface
	//! arrangement: the rype_*_arrow_ex entry points take ownership of the input
	//! stream *once they begin consuming it*, and a -1 return does not say
	//! whether that happened — argument validation fails before consuming, but a
	//! failure inside ArrowArrayStreamReader::from_raw fails after (see the
	//! ownership contract in ext/rype/rype.h). A caller that deletes on -1 double
	//! -frees in the second case, and one that does not leaks in the first. The
	//! struct cannot be probed afterwards either, since it lives inside the
	//! object a deleting release callback would have freed. Keeping ownership on
	//! our side removes the question: exactly one destruction happens, at
	//! GlobalState teardown, whatever RYpe did.
	ArrowArrayStream stream;

	RypeInputStream(unique_ptr<QueryResult> result, RypeIdMap &id_map, RypeInputStreamOptions options);

	RypeInputStream(const RypeInputStream &) = delete;
	RypeInputStream &operator=(const RypeInputStream &) = delete;

	//! Drop the scan (streaming QueryResult and its chunk state) without
	//! destroying the object. Called from the release callback so the scan's
	//! resources go as soon as RYpe is done, rather than at GlobalState
	//! teardown. Idempotent.
	void ReleaseScan();

private:
	static int StreamGetSchema(ArrowArrayStream *stream, ArrowSchema *out);
	static int StreamGetNext(ArrowArrayStream *stream, ArrowArray *out);
	static void StreamRelease(ArrowArrayStream *stream);
	static const char *StreamGetLastError(ArrowArrayStream *stream);

	//! Build one record batch; leaves out->release null when the input is drained.
	void FetchBatch(ArrowArray *out);
	//! Append rows of `chunk` starting at `from` to `appender`, recording
	//! identifiers, stopping at `row_limit` or when `batch_bytes` would exceed
	//! the byte ceiling — whichever comes first. `batch_bytes` is carried across
	//! calls within one batch and updated in place. `require_progress` forces the
	//! first row to be taken even if it alone exceeds the ceiling, so an
	//! oversized read cannot stall the stream. Returns the exclusive end of the
	//! range actually appended.
	idx_t AppendSlice(ArrowAppender &appender, DataChunk &chunk, idx_t from, idx_t row_limit, bool require_progress,
	                  idx_t &batch_bytes);

	unique_ptr<QueryResult> result;
	unique_ptr<ChunkScanState> scan_state;
	RypeIdMap &id_map;
	RypeInputStreamOptions options;

	ClientProperties client_properties;
	vector<LogicalType> arrow_types;
	vector<string> arrow_names;
	unordered_map<idx_t, const shared_ptr<ArrowTypeExtensionData>> extension_types;

	//! Scratch chunk handed to ArrowAppender: column 0 is the surrogate id we
	//! synthesize, the sequence columns reference the source chunk's vectors.
	DataChunk transformed;

	//! Row capacity hint for the next ArrowAppender. ArrowVarcharData::Initialize
	//! reserves (capacity + 1) * 8 bytes of offsets per sequence column, so under
	//! a byte ceiling — which usually closes a batch well short of batch_rows —
	//! passing batch_rows would reserve far more than any batch uses. Start small
	//! and let each batch's actual row count seed the next.
	idx_t appender_capacity;

	ErrorData last_error;
};

//! Start the single scan of the caller's relation and wrap it in a stream.
//!
//! The returned object must be kept alive by the caller until RYpe is finished
//! with the stream — i.e. until after the RYpe output stream has been released.
//! In practice it is stored in the table function's GlobalState and reset in
//! the destructor body, after `output_stream.release` and before the
//! sub-connection whose QueryResult it holds.
//!
//! `id_map` must outlive the stream, for the same reason and by the same means.
unique_ptr<RypeInputStream> BuildRypeInputStream(Connection &conn, RypeIdMap &id_map, RypeInputStreamOptions options);

} // namespace duckdb
