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
// Threading. get_next runs on whichever thread pulls RYpe's output stream. A
// stream with an id map is pulled by its table function's Execute, which pins
// MaxThreads() to 1, so the RypeIdMap is written and read from one thread. A
// filter stream keeps no id map, and one classify job at a time pulls it. RYpe's
// internal rayon parallelism happens inside classification, after the pull
// returns.
//
// Filter streams (BuildRypeFilterStream) read the caller's own columns instead,
// keep no id map, and end after a set number of rows so that one RYpe call sees
// exactly one group. The rows a group sent are held, shared rather than copied,
// for the caller to take once RYpe has classified them. A filter stream can also
// read ahead (SetReadAhead): another thread fetches the scan's chunks while RYpe
// classifies.

#include "rype_id_map.hpp"

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_appender.hpp"
#include "duckdb/common/error_data.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/function/table/arrow/arrow_duck_schema.hpp"
#include "miint_streaming_query.hpp"
#include "duckdb/main/client_properties.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/query_result.hpp"

#include <atomic>
#include <condition_variable>
#include <deque>
#include <functional>
#include <mutex>
#include <string>

namespace duckdb {

//! Default byte ceiling for one Arrow record batch's sequence payload, and the
//! default of the `miint_rype_arrow_batch_bytes` setting registered in
//! LoadInternal.
//!
//! Why a power of two: ArrowBuffer::reserve rounds every growth to
//! NextPowerOfTwo and grows by realloc
//! (duckdb/src/include/duckdb/common/arrow/arrow_buffer.hpp). A single-column
//! batch that stops just short of 2^k therefore has capacity exactly 2^k, so
//! the over-allocation is at most one read length rather than up to a full
//! doubling. Sizing by rows instead — which is all rype_recommend_batch_size
//! can give us, from a sampled mean — leaves the byte size of a batch unbounded,
//! and that byte size is what actually drives peak memory
//! (the-miint/Qiita#459).
//!
//! The exact-2^k argument holds for one buffer. A paired batch splits the
//! budget across two aux buffers that ArrowBuffer rounds up independently, so
//! the worst case is 1.5x the ceiling: mates of unequal length landing at, say,
//! 160 MiB and 96 MiB round to 256 + 128 MiB. That is bounded and constant,
//! which is the property being bought here; it is not the "one read length"
//! bound that a single column enjoys.
//!
//! 256 MiB is large enough that per-batch fixed costs stay negligible and
//! small enough that the ceiling is a constant rather than a fraction of the
//! corpus.
constexpr idx_t RYPE_ARROW_BATCH_BYTES = 256ULL * 1024 * 1024;

//! Read the session-scoped `miint_rype_arrow_batch_bytes` setting: the byte
//! ceiling above, 0 to disable it entirely.
//!
//! A DuckDB setting rather than an environment variable because this has to be
//! reachable from a SQL test in every lane. The multi-batch path — surrogate ids
//! and id-map entries staying in lockstep across Arrow batch boundaries — is the
//! load-bearing invariant of this stream, and at 256 MiB it would only ever run
//! on corpora far too large to keep as fixtures. Gating that test on an env var
//! kept it out of `make test` and out of every CI lane but one; and getenv on
//! the ingest path meant a typo like `1G` parsed as 1 byte with no diagnostic.
//! DuckDB parses and range-checks the value for us.
idx_t GetRypeArrowBatchBytes(ClientContext &context);

//! A BIGINT extension option, registered in LoadInternal, that must not be
//! negative: `if_null` when it is NULL. `zero_means` says, in the error for a
//! negative value, what 0 does.
idx_t GetNonNegativeOption(ClientContext &context, const char *name, idx_t if_null, const char *zero_means);

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

//! A run of one scanned chunk's rows that a filter stream sent to Arrow, held for
//! the caller. The chunk is shared, not copied, so a chunk that straddles two
//! groups appears in both, once in each.
struct RypeHeldRows {
	//! Every column of the scan, in the relation's own types.
	shared_ptr<DataChunk> chunk;
	//! Rows [from, to) of `chunk`, in scan order.
	idx_t from;
	idx_t to;
	//! The Arrow id of row `from`; row i's id is first_id + (i - from).
	idx_t first_id;
};

//! Where the scan puts the columns a RypeInputStream reads. The defaults are the
//! (read_id, sequence[, pair_sequence]) scan BuildRypeInputStream writes.
struct RypeSourceColumns {
	//! Read only when there is an id map.
	idx_t id = 0;
	idx_t sequence = 1;
	//! Read only when options.include_pair_column is set.
	idx_t pair_sequence = 2;
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

	//! A null `id_map` makes a filter stream: no identifier is read, a row's Arrow
	//! id is its position in the scan, and rows go out in groups (BeginGroup).
	RypeInputStream(unique_ptr<StreamingQuery> result, optional_ptr<RypeIdMap> id_map, RypeInputStreamOptions options,
	                RypeSourceColumns source = RypeSourceColumns());

	RypeInputStream(const RypeInputStream &) = delete;
	RypeInputStream &operator=(const RypeInputStream &) = delete;

	//! What RYpe's batch sizing needs to know about the input, read from the
	//! front of this stream's own scan.
	struct Sizing {
		size_t avg_read_length;
		bool is_paired;
		//! False when the relation turned out to be empty, so the caller can tell
		//! a real measurement from the fallback.
		bool sampled;
		//! What holding one scanned row costs, every column counted: the first
		//! chunk's allocation over its rows. 0 when not sampled. A filter stream
		//! holds a group of these for its caller.
		idx_t held_bytes_per_row = 0;
	};

	//! Measure `avg_read_length` and `is_paired` from the first source chunk
	//! WITHOUT consuming it — the chunk stays queued for the first get_next.
	//!
	//! These used to be two extra queries against the caller's relation
	//! (SampleAvgReadLength's LIMIT 1000, and a full scan for sequence2 content).
	//! That is forbidden: a relation passed by name is not guaranteed to return
	//! the same rows twice, and a registered Arrow relation returns *zero* rows on
	//! a second read, which is indistinguishable from end-of-input (#229, and
	//! docs/internals/reading-tables-views.md "Read the relation ONCE"). The old
	//! TEMP-table design could afford those probes because it had already taken
	//! its one read and was querying the snapshot; removing the snapshot without
	//! moving the probes turned one read of the caller's relation into three.
	//!
	//! Both facts are therefore estimates taken from the head of the single pass.
	//! For avg_read_length that matches what the LIMIT 1000 sample always did. For
	//! is_paired it is weaker than the full scan #199 introduced: a relation whose
	//! first chunk is single-ended but which pairs up later is under-budgeted, the
	//! direction that OOMs. That is the deliberate trade — correctness of the data
	//! outranks accuracy of a sizing hint — and `max_memory` is the explicit lever
	//! for anyone who hits it. When the relation has no sequence2 column at all,
	//! is_paired is known false without sampling and nothing is estimated.
	Sizing SampleSizing(size_t fallback_read_length = 300);

	//! Set the row ceiling for one Arrow record batch.
	//!
	//! Deferred configuration rather than a constructor argument because the
	//! ceiling comes from RYpe's batch sizing, which needs SampleSizing's numbers,
	//! which need this stream to exist. Must be called before RYpe pulls the first
	//! batch; SampleSizing does not consume anything, so that ordering is free.
	void SetBatchRows(idx_t rows);

	//! Filter streams only: the most rows one group sends before the stream reports
	//! end-of-stream. Deferred for the same reason as SetBatchRows; must be set
	//! before the first BeginGroup.
	void SetGroupRows(idx_t rows);

	//! Filter streams only. Arm the stream for the next group and return true, or
	//! return false when the scan has no rows left, beginning nothing. Call it
	//! before every rype_*_arrow_ex call, the first one included: that call moves
	//! the struct out of `stream` (arrow-rs's from_raw leaves an empty one behind)
	//! and releases its copy, and this puts the callbacks back. Releasing ends the
	//! group but keeps the scan, which carries on, mid-chunk if need be, at the
	//! next group. The previous group's rows must have been taken.
	bool BeginGroup();

	//! Filter streams only: the rows the current group has sent to Arrow, in scan
	//! order, moved out to the caller.
	vector<RypeHeldRows> TakeGroup();

	//! The types of the scan's columns, as the query bound them. Only while the
	//! scan is open (before ReleaseScan).
	const vector<LogicalType> &ScanTypes() const {
		return result->GetTypes();
	}

	//! Drop the scan (streaming QueryResult and its chunk state) without
	//! destroying the object. Called from the release callback so the scan's
	//! resources go as soon as RYpe is done, rather than at GlobalState
	//! teardown. Idempotent.
	void ReleaseScan();

	//! Filter streams only: let the scan run up to `bytes` of chunks ahead of the
	//! batches, fetched by ReadAhead. `start` is called, on the thread taking a
	//! chunk, when what was read ahead has fallen to half of `bytes` and no
	//! ReadAhead is under way or called for; it should have ReadAhead run on another
	//! thread. It is not called again until that ReadAhead returns.
	//!
	//! Why: DuckDB v2.0's streaming result produces as it is fetched from (the
	//! fetching thread runs the scan's pipeline), not ahead of a consumer that is
	//! busy. Fetched only from get_next, the scan stops whenever RYpe classifies.
	//! Raising max_streaming_buffer_size does not change that.
	void SetReadAhead(idx_t bytes, std::function<void()> start);

	//! Fetch the scan's chunks ahead of the batches until those not yet taken come
	//! to the read-ahead's bytes, the scan ends, or the stream is aborted. Runs
	//! alongside get_next on another thread. One fetch is made at a time, so the
	//! chunks keep scan order; and a pull that finds nothing read ahead, and no
	//! fetch under way, fetches for itself rather than wait for a ReadAhead that
	//! may not have started. A scan error is kept for the pull that reaches it.
	void ReadAhead();

	//! The bytes of the chunks read ahead and not yet taken.
	idx_t ReadAheadBytes();

	//! Make every later get_next fail, so a RYpe call reading the stream on
	//! another thread stops at its next batch, and stop any ReadAhead. Safe to call
	//! from any thread.
	void Abort();

private:
	static int StreamGetSchema(ArrowArrayStream *stream, ArrowSchema *out);
	static int StreamGetNext(ArrowArrayStream *stream, ArrowArray *out);
	static void StreamRelease(ArrowArrayStream *stream);
	static const char *StreamGetLastError(ArrowArrayStream *stream);

	//! Make a source chunk with unread rows current, loading the next one if the
	//! current is spent. False at end of stream; throws only when the scan
	//! reported an actual error.
	bool EnsureChunk();
	//! The scan's next chunk: the oldest read ahead, else one fetched here unless a
	//! fetch is under way on another thread. Null at end of stream.
	shared_ptr<DataChunk> NextChunk();
	//! With `fetching` taken and `guard` held: fetch one chunk into `ahead`, or end
	//! the scan, releasing `guard` meanwhile.
	void FetchLocked(std::unique_lock<std::mutex> &guard);
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

	//! Point `stream` at this object's callbacks.
	void ArmStream();

	unique_ptr<StreamingQuery> result;
	//! Cursor over the stream (v1.5 used duckdb's ChunkScanState, removed in v2.0).
	//! Shared so that a filter stream's held rows can keep the chunk.
	shared_ptr<DataChunk> current_chunk;
	idx_t chunk_offset = 0;
	optional_ptr<RypeIdMap> id_map;
	RypeInputStreamOptions options;
	RypeSourceColumns source;
	//! The Arrow id the next row gets; equals id_map->size() when there is a map.
	idx_t next_id = 0;

	//! Filter streams: the group size, the rows the current group may still send,
	//! and the runs it has sent.
	idx_t group_rows = 0;
	idx_t group_remaining = 0;
	vector<RypeHeldRows> held;

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
	//! Set by Abort, from any thread.
	std::atomic<bool> aborted {false};

	//! The scan's fetched chunks not yet taken, and the read-ahead's state. `result`
	//! is fetched from by whichever thread holds `fetching`, so by one at a time.
	struct AheadChunk {
		shared_ptr<DataChunk> chunk;
		idx_t bytes;
	};
	std::mutex ahead_lock;
	//! Signalled when a fetch finishes and when the stream is aborted.
	std::condition_variable ahead_progress;
	// Guarded by ahead_lock.
	std::deque<AheadChunk> ahead;
	idx_t ahead_bytes = 0;
	bool fetching = false;
	//! The scan has no chunk left to fetch; scan_error says whether it failed.
	bool scan_ended = false;
	ErrorData scan_error;
	//! A ReadAhead is called for or under way.
	bool reading_ahead = false;
	//! Set once, before the batches are pulled.
	idx_t read_ahead_bytes = 0;
	std::function<void()> start_read_ahead;
};

//! Start the single scan of the caller's relation and wrap it in a stream.
//! sequence1, and sequence2 when the pair column is included, must be VARCHAR
//! or BLOB; a VARCHAR is reinterpreted as BLOB in place, not cast.
//!
//! The returned object must be kept alive by the caller until RYpe is finished
//! with the stream — i.e. until after the RYpe output stream has been released.
//! In practice it is stored in the table function's GlobalState and reset in
//! the destructor body, after `output_stream.release` and before the
//! sub-connection whose QueryResult it holds.
//!
//! `id_map` must outlive the stream, for the same reason and by the same means.
unique_ptr<RypeInputStream> BuildRypeInputStream(Connection &conn, RypeIdMap &id_map, RypeInputStreamOptions options);

//! Start the single scan for a filter stream: `columns`, the relation's columns
//! as bind saw them, in their own types, so each group's rows can be handed back
//! whole. sequence1 must be among them and goes to Arrow as `sequence`;
//! sequence2, when present, goes as `pair_sequence`. Both must be VARCHAR or
//! BLOB; a VARCHAR is reinterpreted as BLOB in place, not cast. Of `options`,
//! only relation_quoted, source_name, batch_rows and batch_bytes are read. Same
//! lifetime rule as above.
unique_ptr<RypeInputStream> BuildRypeFilterStream(Connection &conn, const vector<string> &columns,
                                                  RypeInputStreamOptions options);

} // namespace duckdb
