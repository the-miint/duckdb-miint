#include "align_minimap2.hpp"
#include "align_common.hpp"
#include "shard_debug.hpp"
#include "duckdb/common/allocator.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include <exception>

namespace duckdb {

unique_ptr<FunctionData> AlignMinimap2TableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                          vector<LogicalType> &return_types,
                                                          vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required: query_table (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("align_minimap2 requires query_table parameter");
	}

	data->query_table = input.inputs[0].ToString();

	// Parse subject_table named parameter
	auto subject_param = input.named_parameters.find("subject_table");
	if (subject_param != input.named_parameters.end() && !subject_param->second.IsNull()) {
		data->subject_table = subject_param->second.ToString();
	}

	// Parse index_path named parameter
	auto index_path_param = input.named_parameters.find("index_path");
	if (index_path_param != input.named_parameters.end() && !index_path_param->second.IsNull()) {
		data->index_path = index_path_param->second.ToString();
	}

	// VALIDATION: Exactly one of subject_table or index_path must be provided
	bool has_subject = !data->subject_table.empty();
	bool has_index = data->using_prebuilt_index();

	if (!has_subject && !has_index) {
		throw BinderException("align_minimap2 requires either subject_table or index_path parameter");
	}
	if (has_subject && has_index) {
		throw BinderException(
		    "align_minimap2: Cannot specify both subject_table and index_path. "
		    "Use subject_table to build index from sequences, or index_path to load pre-built index.");
	}

	// Validate query table/view exists (BIGINT read_id is opt-in for PR 1).
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);

	// Minimap2Aligner::align never reads quality scores -- alignment doesn't use
	// them. Dropping the flags here (rather than passing has_qual1/has_qual2
	// through unchanged) means every read of query_table downstream
	// (MaterializeQueryReads' multi-part snapshot and every QuerySequenceStream
	// over query_table) projects only the columns alignment actually consumes.
	// For the multi-part snapshot in particular this roughly halves its size for
	// typical short-read FASTQ input, where qual1/qual2 are comparable in size
	// to sequence1/sequence2.
	data->query_schema.has_qual1 = false;
	data->query_schema.has_qual2 = false;

	// Parse optional named parameters
	auto per_subject_param = input.named_parameters.find("per_subject_database");
	if (per_subject_param != input.named_parameters.end() && !per_subject_param->second.IsNull()) {
		data->per_subject_database = per_subject_param->second.GetValue<bool>();

		// Validate per_subject_database incompatible with index_path
		if (data->per_subject_database && data->using_prebuilt_index()) {
			throw BinderException("per_subject_database mode is incompatible with index_path. "
			                      "Pre-built indexes contain all subjects.");
		}
	}

	// Parse debug parameter
	auto debug_param = input.named_parameters.find("debug");
	if (debug_param != input.named_parameters.end() && !debug_param->second.IsNull()) {
		data->debug = debug_param->second.GetValue<bool>();
	}

	// Parse minimap2 config parameters (preset, max_secondary, k, w, eqx)
	ParseMinimap2ConfigParams(input.named_parameters, data->config, data->using_prebuilt_index());

	// per_subject_database re-aligns every query against every subject in turn, so "this query
	// produced no row" is only ever a statement about ONE subject. Emitting an unmapped row there
	// would assert "did not align" about a query that aligns perfectly well to another subject --
	// with N subjects, a query matching one of them would produce N-1 false negatives. This is the
	// same reason align_minimap2_sharded does not offer the flag at all; here the mode is a runtime
	// choice, so reject the combination rather than return confidently wrong rows.
	if (data->per_subject_database && data->config.include_unmapped) {
		throw BinderException("include_unmapped cannot be combined with per_subject_database: each query is aligned "
		                      "against every subject separately, so an unmapped row would claim the query did not "
		                      "align when it may align to a different subject");
	}

	// Handle subject_table vs index_path modes
	if (data->using_prebuilt_index()) {
		// Validate index file exists
		auto &fs = FileSystem::GetFileSystem(context);
		if (!fs.FileExists(data->index_path)) {
			throw BinderException("Index file does not exist: %s", data->index_path);
		}

		// Validate it's a valid minimap2 index
		// NOTE: This is advisory only (TOCTOU) - actual load may still fail
		if (!miint::Minimap2Aligner::is_index_file(data->index_path)) {
			throw BinderException("File is not a valid minimap2 index: %s", data->index_path);
		}

		// Note: subjects vector remains empty in this mode
		// Subject names in the .mmi file are opaque bytes — default to VARCHAR.
		data->subject_id_type = LogicalType::VARCHAR;
	} else {
		// Traditional mode: validate subject schema, then load subjects.
		// The subject schema's id_type drives the subject-side BIGINT/VARCHAR
		// dispatch in ReadSubjectTable AND the output `reference` /
		// `mate_reference` column types.
		auto subject_schema = ValidateSequenceTableSchema(context, data->subject_table, /*allow_bigint=*/true);
		data->subject_id_type = subject_schema.id_type;
		data->subjects = ReadSubjectTable(context, data->subject_table, subject_schema);
	}

	// Output schema: read_id mirrors the query column type; reference and
	// mate_reference mirror the subject side.
	data->types = GetAlignmentOutputTypes(data->query_schema.id_type, data->subject_id_type);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMinimap2TableFunction::InitGlobal(ClientContext &context,
                                                                            TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->debug = data.debug;
	gstate->start_time = std::chrono::steady_clock::now();

	if (data.per_subject_database) {
		// Per-subject mode: single-threaded, builds index per-subject in Execute()
		gstate->per_subject_mode = true;
		gstate->num_threads = 1;
		auto ps = std::make_unique<PerSubjectModeState>();
		ps->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
		gstate->per_subject = std::move(ps);
		SHARD_DBG(*gstate, "InitGlobal: per_subject_mode, num_threads=1");
	} else if (data.using_prebuilt_index()) {
		// Prebuilt index: open it and read the first part to find out whether
		// it's single- or multi-part. A single-part .mmi (the common case, and
		// every index built by save_minimap2_index) takes the original path
		// unchanged: SharedMinimap2Index for multi-threaded access, no snapshot,
		// no reader retained. A multi-part .mmi (built externally via
		// `minimap2 -I <batch>` smaller than the reference) streams one part at
		// a time so peak memory is one part instead of the whole index — see
		// StandardModeState's doc comment and AdvanceMinimap2Part below.
		auto st = std::make_unique<StandardModeState>();
		std::unique_ptr<miint::Minimap2IndexReader> reader;
		std::shared_ptr<miint::SharedMinimap2Index> first_part;
		bool at_eof = false;
		try {
			reader = std::make_unique<miint::Minimap2IndexReader>(data.index_path, data.config);
			first_part = reader->ReadNextPart();
			if (first_part) {
				// Inside the same try as the load above: AtEof()'s probe read can
				// throw std::runtime_error (fgetpos/fsetpos failure — see its doc
				// comment), and that should get the same IOException wrapping as
				// every other index-load failure here, not escape raw.
				at_eof = reader->AtEof();
			}
		} catch (const std::exception &e) {
			throw IOException("Failed to load minimap2 index from '%s': %s", data.index_path, e.what());
		}
		if (!first_part) {
			throw IOException("Failed to load minimap2 index from '%s': index file contains no parts", data.index_path);
		}

		if (at_eof) {
			st->shared_index = std::move(first_part);
			gstate->standard = std::move(st);
			gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
			SHARD_DBG_MEM(*gstate, "InitGlobal: single-part prebuilt index loaded, MaxThreads()=%zu",
			              static_cast<size_t>(gstate->num_threads));
		} else {
			// include_unmapped cannot be validated at Bind: whether the index is
			// multi-part is only knowable once it's opened here. Same reasoning
			// as the per_subject_database rejection above and the sharded
			// function's outright refusal to offer the parameter (#185): each
			// part is aligned independently, so a per-part "no chain" is only a
			// statement about that part, not the whole reference.
			if (data.config.include_unmapped) {
				throw InvalidInputException(
				    "align_minimap2: include_unmapped cannot be combined with a multi-part index ('%s'). Each part "
				    "is aligned independently, so a per-part unmapped row would claim a query did not align when "
				    "it may align in a later part.",
				    data.index_path);
			}
			st->shared_index = std::move(first_part);
			st->index_reader = std::move(reader);

			// The query relation must be replayed once per part (#229 — see
			// docs/internals/reading-tables-views.md § "Read the relation
			// ONCE"): re-reading it directly would silently drop rows for any
			// relation not stable across re-evaluation. Same snapshot pattern
			// as align_minimap2_sharded's multi-shard case.
			gstate->snapshot_conn = make_uniq<Connection>(DatabaseInstance::GetDatabase(context));
			InheritTempObjects(context, *gstate->snapshot_conn);
			idx_t snapshot_row_count = 0;
			gstate->query_snapshot =
			    MaterializeQueryReads(*gstate->snapshot_conn, data.query_table, data.query_schema, &snapshot_row_count);

			// MaterializeQueryReads reads and frees a corpus-sized amount of memory
			// on THIS thread (the query's calling thread, not one of DuckDB's
			// TaskScheduler worker threads). DuckDB only returns a thread's freed
			// jemalloc pages to the OS from TaskScheduler::ExecuteForever's
			// idle-timeout path (duckdb/src/parallel/task_scheduler.cpp) -- a
			// mechanism this thread never runs through. Without an explicit flush
			// here, the freed memory is correctly accounted as released by
			// DuckDB's own bookkeeping (visible in duckdb_memory()) but stays
			// physically resident, invisible to memory_limit and to the OS/cgroup
			// ceiling this table function is trying to stay under.
			if (Allocator::SupportsFlush()) {
				Allocator::FlushAll();
			}

			// No queries at all means every remaining part would be loaded and
			// decoded (potentially the whole multi-GB index) only to align zero
			// rows against it. part_generation starts at 0 and the first fetch
			// against an empty snapshot returns empty immediately, so marking
			// parts_exhausted now makes that first AdvanceMinimap2Part call
			// return false without touching the reader again.
			if (snapshot_row_count == 0) {
				st->parts_exhausted = true;
			}

			gstate->standard = std::move(st);
			gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
			SHARD_DBG_MEM(
			    *gstate, "InitGlobal: multi-part prebuilt index '%s', snapshot '%s' materialized, MaxThreads()=%zu",
			    data.index_path.c_str(), gstate->query_snapshot.c_str(), static_cast<size_t>(gstate->num_threads));
		}
	} else {
		// Standard mode with subject table: build shared index
		auto st = std::make_unique<StandardModeState>();
		st->shared_index = miint::Minimap2Aligner::BuildSharedIndex(data.subjects, data.config);
		gstate->standard = std::move(st);
		gstate->num_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
		SHARD_DBG_MEM(*gstate, "InitGlobal: shared index built from %zu subjects, MaxThreads()=%zu",
		              data.subjects.size(), static_cast<size_t>(gstate->num_threads));
	}

	// Set up lazy streaming reader for standard mode.
	// Sub-batches are fetched on demand in Execute(), overlapping I/O with alignment.
	if (!gstate->per_subject_mode) {
		if (gstate->standard->index_reader) {
			// Multi-part: replay from the snapshot materialized above.
			gstate->standard->query_stream = std::make_shared<QuerySequenceStream>(
			    *gstate->snapshot_conn, gstate->query_snapshot, data.query_schema);
		} else {
			gstate->standard->query_stream =
			    std::make_shared<QuerySequenceStream>(context, data.query_table, data.query_schema);
		}
		SHARD_DBG(*gstate, "InitGlobal: query stream initialized for lazy sub-batching");
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> AlignMinimap2TableFunction::InitLocal(ExecutionContext &context,
                                                                          TableFunctionInitInput &input,
                                                                          GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();

	if (!gstate.per_subject_mode) {
		// Standard mode: create per-thread aligner
		auto &data = input.bind_data->Cast<Data>();
		lstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
		if (gstate.standard->index_reader) {
			// Multi-part: attach lazily on first real use in
			// ExecuteStandardMultiPart — see LocalState::attached_to_part.
		} else {
			lstate->aligner->attach_shared_index(gstate.standard->shared_index);
			lstate->part_generation = gstate.standard->part_generation;
			lstate->attached_to_part = true;
		}
		auto thread_num = gstate.init_local_count.fetch_add(1) + 1;
		SHARD_DBG(gstate, "InitLocal: thread %zu of %zu initialized", static_cast<size_t>(thread_num),
		          static_cast<size_t>(gstate.num_threads));
	}

	return lstate;
}

// Per-subject mode: single-threaded with global mutex
static void ExecutePerSubject(ClientContext &context, const AlignMinimap2TableFunction::Data &bind_data,
                              AlignMinimap2TableFunction::PerSubjectModeState &ps, DataChunk &output) {
	std::lock_guard<std::mutex> lock(ps.lock);

	if (ps.done) {
		output.SetCardinality(0);
		return;
	}

	idx_t available = ps.result_buffer.size() - ps.buffer_offset;

	while (available == 0) {
		ps.result_buffer.clear();
		ps.buffer_offset = 0;

		// Load all queries into memory on first use.
		//
		// ONE streaming pass, deliberately not LIMIT/OFFSET paging (#229). The
		// paging reader re-issued a fresh query per batch, so any relation that is
		// not stable across re-evaluation — a volatile view, a view over a
		// changing table, a registered Arrow stream — silently returned a
		// *different* row set for the second batch, and a zero-row batch was
		// indistinguishable from end-of-input. Measured: 2000 queries behind a
		// nextval() view aligned ids 1..1024 and 3025..4000, never touching
		// 1025..2000, with no error. Every query is needed here regardless, so a
		// single pass is also strictly less work than paging was.
		if (!ps.queries_loaded) {
			ps.all_queries.is_paired = bind_data.query_schema.has_sequence2;
			QuerySequenceStream stream(context, bind_data.query_table, bind_data.query_schema,
			                           MINIMAP2_QUERY_BATCH_SIZE);
			while (true) {
				auto batch = stream.FetchSubBatch();
				if (batch.empty()) {
					break;
				}
				for (size_t i = 0; i < batch.size(); i++) {
					ps.all_queries.read_ids.push_back(std::move(batch.read_ids[i]));
					ps.all_queries.comments.push_back(std::move(batch.comments[i]));
					ps.all_queries.sequences1.push_back(std::move(batch.sequences1[i]));
					ps.all_queries.quals1.push_back(std::move(batch.quals1[i]));
					if (ps.all_queries.is_paired) {
						ps.all_queries.sequences2.push_back(std::move(batch.sequences2[i]));
						ps.all_queries.quals2.push_back(std::move(batch.quals2[i]));
					}
				}
			}
			ps.queries_loaded = true;
		}

		if (ps.current_subject_idx >= bind_data.subjects.size()) {
			ps.done = true;
			output.SetCardinality(0);
			return;
		}

		ps.aligner->build_single_index(bind_data.subjects[ps.current_subject_idx]);

		if (!ps.all_queries.empty()) {
			ps.aligner->align(ps.all_queries, ps.result_buffer);
		}

		ps.current_subject_idx++;
		available = ps.result_buffer.size() - ps.buffer_offset;
	}

	idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	OutputSAMRecordBatch(output, ps.result_buffer, ps.buffer_offset, output_count, bind_data.query_schema.id_type,
	                     bind_data.subject_id_type);
	ps.buffer_offset += output_count;
}

// Standard mode: multi-threaded with per-thread aligner and lazy sub-batch streaming
static void ExecuteStandard(ClientContext &context, const AlignMinimap2TableFunction::Data &bind_data,
                            AlignMinimap2TableFunction::GlobalState &gstate,
                            AlignMinimap2TableFunction::LocalState &lstate, DataChunk &output) {
	auto &st = *gstate.standard;

	while (true) {
		// 1. If local result_buffer has data, output a chunk
		idx_t available = lstate.result_buffer.size() - lstate.buffer_offset;
		if (available > 0) {
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, lstate.result_buffer, lstate.buffer_offset, output_count,
			                     bind_data.query_schema.id_type, bind_data.subject_id_type);
			lstate.buffer_offset += output_count;
			return;
		}

		// 2. Buffer exhausted — fetch next sub-batch from stream (thread-safe)
		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		auto query_batch = st.query_stream->FetchSubBatch();

		if (query_batch.empty()) {
			SHARD_DBG(gstate, "ExecuteStandard: DONE (stream exhausted)");
			output.SetCardinality(0);
			return;
		}

		SHARD_DBG(gstate, "ExecuteStandard: fetched sub-batch (%zu queries)", query_batch.size());

		// 3. Align using per-thread aligner
		auto align_start = std::chrono::steady_clock::now();
		lstate.aligner->align(query_batch, lstate.result_buffer);
		auto align_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
		        .count();
		SHARD_DBG(gstate, "ExecuteStandard: aligned %zu queries -> %zu results in %ldms", query_batch.size(),
		          lstate.result_buffer.size(), static_cast<long>(align_ms));
	}
}

// Advances StandardModeState past `lstate.part_generation` to the next index
// part, or waits for another thread already doing so.
//
// Called when a thread's fetch from the part it last attached to came back
// empty. Exactly one thread performs the load+swap per transition (the
// `advancing` flag elects a leader); every other thread blocks on `part_cv`.
//
// Only one part is ever resident. That requires dropping BOTH references a
// part can be held by: GlobalState's own `st.shared_index` (reset by the
// leader below) AND every LocalState aligner's shared_ptr copy — including
// this calling thread's. A caller about to become leader or wait must detach
// its own aligner FIRST, before doing either: a thread that instead waited on
// `part_cv` while still attached to the old part (as an earlier version of
// this function did) kept that part alive for the whole load — and with
// enough idle threads sitting in that wait, for the WHOLE REST OF THE QUERY,
// since nothing ever made them detach afterwards either. Measured against a
// human genome multi-part index (2 parts, ~9GB total) with a tiny query
// table and 12 threads: peak RSS was 15.5GB — HIGHER than loading the same
// index single-part (11.75GB) — because only 1-2 of the 12 threads ever did
// real work; the other 10-11 attached to part 1 in InitLocal, got one empty
// fetch, and sat in this function's wait holding that reference for the rest
// of the query, on top of part 2 loading fully. Detaching before waiting (and
// InitLocal no longer attaching eagerly at all — see
// LocalState::attached_to_part) fixes both: an idle thread now holds no part
// reference at all instead of holding the wrong one indefinitely.
//
// Returns true once a part beyond the generation this thread was stuck on
// exists to retry against (whether this thread loaded it or another thread
// did), false once the index reader is exhausted with nothing left to try.
static bool AdvanceMinimap2Part(const AlignMinimap2TableFunction::Data &bind_data,
                                AlignMinimap2TableFunction::GlobalState &gstate,
                                AlignMinimap2TableFunction::LocalState &lstate) {
	auto &st = *gstate.standard;
	idx_t expected_generation = lstate.part_generation;
	lstate.aligner->detach_shared_index();
	lstate.attached_to_part = false;

	std::unique_lock<std::mutex> lock(st.part_lock);

	while (true) {
		if (st.part_generation != expected_generation) {
			// Someone already advanced past the generation this thread was stuck on.
			return true;
		}
		if (st.parts_exhausted) {
			return false;
		}
		if (st.advancing) {
			st.part_cv.wait(lock);
			continue;
		}

		// Become the leader for this transition.
		st.advancing = true;
		st.shared_index.reset(); // free the just-finished part before loading the next
		lock.unlock();

		// Flush this thread's freed jemalloc pages back to the OS before loading
		// the next part. A worker thread driving a multi-part cascade stays
		// continuously busy (fetch/align/advance, part after part) and so never
		// hits the 500ms idle timeout that's the ONLY other place DuckDB flushes
		// a thread's allocator state (TaskScheduler::ExecuteForever) -- without
		// this, every part's freed memory piles up as retained-but-unreturned
		// pages for the whole run instead of actually shrinking peak RSS between
		// parts, even though shared_index.reset() above already dropped the
		// last reference correctly.
		if (Allocator::SupportsFlush()) {
			Allocator::FlushAll();
		}

		// Both the index load AND the query-stream replay happen OUTSIDE the
		// lock: every other thread is either parked on part_cv or about to
		// block on part_lock, so no scheduler thread for this pipeline is free
		// to service a query issued while holding the lock. Both steps are
		// wrapped in the same try/catch so a throw from either one still
		// notifies waiters instead of stranding them — mirrors
		// align_minimap2_sharded's ClaimWork, which guards its equivalent
		// index-load and sequence-prefetch phases the same way.
		std::shared_ptr<miint::SharedMinimap2Index> next_index;
		std::shared_ptr<QuerySequenceStream> next_stream;
		std::exception_ptr load_error;
		try {
			next_index = st.index_reader->ReadNextPart();
			if (next_index) {
				next_stream = std::make_shared<QuerySequenceStream>(*gstate.snapshot_conn, gstate.query_snapshot,
				                                                    bind_data.query_schema);
			}
		} catch (...) {
			load_error = std::current_exception();
		}

		lock.lock();
		st.advancing = false;
		if (load_error) {
			// Stop every other thread too — the reader is now in an unknown
			// state and cannot be trusted for a retry.
			st.parts_exhausted = true;
			st.part_cv.notify_all();
			lock.unlock();
			std::rethrow_exception(load_error);
		}
		if (!next_index) {
			st.parts_exhausted = true;
			st.part_cv.notify_all();
			return false;
		}

		st.shared_index = std::move(next_index);
		st.query_stream = std::move(next_stream);
		st.part_generation++;
		SHARD_DBG_MEM(gstate, "AdvanceMinimap2Part: advanced to part_generation=%zu",
		              static_cast<size_t>(st.part_generation));
		st.part_cv.notify_all();
		return true;
	}
}

// Multi-part prebuilt index: like ExecuteStandard, but when this thread's
// local view of query_stream is exhausted, it coordinates with the other
// worker threads (via AdvanceMinimap2Part) to move to the next index part
// before looping back, rather than treating exhaustion as end-of-results.
static void ExecuteStandardMultiPart(const AlignMinimap2TableFunction::Data &bind_data,
                                     AlignMinimap2TableFunction::GlobalState &gstate,
                                     AlignMinimap2TableFunction::LocalState &lstate, DataChunk &output) {
	auto &st = *gstate.standard;

	while (true) {
		// 1. If local result_buffer has data, output a chunk
		idx_t available = lstate.result_buffer.size() - lstate.buffer_offset;
		if (available > 0) {
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, lstate.result_buffer, lstate.buffer_offset, output_count,
			                     bind_data.query_schema.id_type, bind_data.subject_id_type);
			lstate.buffer_offset += output_count;
			return;
		}

		// 2. Make sure this thread is attached to the current part before fetching.
		// !attached_to_part covers both a fresh thread's lazy first attach and a
		// thread returning from AdvanceMinimap2Part, which always detaches first.
		std::shared_ptr<QuerySequenceStream> current_stream;
		{
			std::lock_guard<std::mutex> lock(st.part_lock);
			if (!lstate.attached_to_part || lstate.part_generation != st.part_generation) {
				// st.shared_index is only null for the instant between a leader
				// resetting it (line ~423) and installing the next part; a thread
				// that reaches here in that window has not yet observed the
				// generation bump, so today it always finds the OLD query_stream
				// already exhausted and immediately loops back into
				// AdvanceMinimap2Part instead of calling align(). Asserted rather
				// than silently tolerated so a future change to that ordering
				// fails loud instead of surfacing as "No index built" mid-query.
				D_ASSERT(st.shared_index != nullptr);
				lstate.aligner->attach_shared_index(st.shared_index);
				lstate.part_generation = st.part_generation;
				lstate.attached_to_part = true;
			}
			current_stream = st.query_stream;
		}

		// 3. Buffer exhausted — fetch next sub-batch from stream (thread-safe)
		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		auto query_batch = current_stream->FetchSubBatch();

		if (!query_batch.empty()) {
			SHARD_DBG(gstate, "ExecuteStandardMultiPart: fetched sub-batch (%zu queries)", query_batch.size());
			auto align_start = std::chrono::steady_clock::now();
			lstate.aligner->align(query_batch, lstate.result_buffer);
			auto align_ms =
			    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
			        .count();
			SHARD_DBG(gstate, "ExecuteStandardMultiPart: aligned %zu queries -> %zu results in %ldms",
			          query_batch.size(), lstate.result_buffer.size(), static_cast<long>(align_ms));
			continue; // loop back to step 1 to output
		}

		// 4. This thread's view of the current part is exhausted.
		if (!AdvanceMinimap2Part(bind_data, gstate, lstate)) {
			SHARD_DBG(gstate, "ExecuteStandardMultiPart: DONE (all parts exhausted)");
			output.SetCardinality(0);
			return;
		}
		// A part beyond lstate.part_generation now exists (whether this thread
		// loaded it or another thread did) — loop back to re-attach and retry.
	}
}

void AlignMinimap2TableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.per_subject_mode) {
		ExecutePerSubject(context, bind_data, *gstate.per_subject, output);
	} else {
		auto &lstate = data_p.local_state->Cast<LocalState>();
		if (gstate.standard->index_reader) {
			ExecuteStandardMultiPart(bind_data, gstate, lstate, output);
		} else {
			ExecuteStandard(context, bind_data, gstate, lstate, output);
		}
	}
}

TableFunction AlignMinimap2TableFunction::GetFunction() {
	// Only query_table is a positional parameter
	auto tf = TableFunction("align_minimap2", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["subject_table"] = LogicalType::VARCHAR;
	tf.named_parameters["index_path"] = LogicalType::VARCHAR;
	tf.named_parameters["per_subject_database"] = LogicalType::BOOLEAN;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["k"] = LogicalType::INTEGER;
	tf.named_parameters["w"] = LogicalType::INTEGER;
	tf.named_parameters["eqx"] = LogicalType::BOOLEAN;
	tf.named_parameters["debug"] = LogicalType::BOOLEAN;
	tf.named_parameters["min_chain_coverage"] = LogicalType::FLOAT;
	// ANY so both `occ_filter := 0` and minimap2's two-value `occ_filter := '1000,5000'` bind;
	// the value is read via ToString() either way. Same approach as read_alignments'
	// reference_lengths.
	tf.named_parameters["occ_filter"] = LogicalType::ANY;
	tf.named_parameters["include_unmapped"] = LogicalType::BOOLEAN;

	// Alignment output order is non-deterministic (depends on thread scheduling),
	// so NO_ORDER lets DuckDB parallelize CTAS pipelines instead of serializing
	// them via preserve_insertion_order.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignMinimap2TableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
