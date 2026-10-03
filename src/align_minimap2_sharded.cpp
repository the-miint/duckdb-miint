#include "align_minimap2_sharded.hpp"
#include "align_common.hpp"
#include "shard_debug.hpp"
#include "shard_progress.hpp"
#include "duckdb/common/file_system.hpp"

namespace duckdb {

// Build minimap2 ShardInfo from raw shard name/counts
// Validates index files exist and are valid minimap2 indexes
static std::vector<ShardInfo> BuildMinimap2ShardInfos(ClientContext &context, const std::string &table_name,
                                                      const std::string &shard_directory, FileSystem &fs) {
	// Get raw shard names and counts from shared utility
	auto raw_shards = ReadShardNameCounts(context, table_name);

	std::vector<ShardInfo> shards;
	shards.reserve(raw_shards.size());

	for (const auto &raw : raw_shards) {
		ShardInfo info;
		info.name = raw.name;
		info.read_count = raw.count;

		// Build index path: shard_directory/shard_name.mmi
		info.index_path = shard_directory;
		if (!info.index_path.empty() && info.index_path.back() != '/') {
			info.index_path += '/';
		}
		info.index_path += info.name + ".mmi";

		// Fail fast: check if .mmi file exists
		if (!fs.FileExists(info.index_path)) {
			throw BinderException("Shard index file does not exist: %s", info.index_path);
		}

		// Validate it's a valid minimap2 index
		if (!miint::Minimap2Aligner::is_index_file(info.index_path)) {
			throw BinderException("File is not a valid minimap2 index: %s", info.index_path);
		}

		shards.push_back(std::move(info));
	}

	return shards;
}

unique_ptr<FunctionData> AlignMinimap2ShardedTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                 vector<LogicalType> &return_types,
                                                                 vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Required: query_table (first positional parameter)
	if (input.inputs.size() < 1) {
		throw BinderException("align_minimap2_sharded requires query_table parameter");
	}
	data->query_table = input.inputs[0].ToString();

	// Required: shard_directory named parameter
	auto shard_dir_param = input.named_parameters.find("shard_directory");
	if (shard_dir_param == input.named_parameters.end() || shard_dir_param->second.IsNull()) {
		throw BinderException("align_minimap2_sharded requires shard_directory parameter");
	}
	data->shard_directory = shard_dir_param->second.ToString();

	// Required: read_to_shard named parameter
	auto read_to_shard_param = input.named_parameters.find("read_to_shard");
	if (read_to_shard_param == input.named_parameters.end() || read_to_shard_param->second.IsNull()) {
		throw BinderException("align_minimap2_sharded requires read_to_shard parameter");
	}
	data->read_to_shard_table = read_to_shard_param->second.ToString();

	// Validate shard_directory exists
	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.DirectoryExists(data->shard_directory)) {
		throw BinderException("Shard directory does not exist: %s", data->shard_directory);
	}

	// Validate query table/view exists. Sharded mode accepts VARCHAR or BIGINT
	// for the query side; the captured id_type drives the output `read_id`
	// column type and how the per-shard read extracts the id column.
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, /*allow_bigint=*/true);

	// Minimap2Aligner::align never reads quality scores. Same as align_minimap2:
	// dropping the flags keeps qual1/qual2 out of the query snapshot and every
	// shard stream, which for long-read FASTQ roughly halves what is held.
	data->query_schema.has_qual1 = false;
	data->query_schema.has_qual2 = false;

	// Validate read_to_shard table schema. Its `read_id` column must share the
	// query table's id type — the strict equality check lets the shard filter in
	// BuildShardReadsSelect compare natively, without implicit casts.
	ValidateReadToShardSchema(context, data->read_to_shard_table, data->query_schema.id_type);

	// Subject side: sharded mode always uses prebuilt .mmi indexes whose
	// subject names are opaque bytes. Output `reference` and `mate_reference`
	// default to VARCHAR — same contract as align_minimap2(index_path:=...).
	data->subject_id_type = LogicalType::VARCHAR;

	// Rebuild output column types with the captured id types. `read_id`
	// mirrors the query side; `reference` / `mate_reference` mirror the
	// subject side (always VARCHAR for prebuilt indexes). The Data() ctor's
	// VARCHAR/VARCHAR placeholder is overwritten here before any caller
	// observes data->types. Must precede the optional shard_name append
	// below — GetAlignmentOutputTypes returns only the 21 alignment columns,
	// so moving this call after the shard_name emplace would clobber it.
	data->types = GetAlignmentOutputTypes(data->query_schema.id_type, data->subject_id_type);

	// Parse minimap2 config parameters (preset, max_secondary, eqx)
	// Always warn about k/w since we use pre-built indexes
	ParseMinimap2ConfigParams(input.named_parameters, data->config, true /* warn_prebuilt_index */);

	// Parse max_threads_per_shard parameter
	auto max_tps_param = input.named_parameters.find("max_threads_per_shard");
	if (max_tps_param != input.named_parameters.end() && !max_tps_param->second.IsNull()) {
		auto val = max_tps_param->second.GetValue<int32_t>();
		if (val < 1 || val > 64) {
			throw BinderException("max_threads_per_shard must be between 1 and 64 (got %d)", val);
		}
		data->max_threads_per_shard = static_cast<idx_t>(val);
	}

	// Parse debug parameter
	auto debug_param = input.named_parameters.find("debug");
	if (debug_param != input.named_parameters.end() && !debug_param->second.IsNull()) {
		data->debug = debug_param->second.GetValue<bool>();
	}

	// Parse progress parameter (opt-in, default false): when true, the function
	// emits clean per-shard progress lines to stderr (see shard_progress.hpp).
	auto progress_param = input.named_parameters.find("progress");
	if (progress_param != input.named_parameters.end() && !progress_param->second.IsNull()) {
		data->progress = progress_param->second.GetValue<bool>();
	}

	// Parse include_shard_name parameter
	auto include_shard_param = input.named_parameters.find("include_shard_name");
	if (include_shard_param != input.named_parameters.end() && !include_shard_param->second.IsNull()) {
		data->include_shard_name = include_shard_param->second.GetValue<bool>();
	}

	// Read shard counts and validate .mmi files exist (fail fast)
	data->shards = BuildMinimap2ShardInfos(context, data->read_to_shard_table, data->shard_directory, fs);

	// Conditionally add shard_name column
	if (data->include_shard_name) {
		data->names.emplace_back("shard_name");
		data->types.emplace_back(LogicalType::VARCHAR);
	}

	// Set output schema
	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> AlignMinimap2ShardedTableFunction::InitGlobal(ClientContext &context,
                                                                                   TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	gstate->shard_count = data.shards.size();
	gstate->max_threads_per_shard = data.max_threads_per_shard;

	// Derive max_active_shards from available threads: ceil(db_threads / max_threads_per_shard)
	// This bounds peak index memory to ceil(threads/tps) * index_size
	idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	idx_t derived = (db_threads + data.max_threads_per_shard - 1) / data.max_threads_per_shard;
	gstate->max_active_shards = std::max<idx_t>(1, std::min(derived, gstate->shard_count));
	gstate->debug = data.debug;
	gstate->progress = data.progress;
	gstate->start_time = std::chrono::steady_clock::now();
	idx_t total = 0;
	for (const auto &shard : data.shards) {
		total += shard.read_count;
	}
	gstate->total_associations.store(total, std::memory_order_relaxed);

	// Read the query and routing relations exactly once each, into snapshots
	// every shard streams from (see GlobalState::snapshots). One copy of the
	// reads, however many shards each is routed to.
	gstate->snapshots = std::make_shared<QuerySnapshots>();
	gstate->snapshots->conn = make_uniq<Connection>(DatabaseInstance::GetDatabase(context));
	InheritTempObjects(context, *gstate->snapshots->conn);
	idx_t snapshot_row_count = 0;
	gstate->snapshots->query_reads =
	    MaterializeQueryReads(*gstate->snapshots->conn, data.query_table, data.query_schema, snapshot_row_count);
	idx_t routing_row_count = 0;
	if (snapshot_row_count == 0) {
		// No reads at all: every shard would load its index only to align nothing.
		// Leaving no shard to claim ends the scan before any index is opened — and
		// with no shard to open a stream, the routing snapshot would be a full copy
		// of a reads x shards relation that nothing ever reads.
		gstate->next_shard_idx = gstate->shard_count;
	} else {
		gstate->snapshots->read_to_shard =
		    MaterializeReadToShard(*gstate->snapshots->conn, data.read_to_shard_table, routing_row_count);
	}
	// This runs on the query's calling thread, never a TaskScheduler worker, so
	// nothing will return the allocation it just freed to the OS on its own —
	// see MakeFreedMemoryFlusher.
	MakeFreedMemoryFlusher(context)();
	SHARD_DBG(*gstate, "InitGlobal: snapshots '%s' (%zu reads) / '%s' (%zu routing rows) materialized",
	          gstate->snapshots->query_reads.c_str(), static_cast<size_t>(snapshot_row_count),
	          gstate->snapshots->read_to_shard.c_str(), static_cast<size_t>(routing_row_count));

	SHARD_DBG_MEM(*gstate, "InitGlobal: shards=%zu db_threads=%zu max_tps=%zu max_active=%zu MaxThreads=%zu",
	              static_cast<size_t>(gstate->shard_count), static_cast<size_t>(db_threads),
	              static_cast<size_t>(gstate->max_threads_per_shard), static_cast<size_t>(gstate->max_active_shards),
	              static_cast<size_t>(gstate->MaxThreads()));
	return gstate;
}

unique_ptr<LocalTableFunctionState>
AlignMinimap2ShardedTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                             GlobalTableFunctionState *global_state) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	// Create per-thread aligner with config
	lstate->aligner = std::make_unique<miint::Minimap2Aligner>(data.config);
	return lstate;
}

// Sub-batch size a shard's stream hands each worker per claim.
static constexpr idx_t SHARD_READ_BATCH_SIZE = 2048;

// Opens a stream over the reads routed to `shard_name`. Every stream of a shard
// (its first part in ClaimWork, each later part in Advance's `prepare`) comes
// through here, so all of them read exactly the same rows.
static std::shared_ptr<QuerySequenceStream>
OpenShardStream(ClientContext &context, const AlignMinimap2ShardedTableFunction::GlobalState &gstate,
                const AlignMinimap2ShardedTableFunction::Data &bind_data, const std::string &shard_name) {
	// The stream holds the snapshots it reads from, so they cannot be dropped
	// while it is open however the scan ends — see QuerySnapshots.
	return std::make_shared<QuerySequenceStream>(
	    context,
	    QuerySequenceStream::SelectSql {BuildShardReadsSelect(gstate.snapshots->query_reads,
	                                                          gstate.snapshots->read_to_shard, bind_data.query_schema,
	                                                          shard_name),
	                                    "reads for shard '" + shard_name + "'"},
	    bind_data.query_schema, SHARD_READ_BATCH_SIZE, gstate.snapshots);
}

std::shared_ptr<ActiveShard> AlignMinimap2ShardedTableFunction::ClaimWork(ClientContext &context, GlobalState &gstate,
                                                                          const Data &bind_data, LocalState &lstate) {
	std::shared_ptr<ActiveShard> active;
	idx_t shard_idx;

	SHARD_DBG(gstate, "ClaimWork: enter");

	{
		std::unique_lock<std::mutex> lock(gstate.lock);

		while (true) {
			// Phase 1: Try to join an existing active shard with capacity
			for (auto &shard : gstate.active_shards) {
				if (shard->ready.load(std::memory_order_acquire) && !shard->exhausted.load(std::memory_order_acquire) &&
				    shard->active_workers.load(std::memory_order_acquire) < gstate.max_threads_per_shard) {
					shard->active_workers.fetch_add(1, std::memory_order_acq_rel);
					auto &info = bind_data.shards[shard->shard_idx];
					SHARD_DBG(gstate, "ClaimWork: JOIN shard %zu '%s' (workers=%zu)",
					          static_cast<size_t>(shard->shard_idx), info.name.c_str(),
					          static_cast<size_t>(shard->active_workers.load(std::memory_order_relaxed)));
					return shard;
				}
			}

			// Phase 2: Try to claim a new shard if under the active shard limit
			bool has_unclaimed = gstate.next_shard_idx < bind_data.shards.size();
			// For small shards (≤ one sub-batch of reads), allow one shard per thread
			// instead of the normal limit, so all threads stay busy
			idx_t max_active = gstate.max_active_shards;
			if (has_unclaimed && bind_data.shards[gstate.next_shard_idx].read_count <= SHARD_READ_BATCH_SIZE) {
				max_active = gstate.max_active_shards * gstate.max_threads_per_shard;
			}
			if (has_unclaimed && gstate.active_shards.size() < max_active) {
				shard_idx = gstate.next_shard_idx++;
				active = std::make_shared<ActiveShard>();
				active->shard_idx = shard_idx;
				active->active_workers.store(1, std::memory_order_release);
				// ready=false (default); set to true after index loaded + stream opened
				gstate.active_shards.push_back(active);
				SHARD_DBG(gstate, "ClaimWork: NEW shard %zu '%s' (active_shards=%zu)", static_cast<size_t>(shard_idx),
				          bind_data.shards[shard_idx].name.c_str(), static_cast<size_t>(gstate.active_shards.size()));
				break; // exit lock to load the index and open the shard's read stream
			}

			// Phase 3: Can't join or start - check if waiting is worthwhile
			bool any_not_exhausted = false;
			for (auto &shard : gstate.active_shards) {
				if (!shard->exhausted.load(std::memory_order_acquire)) {
					any_not_exhausted = true;
					break;
				}
			}

			if (!has_unclaimed && !any_not_exhausted) {
				SHARD_DBG(gstate, "ClaimWork: DONE (no more work)");
				return nullptr; // All shards processed, no active work remaining
			}

			// Wait for: a shard to become ready, capacity to open, or a shard to be removed
			SHARD_DBG(gstate, "ClaimWork: WAIT (active=%zu, unclaimed=%s, any_alive=%s)",
			          static_cast<size_t>(gstate.active_shards.size()), has_unclaimed ? "yes" : "no",
			          any_not_exhausted ? "yes" : "no");
			gstate.cv.wait(lock);
		}
	}
	// Lock released

	// Phase 4: load the index and open this shard's read stream, both OUTSIDE the
	// lock. Either can throw, and either way the cleanup is the same: without it
	// the ActiveShard stays in active_shards with ready=false, so every thread
	// parked in the wait above is never woken and the query hangs instead of
	// failing. Same guard Execute puts around Advance.
	auto &shard_info = bind_data.shards[shard_idx];
	SHARD_DBG(gstate, "ClaimWork: LOADING index '%s'", shard_info.index_path.c_str());
	auto load_start = std::chrono::steady_clock::now();
	try {
		active->parts = std::make_unique<miint::Minimap2PartCursor>(shard_info.index_path, bind_data.config,
		                                                            MakeFreedMemoryFlusher(context));
		auto load_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - load_start)
		        .count();
		SHARD_DBG_MEM(gstate, "ClaimWork: LOADED index shard %zu '%s' (%s) in %ldms", static_cast<size_t>(shard_idx),
		              shard_info.name.c_str(), active->parts->IsMultiPart() ? "multi-part, first part" : "single-part",
		              static_cast<long>(load_ms));
		active->stream = OpenShardStream(context, gstate, bind_data, shard_info.name);
	} catch (...) {
		SHARD_DBG(gstate, "ClaimWork: SETUP FAILED shard %zu", static_cast<size_t>(shard_idx));
		active->exhausted.store(true, std::memory_order_release);
		active->active_workers.fetch_sub(1, std::memory_order_acq_rel);
		{
			std::lock_guard<std::mutex> guard(gstate.lock);
			auto &shards = gstate.active_shards;
			shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
		}
		gstate.cv.notify_all();
		throw;
	}

	// Phase 5: Publish under lock to prevent lost wake-ups with CV
	active->start_time = std::chrono::steady_clock::now();
	{
		std::lock_guard<std::mutex> guard(gstate.lock);
		active->ready.store(true, std::memory_order_release);
	}
	gstate.cv.notify_all();
	if (gstate.progress) {
		// Reads are streamed, so the exact count is only known once the shard is
		// done; the start line reports the read_to_shard rows routed to it.
		shard_progress::Emit("minimap2", shard_progress::FormatShardStart(static_cast<uint64_t>(shard_idx) + 1,
		                                                                  static_cast<uint64_t>(gstate.shard_count),
		                                                                  shard_info.name,
		                                                                  static_cast<int64_t>(shard_info.read_count)));
	}
	return active;
}

void AlignMinimap2ShardedTableFunction::ReleaseWork(GlobalState &gstate, LocalState &lstate) {
	auto active = lstate.current_active_shard;
	lstate.aligner->detach_shared_index();
	lstate.part = miint::Minimap2PartCursor::Attachment {};
	lstate.has_shard = false;

	auto prev_workers = active->active_workers.fetch_sub(1, std::memory_order_acq_rel);
	// prev_workers is the value BEFORE decrement

	SHARD_DBG(gstate, "ReleaseWork: shard %zu (workers %zu->%zu, exhausted=%s)", static_cast<size_t>(active->shard_idx),
	          static_cast<size_t>(prev_workers), static_cast<size_t>(prev_workers - 1),
	          active->exhausted.load(std::memory_order_relaxed) ? "yes" : "no");

	{
		std::lock_guard<std::mutex> guard(gstate.lock);
		if (prev_workers == 1 && active->exhausted.load(std::memory_order_acquire)) {
			// Last worker on an exhausted shard - remove from active list
			auto &shards = gstate.active_shards;
			shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
			SHARD_DBG_MEM(gstate, "ReleaseWork: REMOVED shard %zu (active_shards=%zu)",
			              static_cast<size_t>(active->shard_idx), static_cast<size_t>(shards.size()));
			if (gstate.progress) {
				// The relaxed load sees every worker's relaxed alignments_emitted
				// fetch_add: each fetch_add is sequenced-before that worker's
				// acq_rel active_workers.fetch_sub, and this REMOVE runs only for
				// the last worker (its fetch_sub read 1), whose acquire chains
				// back through the prior decrements' releases — so all fetch_adds
				// happen-before this load. (A miscount would only mis-state this
				// diagnostic line; alignment results are unaffected.)
				const double elapsed_s = std::chrono::duration_cast<std::chrono::duration<double>>(
				                             std::chrono::steady_clock::now() - active->start_time)
				                             .count();
				shard_progress::Emit("minimap2",
				                     shard_progress::FormatShardDone(
				                         static_cast<uint64_t>(active->shard_idx) + 1,
				                         static_cast<uint64_t>(gstate.shard_count), lstate.current_shard_name,
				                         active->stream->RowsDelivered(),
				                         active->alignments_emitted.load(std::memory_order_relaxed), elapsed_s));
			}
		}
		lstate.current_active_shard = nullptr;
		// Notify under lock to prevent lost wake-ups with CV
		gstate.cv.notify_all();
	}
}

void AlignMinimap2ShardedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	while (true) {
		// Check if we have buffered results to output
		idx_t available = local_state.result_buffer.size() - local_state.buffer_offset;

		if (available > 0) {
			// Output up to STANDARD_VECTOR_SIZE results. Id-column types come
			// from the bind data: query side may be VARCHAR or BIGINT; subject
			// side is always VARCHAR for sharded mode (prebuilt .mmi indexes
			// store subject names as opaque bytes).
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, local_state.result_buffer, local_state.buffer_offset, output_count,
			                     bind_data.query_schema.id_type, bind_data.subject_id_type);
			if (bind_data.include_shard_name) {
				auto shard_col_idx = output.ColumnCount() - 1;
				auto &shard_vec = output.data[shard_col_idx];
				for (idx_t i = 0; i < output_count; i++) {
					FlatVector::GetData<string_t>(shard_vec)[i] =
					    StringVector::AddString(shard_vec, local_state.current_shard_name);
				}
			}
			local_state.buffer_offset += output_count;
			return;
		}

		// Buffer is empty, need to get more results

		// Claim a shard if we don't have one
		if (!local_state.has_shard) {
			auto active = ClaimWork(context, global_state, bind_data, local_state);
			if (!active) {
				// No more shards to process
				output.SetCardinality(0);
				return;
			}
			local_state.current_active_shard = active;
			local_state.has_shard = true;
			local_state.current_shard_name = bind_data.shards[active->shard_idx].name;
			// Attached to the shard's current part under the cursor lock below.
		}

		// Pick up the stream of the part this thread is attached to. Attaching and
		// reading `stream` happen together inside the cursor, so the reads this
		// thread drains always belong to the part its aligner holds — see
		// Minimap2PartCursor's WithCurrentPart.
		auto &active = local_state.current_active_shard;
		auto stream =
		    active->parts->WithCurrentPart(local_state.part, *local_state.aligner, [&]() { return active->stream; });

		miint::SequenceRecordBatch query_batch;
		try {
			query_batch = stream->FetchSubBatch();
		} catch (...) {
			// The shard's read query failed mid-stream (e.g. out of memory). Same
			// guard, for the same reason, as the one around Advance below.
			active->exhausted.store(true, std::memory_order_release);
			ReleaseWork(global_state, local_state);
			throw;
		}

		if (query_batch.empty()) {
			// This thread has no work left against the current part. For a
			// multi-part shard, move on to the next part (the leader loads it and
			// opens its stream; everyone else waits) and stream the shard's reads
			// again against it. Only when no part is left is the shard exhausted.
			//
			// The stream has returned empty, so RowsDelivered() is now the exact
			// number of reads this shard aligns against every part.
			const idx_t delivered = stream->RowsDelivered();
			SHARD_DBG(global_state, "Execute: shard %zu current part exhausted after %zu reads",
			          static_cast<size_t>(active->shard_idx), static_cast<size_t>(delivered));
			// Progress was estimated from read_to_shard's row count, which also
			// counts reads absent from query_table. Correct it to the exact count,
			// once per shard (every part streams the same rows).
			if (!active->progress_reconciled.exchange(true, std::memory_order_acq_rel)) {
				const idx_t estimated = bind_data.shards[active->shard_idx].read_count;
				if (delivered > estimated) {
					global_state.total_associations.fetch_add(delivered - estimated, std::memory_order_relaxed);
				} else {
					global_state.total_associations.fetch_sub(estimated - delivered, std::memory_order_relaxed);
				}
			}
			// A shard can legitimately match zero reads (read_to_shard may name
			// read_ids absent from query_table). Without this, Advance walks and
			// fully decodes every remaining part of this shard's index, potentially
			// many GB, to align nothing against them. Mirrors align_minimap2's
			// zero-row guard.
			if (delivered == 0) {
				active->parts->MarkExhausted();
			}
			std::shared_ptr<QuerySequenceStream> next_stream;
			bool advanced;
			try {
				advanced = active->parts->Advance(
				    local_state.part, *local_state.aligner,
				    /*prepare=*/
				    [&]() {
					    next_stream = OpenShardStream(context, global_state, bind_data, local_state.current_shard_name);
				    },
				    /*publish=*/
				    [&]() {
					    active->stream = std::move(next_stream);
					    // Every read is aligned again against the incoming part, so
					    // the denominator has to grow with it, or Progress() saturates
					    // at 100% the moment part 1 finishes and sits there for the
					    // rest of the shard. Part count isn't knowable up front (the
					    // cursor discovers parts as it walks the file), so the estimate
					    // grows as each part is found.
					    global_state.total_associations.fetch_add(delivered, std::memory_order_relaxed);
					    SHARD_DBG_MEM(global_state, "Execute: shard %zu next part loaded, publishing",
					                  static_cast<size_t>(active->shard_idx));
				    });
			} catch (...) {
				// Advance loads an index part and can throw for reasons that are
				// realistic in exactly the low-memory regime this streaming exists
				// for (bad_alloc on a later part, an unnamed sequence, a seek
				// failure), and `prepare` runs a query that can fail. Letting that
				// escape Execute would skip ReleaseWork: this worker's
				// active_workers count would stay up and the shard would never be
				// marked exhausted, so a thread parked in ClaimWork's wait is never
				// woken and the query HANGS instead of failing. Same guard ClaimWork
				// puts around its own index load.
				active->exhausted.store(true, std::memory_order_release);
				ReleaseWork(global_state, local_state);
				throw;
			}
			if (advanced) {
				continue; // re-attach to the newer part and stream from its start
			}
			active->exhausted.store(true, std::memory_order_release);
			ReleaseWork(global_state, local_state);
			continue;
		}

		// Track progress by sequences claimed (before align, so progress updates during I/O)
		global_state.associations_processed.fetch_add(query_batch.size(), std::memory_order_relaxed);
		SHARD_DBG(global_state, "Execute: shard %zu fetched sub-batch (%zu reads)",
		          static_cast<size_t>(active->shard_idx), static_cast<size_t>(query_batch.size()));

		// Align batch — release old buffer capacity before allocating new results
		local_state.result_buffer.clear();
		local_state.result_buffer.shrink_to_fit();
		SHARD_DBG_MEM(global_state, "Execute: shard %zu buffer cleared+shrunk before align",
		              static_cast<size_t>(active->shard_idx));
		local_state.buffer_offset = 0;

		auto align_start = std::chrono::steady_clock::now();
		local_state.aligner->align(query_batch, local_state.result_buffer);
		auto align_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
		        .count();
		// Filter out unmapped reads
		FilterMappedOnly(local_state.result_buffer);
		if (global_state.progress) {
			active->alignments_emitted.fetch_add(local_state.result_buffer.size(), std::memory_order_relaxed);
		}
		SHARD_DBG_MEM(global_state, "Execute: shard %zu ALIGN %zu reads -> %zu results in %ldms",
		              static_cast<size_t>(active->shard_idx), static_cast<size_t>(query_batch.size()),
		              static_cast<size_t>(local_state.result_buffer.size()), static_cast<long>(align_ms));

		// Results (if any) are output on the next iteration; an empty result loops
		// straight back to fetch the next sub-batch.
	}
}

double AlignMinimap2ShardedTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                                   const GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto total = gstate.total_associations.load(std::memory_order_relaxed);
	if (total == 0) {
		return 100.0;
	}
	auto processed = gstate.associations_processed.load(std::memory_order_relaxed);
	return std::min(100.0, 100.0 * static_cast<double>(processed) / static_cast<double>(total));
}

TableFunction AlignMinimap2ShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_minimap2_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	// Named parameters
	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["eqx"] = LogicalType::BOOLEAN;
	tf.named_parameters["max_threads_per_shard"] = LogicalType::INTEGER;
	tf.named_parameters["debug"] = LogicalType::BOOLEAN;
	tf.named_parameters["progress"] = LogicalType::BOOLEAN;
	tf.named_parameters["min_chain_coverage"] = LogicalType::FLOAT;
	tf.named_parameters["include_shard_name"] = LogicalType::BOOLEAN;
	// occ_filter is per-index, so it applies cleanly to each shard independently (#187).
	//
	// include_unmapped is deliberately NOT offered here (#185). A query that finds no chain in
	// shard A routinely maps in shard B, so a per-shard synthetic row would assert "did not
	// align" about a query that did — the opposite of the guarantee the flag exists to give.
	// Doing it correctly needs cross-shard reconciliation, emitting a row only for queries that
	// mapped in no shard at all, which is a global aggregation this per-shard pipeline has no
	// place to hang. Leaving it unregistered makes DuckDB reject the parameter outright rather
	// than silently returning wrong rows.
	tf.named_parameters["occ_filter"] = LogicalType::ANY;

	tf.table_scan_progress = Progress;

	// Alignment output order is non-deterministic — NO_ORDER enables parallel CTAS.
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void AlignMinimap2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
