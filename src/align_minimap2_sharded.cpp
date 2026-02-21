#include "align_minimap2_sharded.hpp"
#include "align_common.hpp"
#include "duckdb/common/file_system.hpp"

#include <cstdio>

namespace duckdb {

// Read current process RSS from /proc/self/status (Linux only, returns MB)
static long GetRSSMB() {
	FILE *f = fopen("/proc/self/status", "r");
	if (!f) {
		return -1;
	}
	long rss_kb = -1;
	char line[256];
	while (fgets(line, sizeof(line), f)) {
		if (strncmp(line, "VmRSS:", 6) == 0) {
			sscanf(line + 6, "%ld", &rss_kb);
			break;
		}
	}
	fclose(f);
	return rss_kb > 0 ? rss_kb / 1024 : -1;
}

// NOLINTNEXTLINE: macro for debug logging with elapsed time and thread ID
#define SHARD_DBG(gstate, ...)                                                                                         \
	do {                                                                                                               \
		if ((gstate).debug) {                                                                                          \
			auto _elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -   \
			                                                                      (gstate).start_time);                \
			auto _tid = std::hash<std::thread::id> {}(std::this_thread::get_id()) % 10000;                             \
			fprintf(stderr, "[%7ldms T%04zu] ", static_cast<long>(_elapsed.count()), _tid);                            \
			fprintf(stderr, __VA_ARGS__);                                                                              \
			fprintf(stderr, "\n");                                                                                     \
		}                                                                                                              \
	} while (0)

// Like SHARD_DBG but appends current RSS
#define SHARD_DBG_MEM(gstate, ...)                                                                                     \
	do {                                                                                                               \
		if ((gstate).debug) {                                                                                          \
			auto _elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() -   \
			                                                                      (gstate).start_time);                \
			auto _tid = std::hash<std::thread::id> {}(std::this_thread::get_id()) % 10000;                             \
			fprintf(stderr, "[%7ldms T%04zu] ", static_cast<long>(_elapsed.count()), _tid);                            \
			fprintf(stderr, __VA_ARGS__);                                                                              \
			fprintf(stderr, " [RSS=%ldMB]\n", GetRSSMB());                                                             \
		}                                                                                                              \
	} while (0)

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

	// Validate query table/view exists
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table, true /* allow_paired */);

	// Validate read_to_shard table schema
	ValidateReadToShardSchema(context, data->read_to_shard_table);

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

	// Read shard counts and validate .mmi files exist (fail fast)
	data->shards = BuildMinimap2ShardInfos(context, data->read_to_shard_table, data->shard_directory, fs);

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
	gstate->start_time = std::chrono::steady_clock::now();
	idx_t total = 0;
	for (const auto &shard : data.shards) {
		total += shard.read_count;
	}
	gstate->total_associations = total;
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

std::shared_ptr<ActiveShard> AlignMinimap2ShardedTableFunction::ClaimWork(GlobalState &gstate, const Data &bind_data,
                                                                          LocalState &lstate) {
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
			if (has_unclaimed && gstate.active_shards.size() < gstate.max_active_shards) {
				shard_idx = gstate.next_shard_idx++;
				active = std::make_shared<ActiveShard>();
				active->shard_idx = shard_idx;
				// Compute per-shard batch size: divide reads evenly across threads
				auto &si = bind_data.shards[shard_idx];
				active->batch_size = std::max<idx_t>(1, si.read_count / gstate.max_threads_per_shard);
				active->active_workers.store(1, std::memory_order_release);
				// ready=false (default); set to true after index loads
				gstate.active_shards.push_back(active);
				SHARD_DBG(gstate, "ClaimWork: NEW shard %zu '%s' (active_shards=%zu, batch_size=%zu)",
				          static_cast<size_t>(shard_idx), si.name.c_str(),
				          static_cast<size_t>(gstate.active_shards.size()), static_cast<size_t>(active->batch_size));
				break; // exit lock to load index
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

	// Phase 4: Load index OUTSIDE lock
	auto &shard_info = bind_data.shards[shard_idx];
	SHARD_DBG(gstate, "ClaimWork: LOADING index '%s'", shard_info.index_path.c_str());
	auto load_start = std::chrono::steady_clock::now();
	try {
		auto shared_idx = std::make_shared<miint::SharedMinimap2Index>(shard_info.index_path, bind_data.config);
		active->index = std::move(shared_idx);
	} catch (...) {
		// Remove failed shard from active list and notify waiters
		SHARD_DBG(gstate, "ClaimWork: LOAD FAILED shard %zu", static_cast<size_t>(shard_idx));
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
	auto load_ms =
	    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - load_start).count();
	SHARD_DBG_MEM(gstate, "ClaimWork: LOADED shard %zu '%s' in %ldms", static_cast<size_t>(shard_idx),
	              shard_info.name.c_str(), static_cast<long>(load_ms));

	// Phase 5: Publish under lock to prevent lost wake-ups with CV
	{
		std::lock_guard<std::mutex> guard(gstate.lock);
		active->ready.store(true, std::memory_order_release);
	}
	gstate.cv.notify_all();
	return active;
}

void AlignMinimap2ShardedTableFunction::ReleaseWork(GlobalState &gstate, LocalState &lstate) {
	auto active = lstate.current_active_shard;
	lstate.aligner->detach_shared_index();
	lstate.has_shard = false;

	auto prev_workers = active->active_workers.fetch_sub(1, std::memory_order_acq_rel);
	// prev_workers is the value BEFORE decrement

	SHARD_DBG(gstate, "ReleaseWork: shard %zu (workers %zu->%zu, exhausted=%s)", static_cast<size_t>(active->shard_idx),
	          static_cast<size_t>(prev_workers), static_cast<size_t>(prev_workers - 1),
	          active->exhausted.load(std::memory_order_relaxed) ? "yes" : "no");

	if (prev_workers == 1 && active->exhausted.load(std::memory_order_acquire)) {
		// Last worker on an exhausted shard - remove from active list
		std::lock_guard<std::mutex> guard(gstate.lock);
		auto &shards = gstate.active_shards;
		shards.erase(std::remove(shards.begin(), shards.end(), active), shards.end());
		SHARD_DBG_MEM(gstate, "ReleaseWork: REMOVED shard %zu (active_shards=%zu)",
		              static_cast<size_t>(active->shard_idx), static_cast<size_t>(shards.size()));
	}

	lstate.current_active_shard = nullptr;

	// Always notify: capacity may have opened or a shard slot freed
	gstate.cv.notify_all();
}

void AlignMinimap2ShardedTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	while (true) {
		// Check if we have buffered results to output
		idx_t available = local_state.result_buffer.size() - local_state.buffer_offset;

		if (available > 0) {
			// Output up to STANDARD_VECTOR_SIZE results
			idx_t output_count = std::min(available, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSAMRecordBatch(output, local_state.result_buffer, local_state.buffer_offset, output_count);
			local_state.buffer_offset += output_count;
			return;
		}

		// Buffer is empty, need to get more results

		// Claim a shard if we don't have one
		if (!local_state.has_shard) {
			auto active = ClaimWork(global_state, bind_data, local_state);
			if (!active) {
				// No more shards to process
				output.SetCardinality(0);
				return;
			}
			local_state.current_active_shard = active;
			local_state.has_shard = true;
			local_state.aligner->attach_shared_index(active->index);
		}

		// Atomically claim batch offset
		auto &active = local_state.current_active_shard;
		auto &shard_info = bind_data.shards[active->shard_idx];
		idx_t my_offset = active->next_batch_offset.fetch_add(active->batch_size, std::memory_order_acq_rel);

		if (my_offset >= shard_info.read_count) {
			// All batches claimed already
			SHARD_DBG(global_state, "Execute: shard %zu offset %zu >= read_count %zu, exhausted",
			          static_cast<size_t>(active->shard_idx), static_cast<size_t>(my_offset),
			          static_cast<size_t>(shard_info.read_count));
			active->exhausted.store(true, std::memory_order_release);
			ReleaseWork(global_state, local_state);
			continue;
		}

		// Read query batch at claimed offset.
		// ReadShardQueryBatch takes offset by reference and updates it; we use a
		// throwaway copy because the atomic next_batch_offset drives pagination.
		idx_t batch_offset = my_offset;
		miint::SequenceRecordBatch query_batch;
		SHARD_DBG(global_state,
		          "Execute: shard %zu READ offset=%zu batch_size=%zu (query: SELECT ... FROM %s q JOIN %s r "
		          "ON q.read_id=r.read_id WHERE r.shard_name='%s' ORDER BY %s LIMIT %zu OFFSET %zu)",
		          static_cast<size_t>(active->shard_idx), static_cast<size_t>(my_offset),
		          static_cast<size_t>(active->batch_size), bind_data.query_table.c_str(),
		          bind_data.read_to_shard_table.c_str(), shard_info.name.c_str(),
		          bind_data.query_schema.is_physical_table ? "q.rowid" : "q.read_id",
		          static_cast<size_t>(active->batch_size), static_cast<size_t>(batch_offset));
		auto read_start = std::chrono::steady_clock::now();
		bool has_more =
		    ReadShardQueryBatch(context, bind_data.query_table, bind_data.read_to_shard_table, shard_info.name,
		                        bind_data.query_schema, active->batch_size, batch_offset, query_batch);
		auto read_ms =
		    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - read_start)
		        .count();
		SHARD_DBG_MEM(global_state, "Execute: shard %zu READ done: %zu reads in %ldms (has_more=%s)",
		              static_cast<size_t>(active->shard_idx), static_cast<size_t>(query_batch.size()),
		              static_cast<long>(read_ms), has_more ? "yes" : "no");

		if (query_batch.empty() && !has_more) {
			// Actual data exhausted (read_count was approximate)
			active->exhausted.store(true, std::memory_order_release);
			ReleaseWork(global_state, local_state);
			continue;
		}

		if (!has_more) {
			// Last batch - mark exhausted so no new workers join
			active->exhausted.store(true, std::memory_order_release);
		}

		// Align batch — release old buffer capacity before allocating new results
		local_state.result_buffer.clear();
		local_state.result_buffer.shrink_to_fit();
		SHARD_DBG_MEM(global_state, "Execute: shard %zu buffer cleared+shrunk before align",
		              static_cast<size_t>(active->shard_idx));
		local_state.buffer_offset = 0;

		if (!query_batch.empty()) {
			auto align_start = std::chrono::steady_clock::now();
			local_state.aligner->align(query_batch, local_state.result_buffer);
			auto align_ms =
			    std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::steady_clock::now() - align_start)
			        .count();
			// Filter out unmapped reads
			FilterMappedOnly(local_state.result_buffer);
			SHARD_DBG_MEM(global_state, "Execute: shard %zu ALIGN %zu reads -> %zu results in %ldms",
			              static_cast<size_t>(active->shard_idx), static_cast<size_t>(query_batch.size()),
			              static_cast<size_t>(local_state.result_buffer.size()), static_cast<long>(align_ms));
			// Track progress by associations processed
			global_state.associations_processed.fetch_add(query_batch.size(), std::memory_order_relaxed);
		}

		// If this was the last batch and we got results, we'll output them in the next iteration.
		// If no results, we'll loop back and ReleaseWork will happen when we try to claim the next batch.
	}
}

double AlignMinimap2ShardedTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                                   const GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	if (gstate.total_associations == 0) {
		return 100.0;
	}
	return 100.0 * static_cast<double>(gstate.associations_processed.load(std::memory_order_relaxed)) /
	       static_cast<double>(gstate.total_associations);
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

	tf.table_scan_progress = Progress;

	return tf;
}

void AlignMinimap2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
