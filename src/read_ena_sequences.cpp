#include "read_ena_sequences.hpp"
#include "SequenceRecord.hpp"
#include "ena_resolver_cache.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/vector_size.hpp"
#include <cerrno>
#include <fstream>

namespace duckdb {

// ---- Data ----

ReadENASequencesTableFunction::Data::Data(std::vector<miint::ENARunInfo> runs, bool include_fp, uint8_t offset,
                                          bool trim, const std::string &prefer_format)
    : runs(std::move(runs)), include_filepath(include_fp), qual_offset(offset), use_aspera(false), trim(trim),
      prefer_format(prefer_format), names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1",
                                           "qual2", "run_accession", "sample_accession", "experiment_accession"}),
      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT), LogicalType::LIST(LogicalType::UTINYINT),
             LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}) {
	if (include_filepath) {
		names.emplace_back("filepath");
		types.emplace_back(LogicalType::VARCHAR);
	}
}

// ---- GlobalState ----

ReadENASequencesTableFunction::GlobalState::GlobalState(FileSystem &fs, const std::vector<miint::ENARunInfo> &runs,
                                                        bool use_aspera, bool trim)
    : runs(runs), next_run_idx(0), fs(fs), use_aspera(use_aspera), trim(trim), runs_completed(0),
      total_runs(runs.size()), bytes_completed(0), skipped_warned(false) {
	uint64_t sum_bytes = 0;
	for (const auto &r : runs) {
		sum_bytes += r.total_bytes;
	}
	total_bytes = sum_bytes;
	readers.resize(runs.size());
	run_sequence_counters.resize(runs.size(), 1);
}

// ---- Bind ----

// Resolve accessions to run info via ENA Portal API (batched + cached).
// Preserves input order: multi-run accessions (e.g., studies) are flattened in the
// order ENA returns them; across distinct accessions the caller's input order is kept.
//
// Cache lifetime: the cache is constructed locally and lives only for this call.
// For Phase B (lateral / in_out_function) we need the cache to survive across
// per-outer-row ExecuteInOut invocations — promote ownership to the bind Data
// (or GlobalState) at that point.
static std::vector<miint::ENARunInfo> ResolveRuns(miint::ENAClient &client, const std::vector<std::string> &accessions,
                                                  const std::string &prefer_format) {
	miint::ENAResolverCache cache(256);
	miint::ENAFetcher fetcher = [&client](const std::string &url) {
		return client.FetchURL(url);
	};
	auto resolved = miint::ResolveRunsBatch(fetcher, cache, accessions, prefer_format);

	std::vector<miint::ENARunInfo> runs;
	for (const auto &acc : accessions) {
		auto it = resolved.find(acc);
		if (it == resolved.end()) {
			continue;
		}
		for (const auto &info : it->second) {
			runs.push_back(info);
		}
	}
	return runs;
}

unique_ptr<FunctionData> ReadENASequencesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	// Lateral / subquery dispatch: DuckDB's BindTableInTableOutFunction binds
	// the expressions into a subquery and leaves `input.inputs` empty. An empty
	// inputs vector is therefore the authoritative signal — it does not mean
	// the user called read_ena_sequences() with no args (the function binder
	// rejects that earlier because the overload has exactly one positional
	// argument). A literal NULL is caught separately below with a clear error.
	const bool deferred = input.inputs.empty();

	std::vector<std::string> accessions;
	if (!deferred) {
		if (input.inputs[0].IsNull()) {
			throw InvalidInputException(
			    "read_ena_sequences: accession cannot be NULL (use a VARCHAR literal, VARCHAR[] list, "
			    "or a correlated column / subquery for lateral mode)");
		}
		if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
			accessions.push_back(input.inputs[0].ToString());
		} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
			auto &list_children = ListValue::GetChildren(input.inputs[0]);
			for (const auto &child : list_children) {
				accessions.push_back(child.ToString());
			}
		} else {
			throw InvalidInputException("read_ena_sequences: first argument must be VARCHAR or VARCHAR[]");
		}

		if (accessions.empty()) {
			throw InvalidInputException("read_ena_sequences: at least one accession must be provided");
		}
		for (auto &acc : accessions) {
			// Trim whitespace
			auto start = acc.find_first_not_of(" \t\n\r");
			if (start == std::string::npos) {
				acc.clear();
			} else {
				acc = acc.substr(start, acc.find_last_not_of(" \t\n\r") - start + 1);
			}
			if (acc.empty()) {
				throw InvalidInputException("read_ena_sequences: accession cannot be empty");
			}
		}
	}

	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	uint8_t qual_offset = 33;
	auto qo_param = input.named_parameters.find("qual_offset");
	if (qo_param != input.named_parameters.end() && !qo_param->second.IsNull()) {
		qual_offset = static_cast<uint8_t>(qo_param->second.GetValue<int64_t>());
	}

	std::string download_method = "auto";
	auto dm_param = input.named_parameters.find("download_method");
	if (dm_param != input.named_parameters.end() && !dm_param->second.IsNull()) {
		download_method = dm_param->second.ToString();
		if (download_method != "auto" && download_method != "aspera" && download_method != "http") {
			throw InvalidInputException(
			    "read_ena_sequences: download_method must be 'auto', 'aspera', or 'http' (got '%s')", download_method);
		}
	}

	std::string prefer_format = "auto";
	auto pf_param = input.named_parameters.find("prefer_format");
	if (pf_param != input.named_parameters.end() && !pf_param->second.IsNull()) {
		prefer_format = pf_param->second.ToString();
		if (prefer_format != "auto" && prefer_format != "fastq" && prefer_format != "sff") {
			throw InvalidInputException(
			    "read_ena_sequences: prefer_format must be 'auto', 'fastq', or 'sff' (got '%s')", prefer_format);
		}
	}

	// Named as `trim_sff` rather than `trim` to avoid a parser/binder clash: TRIM
	// is a SQL function keyword, and read_ena_sequences is a dual-path function
	// (`function` + `in_out_function`) where DuckDB's binder routes any non-
	// scalar argument — including `name=value` syntax — through the in-out path,
	// which skips named-parameter extraction. `trim_sff => false` survives both
	// issues. Internal identifiers keep the shorter `trim` name.
	bool trim = true;
	auto trim_param = input.named_parameters.find("trim_sff");
	if (trim_param != input.named_parameters.end() && !trim_param->second.IsNull()) {
		trim = trim_param->second.GetValue<bool>();
	}

	// `max_sequences`: 0 / NULL / absent = unlimited. Negative values are a
	// user error caught here rather than silently coerced, matching the
	// one-named-parameter-rejects-garbage precedent of qual_offset /
	// download_method / prefer_format above.
	uint64_t max_sequences = 0;
	auto ms_param = input.named_parameters.find("max_sequences");
	if (ms_param != input.named_parameters.end() && !ms_param->second.IsNull()) {
		int64_t raw = ms_param->second.GetValue<int64_t>();
		if (raw < 0) {
			throw InvalidInputException(
			    "read_ena_sequences: max_sequences must be >= 0 (got %lld; use 0 for unlimited)",
			    static_cast<long long>(raw));
		}
		max_sequences = static_cast<uint64_t>(raw);
	}

	auto &db = DatabaseInstance::GetDatabase(context);

	std::vector<miint::ENARunInfo> runs;
	if (!deferred) {
		// Resolve accessions to run info via ENA Portal API (batched + cached)
		miint::ENAClient client(db);
		runs = ResolveRuns(client, accessions, prefer_format);

		if (runs.empty()) {
			throw IOException("read_ena_sequences: no sequence data found for the provided accession(s)");
		}
	}

	auto data = make_uniq<Data>(std::move(runs), include_filepath, qual_offset, trim, prefer_format);
	data->deferred_resolution = deferred;
	data->db_ptr = &db;
	data->max_sequences = max_sequences;
	if (deferred) {
		// Per-bind cache + ENAClient: metadata lookups across outer rows dedupe
		// and the ENAClient's rate limiter (~3 req/sec) throttles globally.
		// Capacity 256 matches the eager-resolve ResolveRuns helper.
		data->resolver_cache = std::make_unique<miint::ENAResolverCache>(256);
		data->lateral_client = std::make_unique<miint::ENAClient>(db);
	}

	// Determine whether to use Aspera (only relevant for FASTX runs; SFF uses HTTP only)
	data->use_aspera = false;
#if MIINT_ASPERA_SUPPORTED
	if (deferred) {
		// Lateral path: runs are resolved per outer row, so we can't pre-validate
		// that every run has aspera paths. Fall back to HTTP streaming regardless
		// of download_method; the lateral use case (LIMIT-inside-LATERAL) benefits
		// from short-circuiting, not throughput. Explicit aspera request errors.
		if (download_method == "aspera") {
			throw IOException("read_ena_sequences: download_method='aspera' is not supported with lateral / "
			                  "correlated arguments; use 'http' or omit the parameter");
		}
	} else {
		// Check if there are any FASTX runs that could benefit from Aspera
		bool has_fastx_runs = false;
		for (const auto &run : data->runs) {
			if (run.format == miint::ENASequenceFormat::FASTX) {
				has_fastx_runs = true;
				break;
			}
		}
		if (has_fastx_runs && (download_method == "aspera" || download_method == "auto")) {
			bool aspera_available = miint::AsperaUtils::IsAvailable();
			if (download_method == "aspera" && !aspera_available) {
				throw IOException("read_ena_sequences: download_method='aspera' but ascp not found in PATH. "
				                  "Install IBM Aspera CLI or use download_method='http'");
			}
			if (aspera_available) {
				// Check that all FASTX runs have aspera paths
				bool all_have_aspera = true;
				for (const auto &run : data->runs) {
					if (run.format == miint::ENASequenceFormat::FASTX && run.aspera_paths.empty()) {
						all_have_aspera = false;
						break;
					}
				}
				if (download_method == "aspera" && !all_have_aspera) {
					throw IOException(
					    "read_ena_sequences: download_method='aspera' but not all runs have Aspera paths");
				}
				if (all_have_aspera) {
					std::string ascp_path = miint::AsperaUtils::FindAscp();
					std::string key_path = miint::AsperaUtils::ResolveKey(db, ascp_path, download_method == "aspera");
					if (!key_path.empty()) {
						data->use_aspera = true;
						data->aspera_config = miint::AsperaUtils::BuildConfig(ascp_path, key_path);
					}
				}
			}
		}
	}
#else
	if (download_method == "aspera") {
		throw IOException("read_ena_sequences: Aspera is not supported on this platform");
	}
#endif

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENASequencesTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	FileSystem &fs = FileSystem::GetFileSystem(context);
	auto state = make_uniq<GlobalState>(fs, data.runs, data.use_aspera, data.trim);

#if MIINT_ASPERA_SUPPORTED
	if (data.use_aspera) {
		state->aspera_config = data.aspera_config;
	}
#endif

	return state;
}

unique_ptr<LocalTableFunctionState> ReadENASequencesTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENASequencesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	// Helper: construct a fresh PerRunReader for a given run index.
	// All per-run state (temp files, ascp process, FASTX/SFF reader) lives on
	// the PerRunReader; the GlobalState holds only the run metadata + counters.
	auto make_reader = [&](size_t idx) {
		const miint::AsperaConfig *cfg = nullptr;
#if MIINT_ASPERA_SUPPORTED
		if (global_state.use_aspera) {
			cfg = &global_state.aspera_config;
		}
#endif
		return std::make_unique<miint::PerRunReader>(global_state.fs, global_state.runs[idx], global_state.use_aspera,
		                                             global_state.trim, global_state.open_mutex, cfg,
		                                             bind_data.max_sequences);
	};

	miint::SequenceRecordBatch batch;
	size_t current_run_idx = 0;

	// Claim-read-release loop (same pattern as read_fastx)
	while (true) {
		if (!local_state.has_run) {
			bool all_claimed = false;
			{
				lock_guard<mutex> read_lock(global_state.lock);
				if (global_state.next_run_idx >= global_state.runs.size()) {
					all_claimed = true;
				} else {
					local_state.current_run_idx = global_state.next_run_idx;
					global_state.next_run_idx++;
					local_state.has_run = true;
				}
			} // Lock released before I/O or warning output

			if (all_claimed) {
				// Print skipped-runs warning exactly once across all threads
				if (!global_state.skipped_warned.exchange(true, std::memory_order_acq_rel)) {
					lock_guard<mutex> skip_guard(global_state.skipped_lock);
					if (!global_state.skipped_runs.empty()) {
						std::string msg =
						    "read_ena_sequences: WARNING: " + std::to_string(global_state.skipped_runs.size()) +
						    " run(s) skipped due to download errors:";
						for (const auto &acc : global_state.skipped_runs) {
							msg += " " + acc;
						}
						Printer::Print(msg);
					}
				}
				output.SetCardinality(0);
				return;
			}

			// Open reader — retry once on failure before skipping.
			// On retry, destroy the partially-initialized reader (its dtor runs
			// any teardown) and recreate from scratch. Reset the run's
			// sequence counter so the retry starts numbering from 1.
			try {
				global_state.readers[local_state.current_run_idx] = make_reader(local_state.current_run_idx);
				global_state.readers[local_state.current_run_idx]->Open();
			} catch (const std::exception &e) {
				auto &run = global_state.runs[local_state.current_run_idx];
				Printer::PrintF("read_ena_sequences: warning: run '%s' failed to open (%s), retrying...",
				                run.run_accession, e.what());
				global_state.readers[local_state.current_run_idx].reset();
				global_state.run_sequence_counters[local_state.current_run_idx] = 1;
				try {
					global_state.readers[local_state.current_run_idx] = make_reader(local_state.current_run_idx);
					global_state.readers[local_state.current_run_idx]->Open();
				} catch (const std::exception &e2) {
					Printer::PrintF("read_ena_sequences: warning: run '%s' failed on retry (%s), skipping",
					                run.run_accession, e2.what());
					global_state.readers[local_state.current_run_idx].reset();
					{
						lock_guard<mutex> skip_guard(global_state.skipped_lock);
						global_state.skipped_runs.push_back(run.run_accession);
					}
					global_state.bytes_completed.fetch_add(run.total_bytes, std::memory_order_relaxed);
					global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
					local_state.has_run = false;
					continue;
				}
			}
		}

		current_run_idx = local_state.current_run_idx;

		try {
			batch = global_state.readers[current_run_idx]->ReadBatch(STANDARD_VECTOR_SIZE);
		} catch (const std::exception &e) {
			// Mid-stream read failure.
			//
			// We intentionally do NOT retry from scratch here. The previous
			// iterations of Execute already emitted rows from the first attempt
			// into downstream chunks (DuckDB consumes the output chunk before
			// the next Execute call). Reopening and re-reading would emit the
			// SAME underlying reads a second time, and with the counter reset
			// would give them duplicate sequence_index values too.
			//
			// Instead: skip the remainder of this run with a loud warning. The
			// user gets partial data up to the failure point, plus the run is
			// recorded in skipped_runs so the end-of-scan summary tells them
			// exactly which sample was truncated.
			auto &run = global_state.runs[current_run_idx];
			uint64_t emitted = global_state.run_sequence_counters[current_run_idx] > 0
			                       ? global_state.run_sequence_counters[current_run_idx] - 1
			                       : 0;
			Printer::PrintF("read_ena_sequences: WARNING: run '%s' failed mid-stream after emitting %llu read(s) "
			                "(%s); skipping remainder — downstream sees partial data for this run",
			                run.run_accession, static_cast<unsigned long long>(emitted), e.what());
			global_state.readers[current_run_idx].reset();
			{
				lock_guard<mutex> skip_guard(global_state.skipped_lock);
				global_state.skipped_runs.push_back(run.run_accession);
			}
			global_state.bytes_completed.fetch_add(run.total_bytes, std::memory_order_relaxed);
			global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
			local_state.has_run = false;
			continue;
		}

		if (batch.empty()) {
			// Run completed successfully — wait for ascp (if any) and release.
			// Finish() may throw on ascp non-zero exit; clean up the reader first
			// so the throw doesn't leave dangling state behind.
			auto &reader = *global_state.readers[current_run_idx];
			try {
				reader.Finish();
			} catch (...) {
				global_state.readers[current_run_idx].reset();
				throw;
			}
			global_state.readers[current_run_idx].reset();
			global_state.bytes_completed.fetch_add(global_state.runs[current_run_idx].total_bytes,
			                                       std::memory_order_relaxed);
			global_state.runs_completed.fetch_add(1, std::memory_order_relaxed);
			local_state.has_run = false;
			continue;
		}
		break;
	}

	const auto &run = global_state.runs[current_run_idx];
	auto &seq_counter = global_state.run_sequence_counters[current_run_idx];
	FillOutputFromBatch(output, batch, bind_data, run, seq_counter);
}

// ---- FillOutputFromBatch (shared by Execute + ExecuteInOut) ----
//
// Writes one DuckDB DataChunk from one SequenceRecordBatch. Caller owns the
// sequence counter (it lives on GlobalState for literal-path runs and on
// LocalState for in-out runs), so we take it by reference and post-increment
// one id per row.
void ReadENASequencesTableFunction::FillOutputFromBatch(DataChunk &output, const miint::SequenceRecordBatch &batch,
                                                        const Data &bind_data, const miint::ENARunInfo &run,
                                                        uint64_t &seq_counter) {
	idx_t count = batch.size();
	idx_t field_idx = 0;

	// sequence_index (column 0)
	auto seq_idx_data = FlatVector::GetData<int64_t>(output.data[field_idx++]);
	for (idx_t i = 0; i < count; i++) {
		seq_idx_data[i] = static_cast<int64_t>(seq_counter++);
	}

	// read_id (column 1)
	SetResultVectorString(output.data[field_idx++], batch.read_ids);

	// comment (column 2)
	SetResultVectorStringNullable(output.data[field_idx++], batch.comments);

	// sequence1 (column 3)
	SetResultVectorString(output.data[field_idx++], batch.sequences1);

	// sequence2 (column 4)
	if (batch.is_paired) {
		SetResultVectorString(output.data[field_idx++], batch.sequences2);
	} else {
		SetResultVectorNull(output.data[field_idx++]);
	}

	// qual1 (column 5)
	// SFF stores raw Phred internally with offset 33; user's qual_offset only applies to FASTQ
	uint8_t effective_qual_offset = (run.format == miint::ENASequenceFormat::SFF) ? 33 : bind_data.qual_offset;
	SetResultVectorListUInt8(output.data[field_idx++], batch.quals1, effective_qual_offset);

	// qual2 (column 6)
	if (batch.is_paired) {
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals2, effective_qual_offset);
	} else {
		SetResultVectorNull(output.data[field_idx++]);
	}

	// run_accession (column 7) - constant for all rows in this chunk
	SetResultVectorFilepath(output.data[field_idx++], run.run_accession);

	// sample_accession (column 8) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.sample_accession);

	// experiment_accession (column 9) - constant
	SetResultVectorFilepath(output.data[field_idx++], run.experiment_accession);

	// filepath (column 10) - optional
	if (bind_data.include_filepath) {
		std::string filepath;
		if (run.format == miint::ENASequenceFormat::SFF) {
			filepath = run.sff_url;
		} else if (bind_data.use_aspera && !run.aspera_paths.empty()) {
			filepath = run.aspera_paths[0].host + ":" + run.aspera_paths[0].remote_path;
			for (size_t u = 1; u < run.aspera_paths.size(); u++) {
				filepath += ";" + run.aspera_paths[u].host + ":" + run.aspera_paths[u].remote_path;
			}
		} else if (!run.fastq_urls.empty()) {
			filepath = run.fastq_urls[0];
			for (size_t u = 1; u < run.fastq_urls.size(); u++) {
				filepath += ";" + run.fastq_urls[u];
			}
		}
		SetResultVectorFilepath(output.data[field_idx++], filepath);
	}

	output.SetCardinality(count);
}

// ---- ExecuteInOut ----
//
// State machine for the lateral / in-out path. Mirrors the literal path's
// retry policy (Execute): one-shot retry on open failure, one-shot retry on
// mid-stream read failure, skip-with-warning on second failure. Returning
// `HAVE_MORE_OUTPUT` with cardinality=0 is safe because internal state
// (current_reader reset, pending_runs popped, lateral_sequence_counter reset)
// advances between callbacks — the dispatcher loops until state progresses to
// a data-emitting call or a NEED_MORE_INPUT / terminal return.
//
// CRITICAL user-facing invariant: every skip path prints a loud Printer
// warning via record_skip(). Silent skips in lateral mode are dangerous — the
// outer query just sees zero rows for the failed accession, indistinguishable
// from "no matches". The user depends on these warnings to notice lost
// samples. Both per-attempt retry messages AND final skip notices are emitted.
//
// `lateral_sequence_counter` is reset to 1 when a new outer row is consumed
// (Phase 1) AND on open-retry / mid-stream retry (matching the literal path's
// counter reset invariant #4 in localdocs/PLAN-ena-lateral.md). Within a
// single outer row, the counter increments monotonically across pending_runs
// so sequence_index is globally unique per outer row.
OperatorResultType ReadENASequencesTableFunction::ExecuteInOut(ExecutionContext &context, TableFunctionInput &data_p,
                                                               DataChunk &input, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &local = data_p.local_state->Cast<LocalState>();
	auto &global = data_p.global_state->Cast<GlobalState>();

	// Record a skipped accession and emit a loud warning. Both pieces are
	// user-critical: the user must know which samples they did not get data
	// for, so they can decide whether to retry manually or investigate.
	auto record_skip = [&bind_data](const std::string &outer_accession, const std::string &run_accession,
	                                const std::string &reason) {
		{
			std::lock_guard<std::mutex> guard(bind_data.lateral_skipped_lock);
			bind_data.lateral_skipped_accessions.push_back(run_accession.empty() ? outer_accession : run_accession);
		}
		if (run_accession.empty()) {
			Printer::PrintF("read_ena_sequences: WARNING: skipped outer accession '%s' — %s", outer_accession, reason);
		} else {
			Printer::PrintF("read_ena_sequences: WARNING: skipped run '%s' (from outer accession '%s') — %s",
			                run_accession, outer_accession, reason);
		}
	};

	// Construct a fresh PerRunReader for a given run. Used for initial open
	// and for retry-after-failure. Lateral path always uses HTTP (Bind rejects
	// download_method='aspera' in deferred mode).
	auto make_reader_for = [&](const miint::ENARunInfo &run) {
		return std::make_unique<miint::PerRunReader>(
		    global.fs, run, /*use_aspera=*/false, bind_data.trim, global.open_mutex,
		    static_cast<const miint::AsperaConfig *>(nullptr), bind_data.max_sequences);
	};

	// Phase 1: need a fresh outer row → pull and resolve.
	if (local.row_consumed && local.pending_runs.empty()) {
		if (input.size() == 0) {
			output.SetCardinality(0);
			return OperatorResultType::NEED_MORE_INPUT;
		}
		auto acc_val = input.data[0].GetValue(0);
		if (acc_val.IsNull()) {
			// NULL accession → no rows for this outer row. Not a failure
			// (outer emitted NULL intentionally); skip quietly.
			output.SetCardinality(0);
			return OperatorResultType::NEED_MORE_INPUT;
		}
		local.current_accession = acc_val.ToString();
		local.lateral_sequence_counter = 1;

		// Resolve against the bind-owned cache using the bind-owned ENAClient.
		// Sharing the client across outer rows means its ~3 req/sec rate
		// limiter actually throttles globally (not per-row), and parallel
		// LocalStates serialize on the client's internal mutex.
		miint::ENAFetcher fetcher = [&bind_data](const std::string &url) {
			return bind_data.lateral_client->FetchURL(url);
		};
		try {
			auto resolved = miint::ResolveRunsBatch(fetcher, *bind_data.resolver_cache, {local.current_accession},
			                                        bind_data.prefer_format);
			auto it = resolved.find(local.current_accession);
			if (it == resolved.end() || it->second.empty()) {
				record_skip(local.current_accession, "",
				            "ENA returned no runs for this accession (may not exist, may have no sequence data)");
				local.row_consumed = true;
				output.SetCardinality(0);
				return OperatorResultType::NEED_MORE_INPUT;
			}
			for (const auto &run : it->second) {
				local.pending_runs.push_back(run);
			}
		} catch (const std::exception &e) {
			record_skip(local.current_accession, "", std::string("metadata resolution failed: ") + e.what());
			local.row_consumed = true;
			output.SetCardinality(0);
			return OperatorResultType::NEED_MORE_INPUT;
		}
		local.row_consumed = false;
	}

	// Phase 2: ensure we have an open reader on the current pending run.
	// Retry-once-on-open mirrors Execute (literal path, invariant #3).
	if (!local.current_reader) {
		if (local.pending_runs.empty()) {
			local.row_consumed = true;
			output.SetCardinality(0);
			return OperatorResultType::NEED_MORE_INPUT;
		}
		auto run = local.pending_runs.front();
		local.pending_runs.pop_front();
		local.current_reader = make_reader_for(run);
		try {
			local.current_reader->Open();
		} catch (const std::exception &e) {
			Printer::PrintF("read_ena_sequences: warning: run '%s' failed to open (%s), retrying...", run.run_accession,
			                e.what());
			local.current_reader.reset();
			local.lateral_sequence_counter = 1; // No rows emitted yet; keep counter consistent.
			local.current_reader = make_reader_for(run);
			try {
				local.current_reader->Open();
			} catch (const std::exception &e2) {
				record_skip(local.current_accession, run.run_accession,
				            std::string("open failed on retry: ") + e2.what());
				local.current_reader.reset();
				// State advanced (run popped, reader reset); the dispatcher
				// will call us back to handle the next pending run or yield
				// NEED_MORE_INPUT. This is why HAVE_MORE_OUTPUT + cardinality=0
				// is safe here — forward progress is guaranteed.
				output.SetCardinality(0);
				return OperatorResultType::HAVE_MORE_OUTPUT;
			}
		}
	}

	// Phase 3: read a batch.
	//
	// IMPORTANT: we do NOT retry from scratch on a mid-stream failure. Previous
	// ExecuteInOut calls already streamed rows from this run downstream (the
	// output chunks are consumed between calls), and re-reading from the top
	// would emit the same reads again with the sequence_index counter reset.
	// Downstream would see duplicate (run_accession, sequence_index) pairs with
	// identical — or worse, subtly different — read_id content. That's a data
	// integrity hazard the user can't easily dedupe without guessing which
	// rows were doubled.
	//
	// Instead: skip the rest of this run and report LOUDLY exactly how much
	// was emitted before the failure. The user can decide to re-run the query
	// (which will hit the metadata cache) or accept the partial data.
	miint::SequenceRecordBatch batch;
	try {
		batch = local.current_reader->ReadBatch(STANDARD_VECTOR_SIZE);
	} catch (const std::exception &e) {
		const uint64_t emitted = local.lateral_sequence_counter > 0 ? local.lateral_sequence_counter - 1 : 0;
		const std::string failed_run = local.current_reader->Run().run_accession;
		record_skip(local.current_accession, failed_run,
		            "failed mid-stream after emitting " + std::to_string(emitted) +
		                " read(s); downstream sees partial data for this run (error: " + e.what() + ")");
		local.current_reader.reset();
		output.SetCardinality(0);
		return OperatorResultType::HAVE_MORE_OUTPUT;
	}

	if (batch.empty()) {
		// Reader exhausted — Finish() waits for ascp (no-op for HTTP) and may
		// throw on non-zero exit. In lateral mode we warn but do not abort:
		// the user already has whatever rows we emitted, and killing the
		// entire outer query because one run's completion code was non-zero
		// is overkill. The warning ensures the user still sees the problem.
		try {
			local.current_reader->Finish();
		} catch (const std::exception &e) {
			record_skip(local.current_accession, local.current_reader->Run().run_accession,
			            std::string("finish/transfer completion failed: ") + e.what());
		}
		local.current_reader.reset();
		if (!local.pending_runs.empty()) {
			output.SetCardinality(0);
			return OperatorResultType::HAVE_MORE_OUTPUT;
		}
		local.row_consumed = true;
		output.SetCardinality(0);
		return OperatorResultType::NEED_MORE_INPUT;
	}

	// Phase 4: emit the batch.
	FillOutputFromBatch(output, batch, bind_data, local.current_reader->Run(), local.lateral_sequence_counter);
	return OperatorResultType::HAVE_MORE_OUTPUT;
}

// ---- Progress ----

double ReadENASequencesTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                               const GlobalTableFunctionState *global_state) {
	if (bind_data) {
		auto &bd = bind_data->Cast<Data>();
		if (bd.deferred_resolution) {
			// Lateral mode: total work is unknown (driven by outer side).
			return -1.0;
		}
	}
	if (!global_state) {
		return -1.0;
	}
	auto &state = global_state->Cast<GlobalState>();
	if (state.total_runs == 0) {
		return 100.0;
	}

	// Prefer byte-based progress when ENA provided fastq_bytes metadata.
	// Falls back to run-count progress when byte sizes are unavailable.
	// Returning 0.0 when no work is done is safe — DuckDB's Kalman filter
	// defers initialization until progress > 0.
	if (state.total_bytes > 0) {
		uint64_t done = state.bytes_completed.load(std::memory_order_relaxed);
		return std::min(100.0, static_cast<double>(done) / static_cast<double>(state.total_bytes) * 100.0);
	}

	// Fallback: run-count based progress
	idx_t completed = state.runs_completed.load(std::memory_order_relaxed);
	return std::min(100.0, static_cast<double>(completed) / static_cast<double>(state.total_runs) * 100.0);
}

// ---- Registration ----

TableFunction ReadENASequencesTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_sequences", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	// Dual-path: literal / scalar args take `Execute` (parallel, byte-based
	// progress); correlated / subquery args take `ExecuteInOut` (one outer row
	// at a time, LIMIT-inside-LATERAL short-circuits via LocalState dtor).
	tf.in_out_function = ExecuteInOut;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	tf.named_parameters["download_method"] = LogicalType::VARCHAR;
	tf.named_parameters["prefer_format"] = LogicalType::VARCHAR;
	tf.named_parameters["trim_sff"] = LogicalType::BOOLEAN;
	tf.named_parameters["max_sequences"] = LogicalType::BIGINT;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	tf.table_scan_progress = Progress;
	return tf;
}

void ReadENASequencesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
