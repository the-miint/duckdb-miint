#pragma once

#include "SFFReader.hpp"
#include "SequenceReader.hpp"
#include "aspera_stream.hpp"
#include "aspera_utils.hpp"
#include "duckdb_seq_stream.hpp"
#include "ena_client.hpp"
#include "ena_parser.hpp"
#include "ena_resolver_cache.hpp"
#include "ena_run_info_extractor.hpp"
#include "per_run_reader.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <deque>
#include <memory>
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

class ReadENASequencesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<miint::ENARunInfo> runs;
		bool include_filepath;
		uint8_t qual_offset;
		bool use_aspera;
		bool trim;
		std::string prefer_format;
#if MIINT_ASPERA_SUPPORTED
		miint::AsperaConfig aspera_config;
#endif

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		// Lateral / in-out support: when true, `runs` is empty at Bind time and
		// ExecuteInOut resolves each outer-row's accession via `resolver_cache`.
		// The ENAClient is shared (not per-LocalState) so its ~3 req/sec rate
		// limit actually throttles across outer rows. Its internal mutex
		// serializes concurrent access from parallel LocalStates.
		bool deferred_resolution = false;
		std::unique_ptr<miint::ENAResolverCache> resolver_cache;
		std::unique_ptr<miint::ENAClient> lateral_client;
		DatabaseInstance *db_ptr = nullptr;

		// Skipped-accession tracking for the lateral path.
		// ExecuteInOut has no natural "end of query" callback, so we surface
		// skip warnings per-event rather than as a batched summary. The list
		// is kept for potential future summary use. Both members are mutable
		// so they can be written from ExecuteInOut, which sees bind_data as
		// const (DuckDB shares it across threads during execution).
		mutable std::mutex lateral_skipped_lock;
		mutable std::vector<std::string> lateral_skipped_accessions;

		Data(std::vector<miint::ENARunInfo> runs, bool include_fp, uint8_t offset, bool trim,
		     const std::string &prefer_format);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		mutex open_mutex; // Serializes file opens to avoid overwhelming remote servers
		std::vector<std::unique_ptr<miint::PerRunReader>> readers;
		std::vector<miint::ENARunInfo> runs;
		size_t next_run_idx;
		std::vector<uint64_t> run_sequence_counters;
		FileSystem &fs;
		bool use_aspera;
		bool trim;

		// Progress tracking (byte-based when available, run-count fallback)
		std::atomic<idx_t> runs_completed;
		idx_t total_runs;
		std::atomic<uint64_t> bytes_completed;
		uint64_t total_bytes; // Sum across all runs; 0 means fall back to run-count progress

		// Skipped runs (transient failures after retry)
		mutex skipped_lock;
		std::vector<std::string> skipped_runs;
		std::atomic<bool> skipped_warned;

#if MIINT_ASPERA_SUPPORTED
		miint::AsperaConfig aspera_config;
#endif

		// Literal path: cap at 8 worker threads. In deferred (in-out) mode, runs
		// is empty at Bind time so this returns 0 — harmless because DuckDB's
		// PhysicalTableInOutFunction drives parallelism from the outer side, not
		// from this hint. (See duckdb/src/execution/operator/projection/
		// physical_tableinout_function.cpp — MaxThreads on GlobalOperatorState is
		// not consulted for in-out operators.)
		idx_t MaxThreads() const override {
			return std::min<idx_t>(runs.size(), 8);
		}

		GlobalState(FileSystem &fs, const std::vector<miint::ENARunInfo> &runs, bool use_aspera, bool trim);
	};

	struct LocalState : public LocalTableFunctionState {
		// Standard path (literal args)
		size_t current_run_idx;
		bool has_run;

		// Lateral / in-out path. current_reader owns the active run; when it
		// finishes we pop from pending_runs (a single outer accession may resolve
		// to multiple runs, e.g. SFF studies that flatten one-RunInfo-per-file).
		// row_consumed=true means we need a new outer row on the next call.
		std::unique_ptr<miint::PerRunReader> current_reader;
		std::deque<miint::ENARunInfo> pending_runs;
		std::string current_accession;
		uint64_t lateral_sequence_counter = 1;
		bool row_consumed = true;

		LocalState() : current_run_idx(0), has_run(false) {
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	// In-out (lateral / correlated-argument) execute. Processes one outer row
	// at a time, resolving the accession lazily against Data::resolver_cache.
	static OperatorResultType ExecuteInOut(ExecutionContext &context, TableFunctionInput &data_p, DataChunk &input,
	                                       DataChunk &output);
	// Shared column-fill helper used by both Execute (literal path) and
	// ExecuteInOut (lateral path). Takes the sequence counter by reference so
	// the caller (which owns the counter) can track progress across batches.
	static void FillOutputFromBatch(DataChunk &output, const miint::SequenceRecordBatch &batch, const Data &bind_data,
	                                const miint::ENARunInfo &run, uint64_t &seq_counter);
	static double Progress(ClientContext &context, const FunctionData *bind_data,
	                       const GlobalTableFunctionState *global_state);
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
