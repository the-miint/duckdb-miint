#pragma once

#include "SequenceReader.hpp"
#include "aspera_stream.hpp"
#include "aspera_utils.hpp"
#include "duckdb_seq_stream.hpp"
#include "ena_client.hpp"
#include "ena_parser.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

class ReadENAFastxTableFunction {
public:
	struct RunInfo {
		std::string run_accession;
		std::string sample_accession;
		std::string experiment_accession;
		std::vector<std::string> fastq_urls;         // HTTPS URLs (1 for single-end, 2 for paired-end)
		std::vector<miint::AsperaPath> aspera_paths; // Parsed host + remote_path pairs
		bool is_paired;
		uint64_t total_bytes = 0; // Sum of fastq_bytes for files in this run (0 if unavailable)
	};

	struct Data : public TableFunctionData {
		std::vector<RunInfo> runs;
		bool include_filepath;
		uint8_t qual_offset;
		bool use_aspera;
#if MIINT_ASPERA_SUPPORTED
		miint::AsperaConfig aspera_config;
#endif

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<RunInfo> runs, bool include_fp, uint8_t offset);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		mutex open_mutex; // Serializes file opens to avoid overwhelming remote servers
		std::vector<std::unique_ptr<miint::SequenceReader>> readers;
		std::vector<RunInfo> runs;
		size_t next_run_idx;
		std::vector<uint64_t> run_sequence_counters;
		FileSystem &fs;
		bool use_aspera;

		// Progress tracking (byte-based when fastq_bytes available, run-count fallback)
		std::atomic<idx_t> runs_completed;
		idx_t total_runs;
		std::atomic<uint64_t> bytes_completed;
		uint64_t total_bytes; // Sum across all runs; 0 means fall back to run-count progress

		// Skipped runs (transient failures after retry)
		mutex skipped_lock;
		std::vector<std::string> skipped_runs;
		std::atomic<bool> skipped_warned;

#if MIINT_ASPERA_SUPPORTED
		std::vector<std::unique_ptr<miint::AsperaProcess>> aspera_processes;
		std::vector<std::string> temp_file_paths;
		miint::AsperaConfig aspera_config;
#endif

		idx_t MaxThreads() const override {
			return std::min<idx_t>(runs.size(), 8);
		}

		GlobalState(FileSystem &fs, const std::vector<RunInfo> &runs, bool use_aspera);
		~GlobalState();

		// Open reader for a specific run index (HTTP path).
		void OpenReader(size_t run_idx);

#if MIINT_ASPERA_SUPPORTED
		// Open reader for a specific run index (Aspera path).
		// Reads from the shared aspera_process pipe.
		void OpenReaderAspera(size_t run_idx);
#endif
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_run_idx;
		bool has_run;

		LocalState() : current_run_idx(0), has_run(false) {
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	static double Progress(ClientContext &context, const FunctionData *bind_data,
	                       const GlobalTableFunctionState *global_state);
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
