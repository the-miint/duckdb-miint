#pragma once

#include "SequenceReader.hpp"
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

class ReadENAFastqTableFunction {
public:
	struct RunInfo {
		std::string run_accession;
		std::string sample_accession;
		std::string experiment_accession;
		std::vector<std::string> fastq_urls; // 1 for single-end, 2 for paired-end
		bool is_paired;
	};

	struct Data : public TableFunctionData {
		std::vector<RunInfo> runs;
		bool include_filepath;
		uint8_t qual_offset;

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

		idx_t MaxThreads() const override {
			return std::min<idx_t>(runs.size(), 8);
		}

		GlobalState(FileSystem &fs, const std::vector<RunInfo> &runs);

		// Open reader for a specific run index.
		// Called outside the lock — each thread owns its run_idx exclusively.
		void OpenReader(size_t run_idx);
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
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
