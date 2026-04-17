#pragma once
#include "SequenceReader.hpp"
#include "duckdb_seq_stream.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "SequenceRecord.hpp"
#include "QualScore.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <optional>
#include <thread>
#include <vector>

namespace duckdb {
class ReadFastxTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sequence1_paths;
		std::optional<std::vector<std::string>> sequence2_paths;
		bool include_filepath;
		bool uses_stdin;
		uint8_t qual_offset;
		uint64_t max_batch_bytes;

		std::vector<std::string> names; // field names
		std::vector<LogicalType> types; // field types

		Data(const std::vector<std::string> &r1_paths, const std::optional<std::vector<std::string>> &r2_paths,
		     bool include_fp, bool stdin_used, uint8_t offset, uint64_t max_bytes)
		    : sequence1_paths(r1_paths), sequence2_paths(r2_paths), include_filepath(include_fp),
		      uses_stdin(stdin_used), qual_offset(offset), max_batch_bytes(max_bytes),
		      names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"}),
		      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
		             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT),
		             LogicalType::LIST(LogicalType::UTINYINT)}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		};
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		// Readers are created lazily when a thread claims a file, not all upfront.
		// This avoids opening many HTTP connections simultaneously for remote file arrays.
		std::vector<std::unique_ptr<miint::SequenceReader>> readers;
		std::vector<std::string> sequence1_filepaths;
		std::vector<std::string> sequence2_filepaths;
		size_t next_file_idx; // Next file available for claiming
		bool uses_stdin;
		std::vector<uint64_t>
		    file_sequence_counters; // Per-file sequence counters (no atomic needed - file access is exclusive)
		FileSystem &fs;

		// stdin cannot be read in parallel (no seeking/rewinding).
		// This forces sequential execution, which may be slower than
		// reading from files where DuckDB can parallelize.
		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			return std::min<idx_t>(sequence1_filepaths.size(), std::thread::hardware_concurrency());
		};

		GlobalState(FileSystem &fs, const std::vector<std::string> &sequence1_paths,
		            const std::optional<std::vector<std::string>> &sequence2_paths, bool stdin_used)
		    : next_file_idx(0), uses_stdin(stdin_used), fs(fs) {
			sequence1_filepaths = sequence1_paths;
			if (sequence2_paths.has_value()) {
				sequence2_filepaths = sequence2_paths.value();
			}

			// Pre-allocate slots but don't open connections yet
			readers.resize(sequence1_paths.size());
			file_sequence_counters.resize(sequence1_paths.size(), 1);
		}

		// Open reader for a specific file index. Called under lock when a thread claims the file.
		void OpenReader(size_t file_idx) {
			if (readers[file_idx]) {
				return; // Already opened
			}
			bool r1_remote = miint::RemoteFileHelper::IsRemotePath(sequence1_filepaths[file_idx]);
			bool r2_remote =
			    !sequence2_filepaths.empty() && miint::RemoteFileHelper::IsRemotePath(sequence2_filepaths[file_idx]);

			if (r1_remote || r2_remote) {
				auto *s1 = CreateDuckDBSeqStream(fs, sequence1_filepaths[file_idx]);
				miint::DuckDBSeqStream *s2 = nullptr;
				if (!sequence2_filepaths.empty()) {
					s2 = CreateDuckDBSeqStream(fs, sequence2_filepaths[file_idx]);
				}
				readers[file_idx] =
				    std::make_unique<miint::SequenceReader>(s1, static_cast<miint::DuckDBSeqStream *>(s2));
			} else if (!sequence2_filepaths.empty()) {
				readers[file_idx] = std::make_unique<miint::SequenceReader>(sequence1_filepaths[file_idx],
				                                                            sequence2_filepaths[file_idx]);
			} else {
				readers[file_idx] = std::make_unique<miint::SequenceReader>(sequence1_filepaths[file_idx]);
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;

		LocalState() : current_file_idx(0), has_file(false) {
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
}; // namespace duckdb
