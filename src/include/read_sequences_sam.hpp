#pragma once
#include "SAMReader.hpp"
#include "QualScore.hpp"
#include "remote_file_helper.hpp"
#include "hfile_duckdb.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <thread>
#include <vector>

namespace duckdb {
class ReadSequencesSamTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;
		bool uses_stdin;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp, bool stdin_used)
		    : file_paths(paths), include_filepath(include_fp), uses_stdin(stdin_used),
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
		std::vector<std::unique_ptr<miint::SAMReader>> readers;
		std::vector<std::string> filepaths; // Original paths (for include_filepath)
		size_t next_file_idx;
		bool uses_stdin;
		std::vector<uint64_t> file_sequence_counters;

		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			auto hw_threads = std::thread::hardware_concurrency();
			if (hw_threads == 0) {
				hw_threads = 1;
			}
			return std::min<idx_t>(readers.size(), std::min<idx_t>(8, hw_threads));
		}

		GlobalState(const std::vector<std::string> &paths, FileSystem &fs, bool stdin_used)
		    : filepaths(paths), next_file_idx(0), uses_stdin(stdin_used) {
			for (const auto &path : paths) {
				try {
					if (miint::RemoteFileHelper::IsRemotePath(path)) {
						hFILE *hf = miint::hfile_duckdb_open(fs, path);
						if (!hf) {
							throw IOException("Failed to open remote file: " + path);
						}
						readers.push_back(std::make_unique<miint::SAMReader>(hf, path, /*include_seq_qual=*/true,
						                                                     /*require_references=*/false));
					} else {
						readers.push_back(std::make_unique<miint::SAMReader>(path, /*include_seq_qual=*/true,
						                                                     /*require_references=*/false));
					}
				} catch (std::exception &e) {
					throw IOException("Error opening '%s': %s", path, e.what());
				}
				file_sequence_counters.emplace_back(1);
			}

			// Multi-threaded BGZF decompression, but ONLY for a single input file.
			// With one file, MaxThreads() == 1 so the scan runs on one core while the
			// rest sit idle; an HTSlib worker pool decompresses blocks ahead of the parser
			// (blocks stay in order, so output is identical). With multiple files we rely
			// on file-level parallelism instead -- giving every reader its own pool would
			// oversubscribe (up to min(files, 8) readers run concurrently, each spawning a
			// pool), so we deliberately leave the multi-file path single-threaded per file.
			if (readers.size() == 1 && !uses_stdin) {
				auto hw = std::thread::hardware_concurrency();
				int decompress_threads = (hw > 1) ? std::min<int>(static_cast<int>(hw) - 1, 4) : 1;
				readers[0]->set_threads(decompress_threads);
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;
		// Reused across Execute calls to accumulate one chunk's raw quality bytes before a
		// single bulk copy into the qual1 LIST child (avoids per-record heap allocation).
		std::vector<uint8_t> qual_scratch;

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
} // namespace duckdb
