#include "SequenceReader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <chrono>
#include <cstdlib>
#include <optional>
#include <thread>
#include <unordered_map>
#include <unordered_set>
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

		std::vector<std::string> names;                 // field names
		std::vector<LogicalType> types;                 // field types
		std::vector<miint::SequenceRecordField> fields; // enum for convenience

		Data(const std::vector<std::string> &r1_paths, const std::optional<std::vector<std::string>> &r2_paths,
		     bool include_fp, bool stdin_used, uint8_t offset)
		    : sequence1_paths(r1_paths), sequence2_paths(r2_paths), include_filepath(include_fp),
		      uses_stdin(stdin_used), qual_offset(offset),
		      names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"}),
		      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
		             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT),
		             LogicalType::LIST(LogicalType::UTINYINT)}),
		      fields({miint::SequenceRecordField::READ_ID, miint::SequenceRecordField::COMMENT,
		              miint::SequenceRecordField::SEQUENCE1, miint::SequenceRecordField::SEQUENCE2,
		              miint::SequenceRecordField::QUAL1, miint::SequenceRecordField::QUAL2}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		};
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::vector<std::unique_ptr<miint::SequenceReader>> readers;
		std::vector<std::string> sequence1_filepaths;
		std::vector<std::string> sequence2_filepaths;
		size_t next_file_idx; // Next file available for claiming
		bool finished;
		bool uses_stdin;
		std::atomic<uint64_t> sequence_index_counter; // Atomic for thread-safe increments

		// Diagnostics
		std::atomic<size_t> total_execute_calls{0};
		std::atomic<size_t> total_rows_read{0};
		std::unordered_map<std::thread::id, std::vector<size_t>> thread_files_claimed;
		std::unordered_map<std::thread::id, size_t> thread_chunk_count;
		std::unordered_set<size_t> files_processed;
		std::chrono::steady_clock::time_point start_time;
		bool diagnostics_enabled;

		// stdin cannot be read in parallel (no seeking/rewinding).
		// This forces sequential execution, which may be slower than
		// reading from files where DuckDB can parallelize.
		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			// Allow up to min(files, 64) threads for per-file parallelism
			return std::max<idx_t>(8, std::min<idx_t>(readers.size(), 64));
		};

		GlobalState(const std::vector<std::string> &sequence1_paths,
		            const std::optional<std::vector<std::string>> &sequence2_paths, bool stdin_used)
		    : next_file_idx(0), finished(false), uses_stdin(stdin_used), sequence_index_counter(1),
		      start_time(std::chrono::steady_clock::now()) {
			sequence1_filepaths = sequence1_paths;
			if (sequence2_paths.has_value()) {
				sequence2_filepaths = sequence2_paths.value();
			}

			for (size_t i = 0; i < sequence1_paths.size(); i++) {
				if (sequence2_paths.has_value()) {
					readers.push_back(
					    std::make_unique<miint::SequenceReader>(sequence1_paths[i], sequence2_paths.value()[i]));
				} else {
					readers.push_back(std::make_unique<miint::SequenceReader>(sequence1_paths[i]));
				}
			}

			diagnostics_enabled = std::getenv("DUCKDB_MIINT_DIAGNOSTICS") != nullptr;
			if (diagnostics_enabled) {
				fprintf(stderr, "[READ_FASTX] Diagnostics enabled. Processing %zu files.\n",
				        sequence1_filepaths.size());
				fprintf(stderr, "[READ_FASTX] MaxThreads = %zu\n", MaxThreads());
			}
		};

		~GlobalState() {
			if (diagnostics_enabled && total_execute_calls.load() > 0) {
				auto duration = std::chrono::steady_clock::now() - start_time;
				auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();

				fprintf(stderr, "\n[READ_FASTX] === FINAL STATISTICS ===\n");
				fprintf(stderr, "[READ_FASTX] Total Execute() calls: %zu\n", total_execute_calls.load());
				fprintf(stderr, "[READ_FASTX] Total rows read: %zu\n", total_rows_read.load());
				fprintf(stderr, "[READ_FASTX] Duration: %ld seconds\n", seconds > 0 ? seconds : 1);
				fprintf(stderr, "[READ_FASTX] Unique threads: %zu\n", thread_files_claimed.size());
				fprintf(stderr, "[READ_FASTX] Files available: %zu\n", sequence1_filepaths.size());
				fprintf(stderr, "[READ_FASTX] Files actually processed: %zu\n", files_processed.size());

				if (seconds > 0) {
					fprintf(stderr, "[READ_FASTX] Throughput: %.2f rows/second\n",
					        total_rows_read.load() / (double)seconds);
				}

				fprintf(stderr, "[READ_FASTX] Thread details:\n");
				for (const auto &[tid, files] : thread_files_claimed) {
					size_t chunks = thread_chunk_count[tid];
					fprintf(stderr, "[READ_FASTX]   Thread %zu: claimed %zu files, processed %zu chunks\n",
					        std::hash<std::thread::id>{}(tid), files.size(), chunks);
					fprintf(stderr, "[READ_FASTX]     Files: ");
					for (size_t i = 0; i < files.size(); i++) {
						fprintf(stderr, "%zu%s", files[i], i < files.size() - 1 ? ", " : "");
					}
					fprintf(stderr, "\n");
				}
				fprintf(stderr, "[READ_FASTX] ========================\n\n");
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

	static void SetResultVector(Vector &result_vector, const miint::SequenceRecordField &field,
	                            const std::vector<miint::SequenceRecord> &records, uint8_t qual_offset);
	static void SetResultVectorNull(Vector &result_vector);
	static void SetResultVectorString(Vector &result_vector, const miint::SequenceRecordField &field,
	                                  const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorStringNullable(Vector &result_vector, const miint::SequenceRecordField &field,
	                                          const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorListUInt8(Vector &result_vector, const miint::SequenceRecordField &field,
	                                     const std::vector<miint::SequenceRecord> &records, uint8_t qual_offset);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath, size_t num_records);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
