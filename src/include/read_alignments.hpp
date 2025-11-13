#pragma once
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
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
class ReadAlignmentsTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sam_paths;
		std::optional<std::string> reference_lengths_table;
		bool include_filepath;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
		std::vector<miint::SAMRecordField> fields;

		explicit Data(const std::vector<std::string> &paths, const std::optional<std::string> &ref_table,
		              bool include_fp)
		    : sam_paths(paths), reference_lengths_table(ref_table), include_filepath(include_fp),
		      names({"read_id", "flags",          "reference",     "position",        "stop_position", "mapq",
		             "cigar",   "mate_reference", "mate_position", "template_length", "tag_as",        "tag_xs",
		             "tag_ys",  "tag_xn",         "tag_xm",        "tag_xo",          "tag_xg",        "tag_nm",
		             "tag_yt",  "tag_md",         "tag_sa"}),
		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::BIGINT,    // position
		             LogicalType::BIGINT,    // stop_position
		             LogicalType::UTINYINT,  // mapq
		             LogicalType::VARCHAR,   // cigar
		             LogicalType::VARCHAR,   // mate_reference
		             LogicalType::BIGINT,    // mate_position
		             LogicalType::BIGINT,    // template_length
		             LogicalType::BIGINT,    // tag_as
		             LogicalType::BIGINT,    // tag_xs
		             LogicalType::BIGINT,    // tag_ys
		             LogicalType::BIGINT,    // tag_xn
		             LogicalType::BIGINT,    // tag_xm
		             LogicalType::BIGINT,    // tag_xo
		             LogicalType::BIGINT,    // tag_xg
		             LogicalType::BIGINT,    // tag_nm
		             LogicalType::VARCHAR,   // tag_yt
		             LogicalType::VARCHAR,   // tag_md
		             LogicalType::VARCHAR}), // tag_sa
		      fields({miint::SAMRecordField::READ_ID,       miint::SAMRecordField::FLAGS,
		              miint::SAMRecordField::REFERENCE,     miint::SAMRecordField::POSITION,
		              miint::SAMRecordField::STOP_POSITION, miint::SAMRecordField::MAPQ,
		              miint::SAMRecordField::CIGAR,         miint::SAMRecordField::MATE_REFERENCE,
		              miint::SAMRecordField::MATE_POSITION, miint::SAMRecordField::TEMPLATE_LENGTH,
		              miint::SAMRecordField::TAG_AS,        miint::SAMRecordField::TAG_XS,
		              miint::SAMRecordField::TAG_YS,        miint::SAMRecordField::TAG_XN,
		              miint::SAMRecordField::TAG_XM,        miint::SAMRecordField::TAG_XO,
		              miint::SAMRecordField::TAG_XG,        miint::SAMRecordField::TAG_NM,
		              miint::SAMRecordField::TAG_YT,        miint::SAMRecordField::TAG_MD,
		              miint::SAMRecordField::TAG_SA}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::vector<std::unique_ptr<miint::SAMReader>> readers;
		std::vector<std::string> filepaths;
		size_t next_file_idx; // Next file available for claiming
		bool finished;

		// Diagnostics
		std::atomic<size_t> total_execute_calls{0};
		std::atomic<size_t> total_rows_read{0};
		std::unordered_map<std::thread::id, std::vector<size_t>> thread_files_claimed; // Track ALL files per thread
		std::unordered_map<std::thread::id, size_t> thread_chunk_count;
		std::unordered_set<size_t> files_processed; // Track which files were actually read
		std::chrono::steady_clock::time_point start_time;
		bool diagnostics_enabled;

		idx_t MaxThreads() const override {
			// Allow up to min(files, 64) threads for per-file parallelism
			// But always allow at least 8 threads even with few files
			return std::max<idx_t>(8, std::min<idx_t>(readers.size(), 64));
		}

		GlobalState(const std::vector<std::string> &paths,
		            std::optional<std::unordered_map<std::string, uint64_t>> ref_lengths)
		    : next_file_idx(0), finished(false), start_time(std::chrono::steady_clock::now()) {
			filepaths = paths;
			for (const auto &path : paths) {
				if (ref_lengths.has_value()) {
					readers.push_back(std::make_unique<miint::SAMReader>(path, ref_lengths.value()));
				} else {
					readers.push_back(std::make_unique<miint::SAMReader>(path));
				}
			}
			diagnostics_enabled = std::getenv("DUCKDB_MIINT_DIAGNOSTICS") != nullptr;
			if (diagnostics_enabled) {
				fprintf(stderr, "[READ_ALIGNMENTS] Diagnostics enabled. Processing %zu files.\n", filepaths.size());
				fprintf(stderr, "[READ_ALIGNMENTS] MaxThreads = %zu\n", MaxThreads());
			}
		}

		~GlobalState() {
			if (diagnostics_enabled && total_execute_calls.load() > 0) {
				auto duration = std::chrono::steady_clock::now() - start_time;
				auto seconds = std::chrono::duration_cast<std::chrono::seconds>(duration).count();

				fprintf(stderr, "\n[READ_ALIGNMENTS] === FINAL STATISTICS ===\n");
				fprintf(stderr, "[READ_ALIGNMENTS] Total Execute() calls: %zu\n", total_execute_calls.load());
				fprintf(stderr, "[READ_ALIGNMENTS] Total rows read: %zu\n", total_rows_read.load());
				fprintf(stderr, "[READ_ALIGNMENTS] Duration: %ld seconds\n", seconds > 0 ? seconds : 1);
				fprintf(stderr, "[READ_ALIGNMENTS] Unique threads: %zu\n", thread_files_claimed.size());
				fprintf(stderr, "[READ_ALIGNMENTS] Files available: %zu\n", filepaths.size());
				fprintf(stderr, "[READ_ALIGNMENTS] Files actually processed: %zu\n", files_processed.size());

				if (seconds > 0) {
					fprintf(stderr, "[READ_ALIGNMENTS] Throughput: %.2f rows/second\n",
					        total_rows_read.load() / (double)seconds);
				}

				fprintf(stderr, "[READ_ALIGNMENTS] Thread details:\n");
				for (const auto &[tid, files] : thread_files_claimed) {
					size_t chunks = thread_chunk_count[tid];
					fprintf(stderr, "[READ_ALIGNMENTS]   Thread %zu: claimed %zu files, processed %zu chunks\n",
					        std::hash<std::thread::id>{}(tid), files.size(), chunks);
					fprintf(stderr, "[READ_ALIGNMENTS]     Files: ");
					for (size_t i = 0; i < files.size(); i++) {
						fprintf(stderr, "%zu%s", files[i], i < files.size() - 1 ? ", " : "");
					}
					fprintf(stderr, "\n");
				}
				fprintf(stderr, "[READ_ALIGNMENTS] ========================\n\n");
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

	static void SetResultVector(Vector &result_vector, const miint::SAMRecordField &field,
	                            const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorString(Vector &result_vector, const miint::SAMRecordField &field,
	                                  const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorStringNullable(Vector &result_vector, const miint::SAMRecordField &field,
	                                          const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorUInt8(Vector &result_vector, const miint::SAMRecordField &field,
	                                 const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorUInt16(Vector &result_vector, const miint::SAMRecordField &field,
	                                  const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorInt64(Vector &result_vector, const miint::SAMRecordField &field,
	                                 const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorInt64Nullable(Vector &result_vector, const miint::SAMRecordField &field,
	                                         const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath, size_t num_records);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
