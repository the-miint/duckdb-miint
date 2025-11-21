#include "SequenceReader.hpp"
#include "duckdb/common/exception.hpp"
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
		bool uses_stdin;
		std::atomic<uint64_t> sequence_index_counter; // Atomic for thread-safe increments

		// stdin cannot be read in parallel (no seeking/rewinding).
		// This forces sequential execution, which may be slower than
		// reading from files where DuckDB can parallelize.
		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			return std::min<idx_t>(readers.size(), std::thread::hardware_concurrency());
		};

		GlobalState(const std::vector<std::string> &sequence1_paths,
		            const std::optional<std::vector<std::string>> &sequence2_paths, bool stdin_used)
		    : next_file_idx(0), uses_stdin(stdin_used), sequence_index_counter(1) {
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

	static void SetResultVectorNull(Vector &result_vector);
	static void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values,
	                                     uint8_t qual_offset);
	static void SetResultVectorListUInt8Nullable(Vector &result_vector, const std::vector<miint::QualScore> &values,
	                                             uint8_t qual_offset, bool is_paired);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
