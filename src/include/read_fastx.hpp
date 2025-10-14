#include "SequenceReader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "SequenceRecord.hpp"
#include <optional>

namespace duckdb {
class ReadFastxTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sequence1_paths;
		std::optional<std::vector<std::string>> sequence2_paths;
		bool include_filepath;
		bool uses_stdin;

		std::vector<std::string> names;                 // field names
		std::vector<LogicalType> types;                 // field types
		std::vector<miint::SequenceRecordField> fields; // enum for convenience

		Data(const std::vector<std::string> &r1_paths, const std::optional<std::vector<std::string>> &r2_paths,
		     bool include_fp, bool stdin_used)
		    : sequence1_paths(r1_paths), sequence2_paths(r2_paths), include_filepath(include_fp),
		      uses_stdin(stdin_used),
		      names({"read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"}),
		      types({LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
		             LogicalType::LIST(LogicalType::UTINYINT), LogicalType::LIST(LogicalType::UTINYINT)}),
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
		size_t current_file_idx;
		bool finished;
		std::string current_filepath;
		bool uses_stdin;

		// stdin cannot be read in parallel (no seeking/rewinding).
		// This forces sequential execution, which may be slower than
		// reading from files where DuckDB can parallelize.
		// Internal threading was removed from read_pe() in favor of DuckDB's
		// parallel execution which provides better resource management.
		idx_t MaxThreads() const override {
			return uses_stdin ? 1 : 8;
		};

		GlobalState(const std::vector<std::string> &sequence1_paths,
		            const std::optional<std::vector<std::string>> &sequence2_paths, bool stdin_used)
		    : current_file_idx(0), finished(false), uses_stdin(stdin_used) {
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

			if (!readers.empty()) {
				current_filepath = sequence1_paths[0];
			}
		};
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void SetResultVector(Vector &result_vector, const miint::SequenceRecordField &field,
	                            const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorNull(Vector &result_vector);
	static void SetResultVectorString(Vector &result_vector, const miint::SequenceRecordField &field,
	                                  const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorStringNullable(Vector &result_vector, const miint::SequenceRecordField &field,
	                                          const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorListUInt8(Vector &result_vector, const miint::SequenceRecordField &field,
	                                     const std::vector<miint::SequenceRecord> &records);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath, size_t num_records);

	static TableFunction GetFunction();
};
}; // namespace duckdb
