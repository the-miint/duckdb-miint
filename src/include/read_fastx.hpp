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
		std::string read1_path;
		std::optional<std::string> read2_path;

		std::vector<std::string> names;                 // field names
		std::vector<LogicalType> types;                 // field types
		std::vector<miint::SequenceRecordField> fields; // enum for convenience

		Data(const std::string &r1_path, const std::string &r2_path)
		    : names({"read_id", "comment", "read1", "read2", "qual1", "qual2"}),
		      types({LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
		             LogicalType::LIST(LogicalType::UTINYINT), LogicalType::LIST(LogicalType::UTINYINT)}),
		      fields({miint::SequenceRecordField::READ_ID, miint::SequenceRecordField::COMMENT,
		              miint::SequenceRecordField::READ1, miint::SequenceRecordField::READ2,
		              miint::SequenceRecordField::QUAL1, miint::SequenceRecordField::QUAL2}) {
			read1_path = r1_path;
			if (!r2_path.empty()) {
				read2_path = r2_path;
			} else {
				read2_path = std::nullopt;
			}
		};
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		miint::SequenceReader reader;

		// A local use at 4 threads exhibited reduced performance, and is likely IO bound on read.
		// TODO: read R1/R2 in parallel within an execution block rather than each execution block
		// being single threaded.
		idx_t MaxThreads() const override {
			return 4;
		};

		GlobalState(const std::string &read1_path, const std::optional<std::string> &read2_path)
		    : reader(miint::SequenceReader(read1_path, read2_path)) {};
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

	static TableFunction GetFunction();
};
}; // namespace duckdb
