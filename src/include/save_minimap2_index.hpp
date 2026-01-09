#pragma once
#include "Minimap2Aligner.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <vector>

namespace duckdb {

class SaveMinimap2IndexTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string subject_table;
		std::string output_path;
		miint::Minimap2Config config;
		std::vector<miint::AlignmentSubject> subjects;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names({"success", "index_path", "num_subjects"}),
		      types({LogicalType::BOOLEAN, LogicalType::VARCHAR, LogicalType::BIGINT}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		bool done;

		idx_t MaxThreads() const override {
			return 1;
		}

		GlobalState() : done(false) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		LocalState() {
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
