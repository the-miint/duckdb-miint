#pragma once

#include "line_file_scan.hpp"
#include "paf_parser.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_paf(path): one row per PAF record, coordinates in miint's convention
// (docs/reading.md, "Coordinate conventions").
class ReadPafTableFunction {
public:
	struct LocalState : public LineScanLocalState {
		miint::PafParser parser;
		miint::PafRecord record;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
