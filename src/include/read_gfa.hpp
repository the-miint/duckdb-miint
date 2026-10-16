#pragma once

#include "gfa_parser.hpp"
#include "line_file_scan.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_gfa(path): one row per GFA record, split at its record type's mandatory field count
// (docs/reading.md, "GFA").
class ReadGfaTableFunction {
public:
	struct LocalState : public LineScanLocalState {
		miint::GfaParser parser;
		miint::GfaRecord record;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
