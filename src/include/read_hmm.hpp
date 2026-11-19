#pragma once

#include "line_file_scan.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// read_hmm(path): one row per profile HMM in a HMMER text file (HMMER3/a to HMMER3/f, or
// HMMER2.0), its header and E-value statistics as typed columns.
class ReadHmmTableFunction {
public:
	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
