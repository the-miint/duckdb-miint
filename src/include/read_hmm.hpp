#pragma once

#include "line_file_scan.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <cstdint>
#include <string>

namespace duckdb {

// read_hmm(path): one row per profile HMM in a HMMER text file (HMMER3/a to HMMER3/f, or
// HMMER2.0), its header and E-value statistics as typed columns.
class ReadHmmTableFunction {
public:
	struct LocalState : public LineScanLocalState {
		// The current file's model being read: its lines so far, each ending in '\n', from the
		// line after the previous model's "//" (or the file's first line).
		std::string record;
		// The file line `record` starts at.
		uint64_t record_start = 1;
		// Whether `record` holds a line that is neither blank nor a comment.
		bool record_has_text = false;
		// The models of the current file read so far.
		uint64_t models_read = 0;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
