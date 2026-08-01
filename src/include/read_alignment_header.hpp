#pragma once

#include "duckdb.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>
#include <vector>

namespace duckdb {

// read_alignment_header(path) -- expose a SAM/BAM header's @SQ lines as
// (tid, reference, length).
//
// A record's tid is its index in the @SQ list, and a BAM is coordinate-sorted by
// TID rather than by reference name. Nothing else in miint exposes it:
// read_alignments/read_sam decode each record's `reference` via sam_hdr_tid2name
// against the file's own header, and name<->tid is a bijection fixed when the
// header is built, so reading names back yields ascending NAMES whatever the @SQ
// order was -- on correct and buggy code alike. That made the @SQ ordering
// contract in copy_sam unverifiable for BAM from SQL (issue #174).
//
// Reads only the header: sam_read1 is never called, so a file with intact @SQ
// lines and a corrupt record is read successfully here even though read_sam
// fails the scan on it.
//
// Local paths only. Remote paths are rejected at bind rather than half-supported;
// read_alignments' hfile_duckdb_open streaming path is deliberately not
// replicated (see docs/reading.md).
class ReadAlignmentHeaderTableFunction {
public:
	struct SQEntry {
		int32_t tid;
		std::string reference;
		int64_t length;
	};

	struct Data : public TableFunctionData {
		std::string path;
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<SQEntry> entries;
		idx_t cursor = 0;
		idx_t MaxThreads() const override {
			return 1;
		}
	};

	struct LocalState : public LocalTableFunctionState {};

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
