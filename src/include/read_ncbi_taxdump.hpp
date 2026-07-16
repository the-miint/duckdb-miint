#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "taxdump_archive.hpp"
#include "taxdump_parser.hpp"

#include <string>
#include <vector>

namespace duckdb {

// read_ncbi_taxdump(source) -> taxonomy tree table
//   node_index (=tax_id), parent_index (NULL at root), name (scientific), rank, is_tip
// `source` is a directory of extracted .dmp files or a taxdump .tar.gz archive.
class ReadNCBITaxdumpTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string source;   // empty => auto-download the default taxdump into the cache
		bool refresh = false; // re-download even if cached (only meaningful for the default source)
	};

	// The whole tree is parsed once in InitGlobal and streamed from here.
	struct GlobalState : public GlobalTableFunctionState {
		std::vector<miint::TaxdumpNode> nodes;
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
	static void Register(ExtensionLoader &loader);
};

// read_ncbi_taxdump_merged(source) -> (old_taxid, new_taxid) from merged.dmp.
// The tree returned by read_ncbi_taxdump holds live nodes only; this exposes the
// retired -> current remap so old taxids can be updated before joining the tree.
class ReadNCBITaxdumpMergedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string source;   // empty => auto-download the default taxdump into the cache
		bool refresh = false; // re-download even if cached (only meaningful for the default source)
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<miint::TaxdumpMerge> merged;
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
	static void Register(ExtensionLoader &loader);
};

// Load the taxdump member files from `source` (a directory of .dmp files or a
// .tar.gz archive), reading through the DuckDB FileSystem. Shared by both table
// functions. Throws IOException on a missing/invalid source.
miint::TaxdumpFiles LoadTaxdumpFiles(ClientContext &context, const std::string &source);

} // namespace duckdb
