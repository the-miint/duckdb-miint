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

// read_ncbi_taxdump_names(source) -> (taxid, name, unique_name, name_class) from
// names.dmp, unfiltered. read_ncbi_taxdump emits only the scientific name, so this
// is how any other name class (notably "genbank common name") is reached.
class ReadNCBITaxdumpNamesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string source;   // empty => auto-download the default taxdump into the cache
		bool refresh = false; // re-download even if cached (only meaningful for the default source)
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<miint::TaxdumpName> names;
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

// read_ncbi_taxdump_deleted(source) -> (taxid) from delnodes.dmp.
// Deleted taxids are absent from the live tree and are not in merged.dmp either, so
// without this a retired taxid is indistinguishable from one that never existed.
class ReadNCBITaxdumpDeletedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string source;   // empty => auto-download the default taxdump into the cache
		bool refresh = false; // re-download even if cached (only meaningful for the default source)
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<int64_t> deleted;
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

// Load the requested taxdump member files from `source` (a directory of .dmp files or
// a .tar.gz archive), reading through the DuckDB FileSystem. Shared by all the table
// functions, each of which requests only the members it parses. Throws IOException on
// a missing/invalid source, or if a requested member is absent.
miint::TaxdumpFiles LoadTaxdumpFiles(ClientContext &context, const std::string &source,
                                     const miint::TaxdumpMemberSet &members);

} // namespace duckdb
