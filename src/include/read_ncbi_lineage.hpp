#pragma once

#include "ncbi_client.hpp"
#include "taxonomy_lineage.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <atomic>
#include <memory>

namespace duckdb {

// read_ncbi_lineage(taxid, api_key := NULL, batch_size := 500)
// Resolves NCBI taxonomy IDs to their rank-collapsed lineage via E-utilities
// (db=taxonomy). Emits the shared lineage schema (identical to the taxonomy_lineage
// macro over a taxdump tree), so the online and offline paths are interchangeable.
class ReadNCBILineageTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> taxids;               // in input order
		std::vector<std::vector<std::string>> batches; // pre-chunked at Bind
		std::string api_key;
		int64_t batch_size;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> taxids, std::vector<std::vector<std::string>> batches, const std::string &api_key,
		     int64_t batch_size);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::unique_ptr<miint::NCBIClient> client;
		std::vector<miint::Lineage> results;
		size_t batch_cursor;
		size_t result_offset;

		std::vector<std::string> missing_taxids;
		std::atomic<bool> summary_emitted {false};

		std::vector<std::vector<std::string>> batches; // held by value (see read_ncbi)

		GlobalState(DatabaseInstance &db, const std::string &api_key, std::vector<std::vector<std::string>> batches);

		idx_t MaxThreads() const override {
			return 1; // single-threaded for rate limiting
		}

		bool FetchNextBatch(ClientContext &context);
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
