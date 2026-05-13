#pragma once

#include "ncbi_client.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>

namespace duckdb {

class ReadNCBITableFunction {
public:
	// Bind-time data: stores accessions and parameters
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions;           // Original, in input order.
		std::vector<std::vector<std::string>> batches; // Pre-chunked at Bind time.
		std::string api_key;
		int64_t batch_size;

		// Schema
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions, std::vector<std::vector<std::string>> batches,
		     const std::string &api_key, int64_t batch_size);
	};

	// Execution-wide state: manages fetched data
	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::unique_ptr<miint::NCBIClient> client;
		std::vector<miint::GenBankMetadata> metadata_results;
		size_t batch_cursor;  // Index of the next batch to fetch.
		size_t result_offset; // Read cursor inside metadata_results.

		// Missing-accession bookkeeping for the end-of-scan summary.
		std::vector<std::string> missing_accessions;
		std::atomic<bool> summary_emitted {false};

		// Note: batches is held by VALUE for the same reason as
		// ReadNCBIFastaTableFunction::GlobalState::work_units (bind/global
		// teardown order is not contractually fixed).
		std::vector<std::vector<std::string>> batches;

		GlobalState(DatabaseInstance &db, const std::string &api_key, std::vector<std::vector<std::string>> batches);

		idx_t MaxThreads() const override {
			return 1; // Single-threaded for rate limiting
		}

		// Fetch the next batch; parse + diff missing; emit per-batch warning.
		// Returns false when no more batches remain.
		bool FetchNextBatch(ClientContext &context);
	};

	// Per-thread state (minimal for single-threaded)
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
