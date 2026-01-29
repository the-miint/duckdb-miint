#pragma once

#include "ncbi_client.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ReadNCBIFastaTableFunction {
public:
	// Bind-time data: stores accessions and parameters
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions;
		std::string api_key;
		bool include_filepath;

		// Schema (matches read_fastx)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions, const std::string &api_key, bool include_filepath);
	};

	// Execution-wide state: manages fetched data
	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::unique_ptr<miint::NCBIClient> client;
		miint::SequenceRecordBatch current_batch;
		size_t next_accession_idx;
		size_t batch_offset;           // Current position within batch
		int64_t sequence_index;        // Global sequence counter
		std::string current_accession; // Track current accession for filepath

		GlobalState(DatabaseInstance &db, const std::string &api_key, const std::vector<std::string> &accessions);

		idx_t MaxThreads() const override {
			return 1; // Single-threaded for rate limiting
		}

		// Fetch next accession's FASTA data
		bool FetchNextAccession();

	private:
		const std::vector<std::string> &accessions;
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
