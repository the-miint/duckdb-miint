#pragma once

#include "blast_client.hpp"
#include "blast_parser.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <mutex>

namespace duckdb {

class BlastSearchTableFunction {
public:
	static constexpr size_t DEFAULT_MAX_BATCH_BYTES = 4 * 1024 * 1024; // 4 MB

	struct Data : public TableFunctionData {
		std::string query_table;
		std::string program = "blastn";
		std::string database = "nt";
		double evalue = 10.0;
		int max_targets = 500;
		bool megablast = true;
		std::string api_key;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::unique_ptr<miint::BlastClient> client;
		std::vector<miint::BlastHit> result_buffer;
		idx_t result_offset = 0;

		std::vector<miint::BlastParser::SequenceBatch> batches;
		idx_t batch_cursor = 0;
		mutex lock;

		// Bind-time params copied for use in FetchNextBatch
		std::string program;
		std::string database;
		double evalue;
		int max_targets;
		bool megablast;

		idx_t MaxThreads() const override {
			return 1;
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
