#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

class WoltkaOguPerSampleBatchedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string relation_name;
		std::string sample_id_field;
		std::string sequence_id_field;
		idx_t batch_size;

		std::vector<std::string> names = {"sample_id", "feature_id", "value"};
		std::vector<LogicalType> types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE};

		Data(const std::string &rel_name, const std::string &sample_field, const std::string &seq_field,
		     idx_t batch_sz)
		    : relation_name(rel_name), sample_id_field(sample_field), sequence_id_field(seq_field),
		      batch_size(batch_sz) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::vector<std::string> all_samples;
		idx_t current_batch_idx;

		// Connection and streaming state for current batch
		unique_ptr<Connection> connection;
		unique_ptr<QueryResult> current_batch_result;

		GlobalState() : current_batch_idx(0) {
		}

		idx_t MaxThreads() const override {
			// Single-threaded execution: process batches sequentially to control memory
			// The woltka macro query within each batch can still use internal parallelism
			return 1;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);

private:
	// Helper to process a fetched chunk into output
	static void ProcessChunkToOutput(const unique_ptr<DataChunk> &chunk, DataChunk &output);

	// Helper to start processing the next batch
	// Returns true if batch has data and first chunk is in output, false if batch was empty
	static bool TryStartNextBatch(GlobalState &gstate, const Data &data, DataChunk &output);
};

} // namespace duckdb
