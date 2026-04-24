#pragma once

#include "ena_client.hpp"
#include "ena_parser.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <mutex>
#include <vector>

namespace duckdb {

class ReadENAAttributesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions;
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		unique_ptr<miint::ENAClient> client;

		// Streaming state. `ResolveAccessions` populates `pending_sample_accs`
		// from the input accessions. Each call to `FetchNextBatch` pops up to
		// BATCH_SIZE samples from `pending_sample_accs`, fetches + parses their
		// XML, and replaces `current_batch`. Execute() emits from
		// `current_batch[current_batch_offset..]` one DuckDB chunk at a time.
		//
		// BATCH_SIZE is the main lever for throughput. ENA accepts comma-
		// separated accession lists up to a fairly generous URL length; 200
		// sample accessions (~14 chars each) stays well under 8 KB. Compared
		// to 50, this cuts HTTP round-trips 4x — relevant for large studies
		// like PRJEB11419 with 33k samples where the resolver rate limit
		// (~3 req/s) is otherwise the wall-clock bottleneck.
		static constexpr size_t BATCH_SIZE = 200;
		bool resolved = false;
		std::vector<std::string> pending_sample_accs;
		std::vector<miint::SampleAttribute> current_batch;
		size_t current_batch_offset = 0;

		// Progress tracking. Byte-based isn't meaningful here (XML size per
		// sample is unknown up front) so we use samples-fetched / samples-total.
		// Both members are atomic so the progress callback (called from
		// DuckDB's query monitor thread, with a const GlobalState&) can read
		// without contending on `lock`. Published via table_scan_progress so
		// users see e.g. "[42%]" during long scans instead of guessing whether
		// the query is hung.
		std::atomic<size_t> total_samples_expected {0};
		std::atomic<size_t> samples_fetched {0};

		GlobalState(DatabaseInstance &db);

		// One-time sample-accession resolution. Populates `pending_sample_accs`.
		void ResolveAccessions(const std::vector<std::string> &accessions);

		// Pop up to BATCH_SIZE sample accessions, fetch + parse their XML, and
		// refill `current_batch`. Sets `current_batch_offset` to 0. Returns
		// true if a batch was fetched, false if no pending samples remain.
		bool FetchNextBatch();

		idx_t MaxThreads() const override {
			return 1;
		}

	private:
		std::vector<std::string> ResolveSampleAccessions(const std::vector<std::string> &accessions);
	};

	struct LocalState : public LocalTableFunctionState {};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	static double Progress(ClientContext &context, const FunctionData *bind_data,
	                       const GlobalTableFunctionState *global_state);
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
