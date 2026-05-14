#pragma once

#include "ncbi_client.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>

namespace duckdb {

class ReadNCBIFastaTableFunction {
public:
	// A unit of fetch work: either one batched epost+efetch of N sequence accessions,
	// or a single Datasets-API assembly download (which can't batch). Pre-built at
	// Bind-time so input order is preserved across mixed inputs — consecutive
	// sequence accessions group into a single batch, an assembly flushes the batch
	// and emits a singleton unit, then the next run of sequences resumes batching.
	enum class WorkUnitKind { SEQUENCE_BATCH, ASSEMBLY };
	struct WorkUnit {
		WorkUnitKind kind;
		std::vector<std::string> accessions; // 1 for ASSEMBLY, up to batch_size for SEQUENCE_BATCH
	};

	// Bind-time data: stores accessions and parameters
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions; // Original, in input order — kept for diagnostics.
		std::vector<WorkUnit> work_units;    // Pre-partitioned fetch plan.
		std::string api_key;
		bool include_filepath;
		int64_t batch_size;

		// Schema (matches read_fastx)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions, std::vector<WorkUnit> work_units, const std::string &api_key,
		     bool include_filepath, int64_t batch_size);
	};

	// Execution-wide state: manages fetched data
	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::unique_ptr<miint::NCBIClient> client;
		miint::SequenceRecordBatch current_batch;
		size_t work_cursor;     // Index of the next WorkUnit to fetch.
		size_t batch_offset;    // Current position within current_batch.
		int64_t sequence_index; // Global sequence counter.

		// Pre-built URL for the filepath column. For batched sequence fetches it
		// reflects the eutils efetch endpoint with the comma-joined ID list; for
		// assemblies it reflects the Datasets API download URL. Empty until first fetch.
		std::string current_filepath_url;

		// Missing-accession bookkeeping. NCBI silently omits invalid IDs from
		// batched responses; we accumulate them here so end-of-scan emits one
		// summary row in miint_warnings() listing every loss across the query.
		std::vector<std::string> missing_accessions;
		std::atomic<bool> summary_emitted {false}; // CAS-guarded once-only emission.

		// Note: work_units is held by VALUE, not a reference into Data. DuckDB's
		// teardown order between bind data and global state is not contractually
		// guaranteed; copying the small vector at scan-init removes any risk
		// of a dangling reference during long-running fetches.
		std::vector<WorkUnit> work_units;

		GlobalState(DatabaseInstance &db, const std::string &api_key, std::vector<WorkUnit> work_units);

		idx_t MaxThreads() const override {
			return 1; // Single-threaded for rate limiting
		}

		// Fetch the next work unit (a sequence batch or one assembly), parse it
		// into current_batch, accumulate any missing accessions, and emit a
		// per-batch warning if NCBI dropped any. Returns false when no more units.
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
