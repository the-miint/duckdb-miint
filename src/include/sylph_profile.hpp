#pragma once
#ifdef MIINT_HAS_SYLPH

#include "SylphDatabase.hpp"
#include "per_sample_table_function.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

// Pull sylph.h directly so SylphProfileParams / SylphSketchParams are
// available as field types on Data below. This keeps sylph as the single
// source of truth for parameter defaults — Bind() seeds these structs via
// sylph_profile_params_default() / sylph_sketch_params_default() and the
// named-parameter dispatch only overrides individual fields the caller
// explicitly set. (Including sylph.h here is cheap; it's a small C header
// with no transitive deps beyond stdint.h.)
#include "sylph.h"

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace duckdb {

// FracMinHash relative-abundance profiling table function.
//
// Wraps sylph (Shaw & Yu 2024, *Nature Biotechnology*) as a DuckDB table
// function: takes a paired-end reads table/view and a `.syldb` reference
// database path, sketches the reads in-memory via the FFI's streaming
// builder, runs sylph_profile, and returns the resulting Arrow batch as a
// DataChunk.
//
// Output schema is fixed at 9 columns — see `SetupOutputSchema()` below.
//
// Sample-input contract (validated at bind time):
//   read_id   VARCHAR (required)
//   sequence1 VARCHAR (required)
//   sequence2 VARCHAR (optional; NULL/empty allowed for SE)
//
// Database lifecycle:
//   The .syldb is loaded once per `sylph_profile` call into the GlobalState
//   and shared across all Execute invocations on that call. A future
//   refinement could promote this to a session-scoped cache so back-to-back
//   calls against the same syldb skip the load — see PLAN-sylph.md.
class SylphProfileTableFunction {
public:
	// Bind-time data: input table name, syldb path, validated input schema,
	// resolved named parameters, and the fixed output schema.
	struct Data : public TableFunctionData {
		std::string source_table;
		std::string syldb_path;

		// Schema of the source table — drives sketch-builder input below.
		// Always has read_id + sequence1; sequence2 is optional.
		SequenceTableSchema schema;

		vector<std::string> output_names;
		vector<LogicalType> output_types;

		// Sylph FFI parameter structs. Bind() seeds both via
		// sylph_profile_params_default() / sylph_sketch_params_default(), so
		// every default flows from sylph itself. Named parameters override
		// individual fields. No miint-side defaults — sylph is the only
		// source of truth.
		SylphProfileParams profile_params;
		SylphSketchParams sketch_params;

		// miint owns scheduler-aware threading (auto-balance against the
		// DuckDB worker pool). 0 = auto.
		uint32_t user_threads = 0;

		// Per-sample partitioning state. When `has_sample_id`, the function
		// runs once per distinct sample value and the output has a prepended
		// `sample_id` column matching the source column's type.
		bool has_sample_id = false;
		PerSampleBindInfo sample_info;
	};

	// Holds the Arrow batch returned by sylph_profile() and the parsed
	// schema metadata DuckDB's converter needs. Single-sample mode keeps
	// one of these in GlobalState; per-sample mode keeps one per worker
	// in LocalState (each sample's batch lives independently).
	//
	// Lifetime contract (mirrors RYpe's): the wrapper's ArrowArray is
	// reference-counted, so DuckDB Vectors that zero-copy fixed-width
	// columns out of it via `array_states[c].owned_data = current_chunk`
	// keep the buffers alive past Execute(). The `arrow_table` holds
	// pointers into `output_schema`'s data, so destruction order matters
	// — see ArrowOutputState's own destructor + Reset().
	struct ArrowOutputState {
		ArrowSchema output_schema {};
		ArrowTableSchema arrow_table;
		shared_ptr<ArrowArrayWrapper> current_chunk;
		idx_t batch_offset = 0;
		std::unordered_map<idx_t, unique_ptr<ArrowArrayScanState>> array_states;

		~ArrowOutputState();

		// Tear down the current sylph_profile batch + schema so a fresh one
		// (next sample in per-sample mode) can be installed. Safe to call
		// when nothing's been installed yet.
		void Reset();

		// Fetch (creating on first use) the per-column scan state DuckDB's
		// Arrow→DuckDB converter threads through.
		ArrowArrayScanState &GetArrayState(idx_t col_idx, ClientContext &context);
	};

	// Global state owns the loaded SylphDatabase, shared read-only across
	// per-sample threads. Inherits from PerSampleGlobalState so the standard
	// `next_sample_idx` claim machinery works.
	//
	// In single-sample mode, the sketch and Arrow output buffer live here;
	// in per-sample mode they migrate to LocalState (one set per thread).
	struct GlobalState : public PerSampleGlobalState {
		std::unique_ptr<miint::SylphDatabaseHandle> db;

		// Single-sample-mode-only fields; unused in per-sample mode.
		::SylphSketch *sketch = nullptr;
		std::unique_ptr<QuerySequenceStream> stream;
		bool profile_done = false;
		ArrowOutputState arrow;

		~GlobalState() override;
	};

	// Per-thread local state for per-sample mode. Each thread owns its own
	// connection (so it can create a TEMP VIEW filtered to one sample),
	// sketch, and Arrow output buffer. Single-sample mode does not allocate
	// this.
	struct LocalState : public LocalTableFunctionState {
		unique_ptr<Connection> conn;
		std::unique_ptr<QuerySequenceStream> stream;
		::SylphSketch *sketch = nullptr;
		ArrowOutputState arrow;
		Value sample_value; // identity of the sample this LocalState is processing

		~LocalState() override;
	};

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

#endif // MIINT_HAS_SYLPH
