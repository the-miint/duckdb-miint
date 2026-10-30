#pragma once
//
// rype_filter(relation, index_path[, minimap2_index_path], threshold := 0.1):
// the rows of `relation` that RYpe does not match and, given a minimap2 index,
// that minimap2 maps neither mate of, every column as the relation has it.
//
// It is the host-depletion steps
//
//   SELECT * FROM reads ANTI JOIN rype_classify(index, 'reads', threshold := t) USING (read_id)
//   ... ANTI JOIN align_minimap2(<what RYpe kept>, index_path := minimap2_index_path) USING (read_id)
//
// in one pass: the relation is read once, by a filter stream
// (rype_input_stream.hpp), in groups of rows. Each group is one RYpe call. The
// rows RYpe returns no hit for are read from the chunks the stream held for
// them, mapped by minimap2 there, and the rows it does not map are emitted from
// those chunks. RYpe needs no read_id, and nothing is read a second time to join
// on it. minimap2 runs only on short reads (see LONG_READ_MEAN_LENGTH in
// rype_filter.cpp).

#include "Minimap2Aligner.hpp"
#include "rype.h"
#include "rype_common.hpp"
#include "rype_input_stream.hpp"

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parallel/task_executor.hpp"

#include <memory>
#include <string>

namespace duckdb {

class RypeFilterTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string relation;
		std::string index_path;
		//! Empty when no minimap2 index was given.
		std::string minimap2_index_path;
		double threshold = RYPE_DEFAULT_THRESHOLD;
		//! The relation's columns at bind: the output schema, and what the scan projects.
		vector<string> column_names;
		vector<LogicalType> column_types;
		//! Where minimap2's inputs are among them. sequence2 and read_id are optional.
		idx_t sequence1_column = DConstants::INVALID_INDEX;
		idx_t sequence2_column = DConstants::INVALID_INDEX;
		idx_t read_id_column = DConstants::INVALID_INDEX;
	};

	//! What DuckDB's threads and the background jobs share; see rype_filter_schedule.hpp.
	struct Schedule;

	struct GlobalState : public GlobalTableFunctionState {
		//! The jobs hold references of their own.
		shared_ptr<Schedule> schedule;
		//! Runs the jobs on DuckDB's async pool. Null when that pool has no threads:
		//! the jobs then run on the DuckDB thread that finds them due.
		unique_ptr<TaskExecutor> executor;
		idx_t max_threads = 1;

		idx_t MaxThreads() const override {
			return max_threads;
		}

		~GlobalState() override;
	};

	struct LocalState : public LocalTableFunctionState {
		//! Set when minimap2 runs and this thread maps the units it emits (no async
		//! threads); attached to the index once it has loaded.
		unique_ptr<miint::Minimap2Aligner> aligner;
		bool attached = false;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
