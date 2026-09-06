#pragma once

#include "KreppPlacer.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace duckdb {

// place_krepp(query_table := ..., index_path := ...[, newick_path := ...])
//
// Places query sequences onto a phylogeny with krepp, one row per candidate
// edge. Column names follow read_jplace where the two overlap, but the schemas
// are not union-compatible - see the note on KreppPlacement.
class PlaceKreppTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string index_path;
		std::string newick_path;
		miint::KreppConfig config;
		SequenceTableSchema schema;

		vector<std::string> names;
		vector<LogicalType> types;

		Data()
		    : names({"fragment", "edge_num", "likelihood", "like_weight_ratio", "distal_length", "pendant_length",
		             "distance"}),
		      types({LogicalType::VARCHAR, LogicalType::BIGINT, LogicalType::DOUBLE, LogicalType::DOUBLE,
		             LogicalType::DOUBLE, LogicalType::DOUBLE, LogicalType::DOUBLE}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::shared_ptr<miint::SharedKreppIndex> index;
		unique_ptr<QuerySequenceStream> stream;
		idx_t max_threads = 1;
		// One short-read warning for the whole scan. Global rather than
		// per-thread so a scan on N workers does not warn N times.
		std::atomic<bool> warned_short {false};

		idx_t MaxThreads() const override {
			return max_threads;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		unique_ptr<miint::KreppPlacer> placer;
		// Placements for the batch in flight. One query yields a variable number
		// of rows, so a batch rarely lands on a DataChunk boundary.
		std::vector<miint::KreppPlacement> pending;
		idx_t emitted = 0;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *gstate);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
