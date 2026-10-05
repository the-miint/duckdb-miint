#pragma once

#include "KreppIndexBuilder.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>
#include <vector>

namespace duckdb {

// krepp_index_create(sequence_table, output_path[, tree_table | newick_path][, k, w, h, ...])
//
// Builds a krepp index from a relation, so a caller never needs krepp's CLI on
// PATH. The counterpart to place_krepp: what this writes, that reads.
//
// Every distinct read_id in `sequence_table` becomes one krepp reference, and
// every row carrying that read_id becomes one FASTA record inside it - so a
// multi-contig genome is just several rows sharing a read_id. Those names are
// what the backbone tree's tips must match.
class KreppIndexCreateTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string sequence_table;
		std::string output_path;
		// Exactly one of these is set; Bind rejects both and neither.
		std::string tree_table;
		std::string newick_path;
		miint::KreppIndexOptions options;
		SequenceTableSchema schema;

		vector<std::string> names;
		vector<LogicalType> types;

		Data()
		    : names({"output_path", "k", "w", "h", "num_references", "status"}),
		      types({LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
		             LogicalType::BIGINT, LogicalType::VARCHAR}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		int64_t num_references = 0;
		// Read back from the index krepp wrote rather than recomputed here, so
		// the row reports what was built and not what was asked for.
		int32_t k = 0;
		int32_t w = 0;
		int32_t h = 0;
		bool done = false;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	struct LocalState : public LocalTableFunctionState {};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<Identifier> &names);
	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);
	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
