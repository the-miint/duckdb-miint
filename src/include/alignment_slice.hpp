#pragma once

#include "AlignmentSlicer.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/query_result.hpp"

#include <string>
#include <vector>

namespace duckdb {

class AlignmentSliceTableFunction {
public:
	struct ColumnInfo {
		std::string name;
		LogicalType type;
		bool required;           // cigar, position, stop_position are required
		bool is_tag;             // tag columns are NULLed after trimming
		bool invalidate_on_trim; // template_length is NULLed on any trim
		bool is_id;              // read_id/reference/mate_reference: preserve input id type, don't coerce to VARCHAR
	};

	// Per-output-column role, precomputed at bind time for fast dispatch in Execute
	enum class ColRole : uint8_t {
		PASS_THROUGH,      // copy from input unchanged
		CIGAR,             // set from slicer result
		POSITION,          // set from slicer result
		STOP_POSITION,     // set from slicer result
		SEQUENCE,          // set from slicer result (with NULL handling)
		QUAL,              // set from slicer result (with NULL handling)
		TAG_OR_INVALIDATE, // NULL when trimmed, pass through otherwise
	};

	struct Data : public TableFunctionData {
		std::string table_name;
		int64_t region_start;
		int64_t region_stop;
		bool include_deletions;

		// Output schema (only recognized columns present in input)
		vector<string> output_names;
		vector<LogicalType> output_types;

		// Per-output-column metadata, parallel to output_names:
		vector<int> output_input_idx; // index into the SELECT result columns
		vector<ColRole> output_roles; // precomputed role for fast dispatch in Execute

		// Indices into the SELECT result for columns needed by the slicer
		int select_cigar_idx = -1;
		int select_pos_idx = -1;
		int select_stop_pos_idx = -1;
		int select_seq_idx = -1;
		int select_qual_idx = -1;
		int select_ref_idx = -1;

		// The SELECT query to execute (built at bind time, selects only needed columns)
		string select_query;
	};

	struct GlobalState : public GlobalTableFunctionState {
		// Connection must outlive the StreamQueryResult (it streams lazily)
		unique_ptr<Connection> conn;
		unique_ptr<QueryResult> query_result;
		unique_ptr<DataChunk> current_chunk;
		idx_t chunk_offset = 0;        // current position within current_chunk
		bool stream_exhausted = false; // true after Fetch() returns null

		// The slicer instance
		unique_ptr<miint::AlignmentSlicer> slicer;

		// Single-reference enforcement state
		string seen_reference;
		bool has_seen_reference = false;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void Register(ExtensionLoader &loader);

private:
	static const vector<ColumnInfo> &GetRecognizedColumns();
};

} // namespace duckdb
