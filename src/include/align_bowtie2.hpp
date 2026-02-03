#pragma once
#include "Bowtie2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <mutex>
#include <vector>

namespace duckdb {

class AlignBowtie2TableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string subject_table;
		miint::Bowtie2Config config;
		SequenceTableSchema query_schema;
		std::vector<miint::AlignmentSubject> subjects; // Pre-loaded at bind time

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names({"read_id", "flags",          "reference",     "position",        "stop_position", "mapq",
		             "cigar",   "mate_reference", "mate_position", "template_length", "tag_as",        "tag_xs",
		             "tag_ys",  "tag_xn",         "tag_xm",        "tag_xo",          "tag_xg",        "tag_nm",
		             "tag_yt",  "tag_md",         "tag_sa"}),
		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::BIGINT,    // position
		             LogicalType::BIGINT,    // stop_position
		             LogicalType::UTINYINT,  // mapq
		             LogicalType::VARCHAR,   // cigar
		             LogicalType::VARCHAR,   // mate_reference
		             LogicalType::BIGINT,    // mate_position
		             LogicalType::BIGINT,    // template_length
		             LogicalType::BIGINT,    // tag_as
		             LogicalType::BIGINT,    // tag_xs
		             LogicalType::BIGINT,    // tag_ys
		             LogicalType::BIGINT,    // tag_xn
		             LogicalType::BIGINT,    // tag_xm
		             LogicalType::BIGINT,    // tag_xo
		             LogicalType::BIGINT,    // tag_xg
		             LogicalType::BIGINT,    // tag_nm
		             LogicalType::VARCHAR,   // tag_yt
		             LogicalType::VARCHAR,   // tag_md
		             LogicalType::VARCHAR})  // tag_sa
		{
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::unique_ptr<miint::Bowtie2Aligner> aligner;
		idx_t current_query_offset;
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset;
		bool done;
		bool finished_aligning; // True after finish() called on aligner

		idx_t MaxThreads() const override {
			// Single-threaded - bowtie2 has internal threading via -p
			return 1;
		}

		GlobalState() : current_query_offset(0), buffer_offset(0), done(false), finished_aligning(false) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		LocalState() {
		}
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

// Scalar function to check if bowtie2 is available in PATH
void RegisterBowtie2AvailableFunction(ExtensionLoader &loader);

} // namespace duckdb
