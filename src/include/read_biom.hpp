#pragma once
#include "BIOMReader.hpp"
#include "remote_file_helper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {
class ReadBIOMTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> biom_paths;
		bool include_filepath;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
		std::vector<miint::BIOMTableField> fields;

		explicit Data(const std::vector<std::string> &paths, bool include_fp)
		    : biom_paths(paths), include_filepath(include_fp), names({"sample_id", "feature_id", "value"}),
		      types({LogicalType::VARCHAR,  // sample_id
		             LogicalType::VARCHAR,  // feature_id
		             LogicalType::DOUBLE}), // value
		      fields(
		          {miint::BIOMTableField::SAMPLE_ID, miint::BIOMTableField::FEATURE_ID, miint::BIOMTableField::VALUE}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock; // Guards current_file_idx for the work-claim cursor
		// HDF5 serialization lives in miint::g_hdf5_mutex (BIOMReader.hpp) so it
		// covers cross-query races; a per-state mutex only serialized within
		// one read_biom() instance and missed the segfault path from issue #72.
		std::vector<std::string> filepaths;    // Original paths (for include_filepath)
		std::vector<std::string> local_paths;  // Resolved local paths (for BIOMReader)
		miint::ResolvedFileSet resolved_files; // RAII cleanup for temp files
		size_t current_file_idx;

		idx_t MaxThreads() const override {
			// Use conservative fixed-thread approach with 4 threads
			// Each thread processes complete files sequentially from shared queue
			// Bounded parallelism prevents excessive memory usage
			return 4;
		}

		GlobalState(const std::vector<std::string> &original_paths, miint::ResolvedFileSet resolved)
		    : filepaths(original_paths), resolved_files(std::move(resolved)), current_file_idx(0) {
			for (const auto &rf : resolved_files.Files()) {
				local_paths.push_back(rf.local_path);
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		std::string path; // Original path (for include_filepath display)
		miint::BIOMTable table;
		size_t current_row = 0;
		size_t total_rows = 0;
		bool done = false;

		bool GetNextFile(GlobalState &global_state) {
			std::string local_path;
			{
				std::lock_guard<std::mutex> guard(global_state.lock);

				if (global_state.current_file_idx >= global_state.filepaths.size()) {
					done = true;
					return false;
				}

				path = global_state.filepaths[global_state.current_file_idx];
				local_path = global_state.local_paths[global_state.current_file_idx];
				global_state.current_file_idx++;
			}

			// Serialize HDF5 calls against every other reader in the process —
			// see miint::g_hdf5_mutex doc in BIOMReader.hpp. The inner scope
			// keeps the BIOMReader destructor (H5Fclose / H5Dclose) inside
			// the lock; releasing earlier would race with another thread's
			// H5Fopen.
			{
				std::lock_guard<std::mutex> hdf5_guard(miint::g_hdf5_mutex);
				{
					auto reader = miint::BIOMReader(local_path);
					table = reader.read();
				}
			}

			total_rows = table.nnz();
			current_row = 0;

			return true;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);
	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void SetResultVector(Vector &result_vector, const miint::BIOMTableField &field,
	                            const miint::BIOMTable &record, size_t current_row, size_t n_rows);
	static void SetResultVectorString(Vector &result_vector, const miint::BIOMTableField &field,
	                                  const miint::BIOMTable &record, size_t current_row, size_t n_rows);
	static void SetResultVectorDouble(Vector &result_vector, const miint::BIOMTableField &field,
	                                  const miint::BIOMTable &record, size_t current_row, size_t n_rows);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath, size_t num_records);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
