#pragma once
#include "MzMLReader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <thread>
#include <vector>

namespace duckdb {
class ReadMzMLTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp)
		    : file_paths(paths), include_filepath(include_fp), names({"spectrum_index",
		                                                              "spectrum_id",
		                                                              "ms_level",
		                                                              "retention_time",
		                                                              "spectrum_type",
		                                                              "polarity",
		                                                              "base_peak_mz",
		                                                              "base_peak_intensity",
		                                                              "total_ion_current",
		                                                              "lowest_mz",
		                                                              "highest_mz",
		                                                              "default_array_length",
		                                                              "precursor_mz",
		                                                              "precursor_charge",
		                                                              "precursor_intensity",
		                                                              "isolation_window_target",
		                                                              "isolation_window_lower",
		                                                              "isolation_window_upper",
		                                                              "activation_method",
		                                                              "collision_energy",
		                                                              "mz_array",
		                                                              "intensity_array",
		                                                              "filter_string",
		                                                              "scan_window_lower",
		                                                              "scan_window_upper",
		                                                              "ms1_scan_index"}),
		      types({LogicalType::INTEGER,
		             LogicalType::VARCHAR,
		             LogicalType::INTEGER,
		             LogicalType::DOUBLE,
		             LogicalType::VARCHAR,
		             LogicalType::VARCHAR,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::INTEGER,
		             LogicalType::DOUBLE,
		             LogicalType::INTEGER,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::VARCHAR,
		             LogicalType::DOUBLE,
		             LogicalType::LIST(LogicalType::DOUBLE),
		             LogicalType::LIST(LogicalType::DOUBLE),
		             LogicalType::VARCHAR,
		             LogicalType::DOUBLE,
		             LogicalType::DOUBLE,
		             LogicalType::INTEGER}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		};
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::vector<std::unique_ptr<miint::MzMLReader>> readers;
		std::vector<std::string> filepaths;
		size_t next_file_idx;

		idx_t MaxThreads() const override {
			auto hw_threads = std::thread::hardware_concurrency();
			if (hw_threads == 0) {
				hw_threads = 1;
			}
			return std::min<idx_t>(filepaths.size(), std::min<idx_t>(8, hw_threads));
		}

		GlobalState(const std::vector<std::string> &paths) : readers(paths.size()), filepaths(paths), next_file_idx(0) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;

		LocalState() : current_file_idx(0), has_file(false) {
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
} // namespace duckdb
