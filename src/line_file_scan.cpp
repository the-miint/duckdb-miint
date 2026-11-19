#include "line_file_scan.hpp"

#include <algorithm>
#include <climits>
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/common/vector/map_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "miint_named_parameter.hpp"

namespace duckdb {

unique_ptr<FunctionData> BindLineScan(ClientContext &context, TableFunctionBindInput &input, vector<Identifier> &names,
                                      vector<LogicalType> &return_types) {
	auto &fs = FileSystem::GetFileSystem(context);
	auto data = make_uniq<LineScanData>();
	const std::string function_name = input.table_function.GetName().GetIdentifierName();
	data->function_name = function_name;

	data->file_paths = ExpandPathArgument(fs, context, input.inputs[0], function_name).paths;
	for (const auto &path : data->file_paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException(function_name + ": stdin is not supported");
		}
	}
	RequireLocalFilesExist(fs, data->file_paths);
	data->include_filepath = ParseIncludeFilepathParameter(input.named_parameters);
	if (data->include_filepath) {
		names.emplace_back("filepath");
		return_types.emplace_back(LogicalType::VARCHAR);
	}
	return std::move(data);
}

namespace {

unique_ptr<GlobalTableFunctionState> InitLineScanGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<LineScanData>();
	idx_t threads = std::max<idx_t>(1, data.file_paths.size());
	return make_uniq<LineScanGlobalState>(threads);
}

// 50 bytes a line is shorter than typical PAF and GFA records, so the estimate errs high, as
// EstimateFileScanCardinality wants.
unique_ptr<NodeStatistics> LineScanCardinality(ClientContext &context, const FunctionData *bind_data) {
	constexpr idx_t BYTES_PER_LINE = 50;
	return EstimateFileScanCardinality(bind_data->Cast<LineScanData>().file_paths, BYTES_PER_LINE);
}

} // namespace

TableFunction MakeLineScanFunction(const std::string &name, table_function_t execute, table_function_bind_t bind,
                                   table_function_init_local_t init_local) {
	TableFunction tf(Identifier(name), {LogicalType::ANY}, execute, bind, InitLineScanGlobal, init_local);
	AddNamedParameter(tf, "include_filepath", LogicalType::BOOLEAN);
	tf.cardinality = LineScanCardinality;
	return tf;
}

LineFile::LineFile(FileSystem &fs, const std::string &path_p, const std::string &function_name)
    : path(path_p), stream(CreateDuckDBSeqStream(fs, path_p), miint::duckdb_seq_close),
      lines([raw = stream.get(), prefix = function_name + ": " + path_p + ": "](char *dst, size_t len) -> size_t {
	      miint::g_seq_read_error.clear();
	      int n = miint::duckdb_seq_read(raw, dst, static_cast<unsigned int>(std::min<size_t>(len, INT_MAX)));
	      if (n < 0) {
		      throw IOException(prefix + miint::g_seq_read_error);
	      }
	      return static_cast<size_t>(n);
      }) {
}

void SetTagMap(Vector &map_vec, idx_t row, const std::vector<miint::SamStyleTag> &tags) {
	auto offset = ListVector::GetListSize(map_vec);
	ListVector::Reserve(map_vec, offset + tags.size());
	auto &keys = MapVector::GetKeys(map_vec);
	auto &values = MapVector::GetValues(map_vec);
	for (size_t i = 0; i < tags.size(); i++) {
		FlatVector::GetDataMutable<string_t>(keys)[offset + i] = StringVector::AddString(keys, tags[i].key);
		FlatVector::GetDataMutable<string_t>(values)[offset + i] = StringVector::AddString(values, tags[i].value);
	}
	auto &entry = FlatVector::GetDataMutable<list_entry_t>(map_vec)[row];
	entry.offset = offset;
	entry.length = tags.size();
	ListVector::SetListSize(map_vec, offset + tags.size());
}

} // namespace duckdb
