#include "read_newick.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <zlib.h>
#include <fstream>
#include <memory>
#include <sstream>

namespace duckdb {

// RAII wrapper for gzFile to prevent resource leaks
struct GzFileDeleter {
	void operator()(gzFile f) const {
		if (f) {
			gzclose(f);
		}
	}
};
using GzFilePtr = std::unique_ptr<gzFile_s, GzFileDeleter>;

// Buffer size for gzip decompression (16KB recommended by zlib)
static constexpr size_t GZIP_BUFFER_SIZE = 16384;

std::string ReadNewickTableFunction::ReadNewickFile(const std::string &path) {
	// Handle stdin
	if (IsStdinPath(path)) {
		std::stringstream buffer;
		buffer << std::cin.rdbuf();
		return buffer.str();
	}

	// Handle gzipped files
	if (IsGzipped(path)) {
		GzFilePtr gz(gzopen(path.c_str(), "rb"));
		if (!gz) {
			throw IOException("Failed to open gzipped file: " + path);
		}

		std::string content;
		char buf[GZIP_BUFFER_SIZE];
		int bytes_read;
		while ((bytes_read = gzread(gz.get(), buf, sizeof(buf))) > 0) {
			content.append(buf, bytes_read);
		}

		if (bytes_read < 0) {
			int err;
			const char *error_msg = gzerror(gz.get(), &err);
			std::string msg = error_msg ? error_msg : "unknown error";
			throw IOException("Error reading gzipped file: " + path + " - " + msg);
		}

		// gz automatically closed by unique_ptr
		return content;
	}

	// Regular file
	std::ifstream file(path, std::ios::binary);
	if (!file) {
		throw IOException("File not found: " + path);
	}

	std::stringstream buffer;
	buffer << file.rdbuf();
	return buffer.str();
}

std::vector<ReadNewickTableFunction::NodeRow> ReadNewickTableFunction::TreeToRows(const miint::NewickTree &tree) {
	std::vector<NodeRow> rows;
	rows.reserve(tree.num_nodes());

	for (uint32_t i = 0; i < tree.num_nodes(); ++i) {
		NodeRow row;
		row.node_index = static_cast<int64_t>(i);
		row.name = tree.name(i);
		row.branch_length = tree.branch_length(i);
		row.edge_id = tree.edge_id(i);

		uint32_t parent = tree.parent(i);
		if (parent == miint::NewickTree::NO_PARENT) {
			row.parent_index = std::nullopt;
		} else {
			row.parent_index = static_cast<int64_t>(parent);
		}

		row.is_tip = tree.is_tip(i);
		rows.push_back(std::move(row));
	}

	return rows;
}

unique_ptr<FunctionData> ReadNewickTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                       vector<LogicalType> &return_types, vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> file_paths;

	// Handle VARCHAR (single path, potentially a glob) or VARCHAR[] (array of literal paths)
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		file_paths = std::move(result.paths);
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			file_paths.push_back(child.ToString());
		}
		if (file_paths.empty()) {
			throw InvalidInputException("read_newick: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_newick: first argument must be VARCHAR or VARCHAR[]");
	}

	// Detect stdin usage
	bool uses_stdin = false;
	if (file_paths.size() == 1 && IsStdinPath(file_paths[0])) {
		uses_stdin = true;
		// Normalize stdin to /dev/stdin
		if (file_paths[0] == "-") {
			file_paths[0] = "/dev/stdin";
		}
	} else if (file_paths.size() > 1) {
		for (const auto &path : file_paths) {
			if (IsStdinPath(path)) {
				throw InvalidInputException(
				    "stdin ('-' or '/dev/stdin') must be a single file path, not part of an array");
			}
		}
	}

	// Validate files exist (skip stdin)
	for (const auto &path : file_paths) {
		if (!IsStdinPath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	// Parse include_filepath parameter
	bool include_filepath = ParseIncludeFilepathParameter(input.named_parameters);

	auto data = duckdb::make_uniq<Data>(file_paths, include_filepath, uses_stdin);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadNewickTableFunction::InitGlobal(ClientContext &context,
                                                                         TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	return duckdb::make_uniq<GlobalState>(data.file_paths, data.uses_stdin);
}

unique_ptr<LocalTableFunctionState> ReadNewickTableFunction::InitLocal(ExecutionContext &context,
                                                                       TableFunctionInitInput &input,
                                                                       GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadNewickTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	idx_t output_idx = 0;

	while (output_idx < STANDARD_VECTOR_SIZE) {
		// If we have remaining rows from current file, output them
		if (local_state.current_row_idx < local_state.current_rows.size()) {
			size_t rows_to_output = std::min<size_t>(STANDARD_VECTOR_SIZE - output_idx,
			                                         local_state.current_rows.size() - local_state.current_row_idx);

			for (size_t i = 0; i < rows_to_output; ++i) {
				const auto &row = local_state.current_rows[local_state.current_row_idx + i];

				// node_index
				FlatVector::GetData<int64_t>(output.data[0])[output_idx + i] = row.node_index;

				// name (empty string if not specified, never NULL)
				auto &name_vec = output.data[1];
				FlatVector::GetData<string_t>(name_vec)[output_idx + i] = StringVector::AddString(name_vec, row.name);

				// branch_length (nullable - NaN becomes NULL)
				auto &bl_vec = output.data[2];
				if (std::isnan(row.branch_length)) {
					FlatVector::Validity(bl_vec).SetInvalid(output_idx + i);
				} else {
					FlatVector::GetData<double>(bl_vec)[output_idx + i] = row.branch_length;
				}

				// edge_id (nullable)
				auto &edge_vec = output.data[3];
				if (row.edge_id.has_value()) {
					FlatVector::GetData<int64_t>(edge_vec)[output_idx + i] = row.edge_id.value();
				} else {
					FlatVector::Validity(edge_vec).SetInvalid(output_idx + i);
				}

				// parent_index (nullable)
				auto &parent_vec = output.data[4];
				if (row.parent_index.has_value()) {
					FlatVector::GetData<int64_t>(parent_vec)[output_idx + i] = row.parent_index.value();
				} else {
					FlatVector::Validity(parent_vec).SetInvalid(output_idx + i);
				}

				// is_tip
				FlatVector::GetData<bool>(output.data[5])[output_idx + i] = row.is_tip;

				// filepath (if included)
				if (bind_data.include_filepath) {
					auto &fp_vec = output.data[6];
					FlatVector::GetData<string_t>(fp_vec)[output_idx + i] =
					    StringVector::AddString(fp_vec, local_state.current_filepath);
				}
			}

			output_idx += rows_to_output;
			local_state.current_row_idx += rows_to_output;
			continue;
		}

		// Need to load next file
		if (!local_state.has_file) {
			std::lock_guard<std::mutex> guard(global_state.lock);

			if (global_state.next_file_idx >= global_state.file_paths.size()) {
				// No more files
				break;
			}

			local_state.current_file_idx = global_state.next_file_idx;
			global_state.next_file_idx++;
			local_state.has_file = true;
		}

		// Load and parse current file
		const std::string &path = global_state.file_paths[local_state.current_file_idx];
		local_state.current_filepath = path;

		try {
			std::string content = ReadNewickFile(path);
			auto tree = miint::NewickTree::parse(content);
			local_state.current_rows = TreeToRows(tree);
			local_state.current_row_idx = 0;
		} catch (const std::exception &e) {
			throw IOException("Error parsing newick file '" + path + "': " + e.what());
		}

		local_state.has_file = false; // Ready to claim next file when done
	}

	output.SetCardinality(output_idx);
}

TableFunction ReadNewickTableFunction::GetFunction() {
	auto tf = TableFunction("read_newick", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadNewickTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
