#include "read_newick.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <zlib.h>
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

static constexpr size_t READ_BUFFER_SIZE = 65536;

// Read entire contents of a FileHandle into a string
static std::string ReadFileHandleToString(FileHandle &handle) {
	std::string content;
	char buf[READ_BUFFER_SIZE];
	while (true) {
		auto n = handle.Read(buf, sizeof(buf));
		if (n <= 0) {
			break;
		}
		content.append(buf, static_cast<size_t>(n));
	}
	return content;
}

// Read gzipped data from a FileHandle via zlib inflate into a string
static std::string InflateFileHandleToString(FileHandle &handle, const std::string &path) {
	z_stream zs = {};
	// 16 + MAX_WBITS: auto-detect gzip/zlib header
	if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
		throw IOException("Failed to initialize zlib for: " + path);
	}

	std::string content;
	char compressed_buf[READ_BUFFER_SIZE];
	char decompressed_buf[READ_BUFFER_SIZE];

	bool input_eof = false;
	while (true) {
		// Refill compressed buffer if needed
		if (zs.avail_in == 0 && !input_eof) {
			auto n = handle.Read(compressed_buf, sizeof(compressed_buf));
			if (n <= 0) {
				input_eof = true;
			} else {
				zs.avail_in = static_cast<uInt>(n);
				zs.next_in = reinterpret_cast<Bytef *>(compressed_buf);
			}
		}

		zs.avail_out = sizeof(decompressed_buf);
		zs.next_out = reinterpret_cast<Bytef *>(decompressed_buf);

		int ret = inflate(&zs, Z_NO_FLUSH);
		size_t produced = sizeof(decompressed_buf) - zs.avail_out;
		if (produced > 0) {
			content.append(decompressed_buf, produced);
		}

		if (ret == Z_STREAM_END) {
			break;
		}
		if (ret != Z_OK) {
			inflateEnd(&zs);
			throw IOException("Error decompressing gzipped file: " + path + " (zlib error " + std::to_string(ret) +
			                  ")");
		}
		if (input_eof && zs.avail_in == 0 && produced == 0) {
			inflateEnd(&zs);
			throw IOException("Truncated gzip stream: " + path);
		}
	}

	inflateEnd(&zs);
	return content;
}

std::string ReadNewickTableFunction::ReadNewickFile(FileSystem &fs, const std::string &path) {
	// Handle stdin
	if (IsStdinPath(path)) {
		std::stringstream buffer;
		buffer << std::cin.rdbuf();
		return buffer.str();
	}

	auto handle = fs.OpenFile(path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));

	if (IsGzipped(path)) {
		return InflateFileHandleToString(*handle, path);
	}

	return ReadFileHandleToString(*handle);
}

void ReadNewickTableFunction::GetSchema(std::vector<std::string> &names, std::vector<LogicalType> &types,
                                        bool include_filepath) {
	names = {"node_index", "name", "branch_length", "edge_id", "parent_index", "is_tip"};
	types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::DOUBLE,
	         LogicalType::BIGINT, LogicalType::BIGINT,  LogicalType::BOOLEAN};
	if (include_filepath) {
		names.emplace_back("filepath");
		types.emplace_back(LogicalType::VARCHAR);
	}
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

	// Validate files exist (skip stdin and remote paths)
	for (const auto &path : file_paths) {
		if (!IsStdinPath(path) && !miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
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
	auto &fs = FileSystem::GetFileSystem(context);
	return duckdb::make_uniq<GlobalState>(data.file_paths, data.uses_stdin, fs);
}

unique_ptr<LocalTableFunctionState> ReadNewickTableFunction::InitLocal(ExecutionContext &context,
                                                                       TableFunctionInitInput &input,
                                                                       GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadNewickTableFunction::EmitNodeRows(const std::vector<NodeRow> &rows, size_t start_idx, size_t count,
                                           DataChunk &output, idx_t output_offset, bool include_filepath,
                                           const std::string &filepath) {
	for (size_t i = 0; i < count; ++i) {
		const auto &row = rows[start_idx + i];

		// node_index
		FlatVector::GetData<int64_t>(output.data[0])[output_offset + i] = row.node_index;

		// name (empty string if not specified, never NULL)
		auto &name_vec = output.data[1];
		FlatVector::GetData<string_t>(name_vec)[output_offset + i] = StringVector::AddString(name_vec, row.name);

		// branch_length (nullable - NaN becomes NULL)
		auto &bl_vec = output.data[2];
		if (std::isnan(row.branch_length)) {
			FlatVector::Validity(bl_vec).SetInvalid(output_offset + i);
		} else {
			FlatVector::GetData<double>(bl_vec)[output_offset + i] = row.branch_length;
		}

		// edge_id (nullable)
		auto &edge_vec = output.data[3];
		if (row.edge_id.has_value()) {
			FlatVector::GetData<int64_t>(edge_vec)[output_offset + i] = row.edge_id.value();
		} else {
			FlatVector::Validity(edge_vec).SetInvalid(output_offset + i);
		}

		// parent_index (nullable)
		auto &parent_vec = output.data[4];
		if (row.parent_index.has_value()) {
			FlatVector::GetData<int64_t>(parent_vec)[output_offset + i] = row.parent_index.value();
		} else {
			FlatVector::Validity(parent_vec).SetInvalid(output_offset + i);
		}

		// is_tip
		FlatVector::GetData<bool>(output.data[5])[output_offset + i] = row.is_tip;

		// filepath (if included)
		if (include_filepath) {
			auto &fp_vec = output.data[6];
			FlatVector::GetData<string_t>(fp_vec)[output_offset + i] = StringVector::AddString(fp_vec, filepath);
		}
	}
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

			EmitNodeRows(local_state.current_rows, local_state.current_row_idx, rows_to_output, output, output_idx,
			             bind_data.include_filepath, local_state.current_filepath);

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
			std::string content = ReadNewickFile(global_state.fs, path);
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
