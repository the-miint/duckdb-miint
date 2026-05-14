#include "read_jplace_newick.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <cctype>

namespace duckdb {

// Convert jplace square-bracket edge numbers [n] to curly-brace syntax {n}
// that the Newick parser understands. The Newick parser treats [...] as comments
// and skips them, but jplace files use [integer] for edge numbering.
// Only converts brackets whose content is a non-negative integer (digits only).
// Non-integer bracket content (e.g., [some comment]) is left as-is for the
// Newick parser to handle as a standard comment.
// Note: jplace edge numbers are non-negative per spec (Matsen et al. 2012).
static std::string ConvertBracketEdgeIds(const std::string &newick) {
	std::string result;
	result.reserve(newick.size());

	size_t i = 0;
	while (i < newick.size()) {
		if (newick[i] == '[') {
			// Check if content is a non-negative integer
			size_t start = i + 1;
			size_t j = start;
			while (j < newick.size() && std::isdigit(static_cast<unsigned char>(newick[j]))) {
				++j;
			}
			if (j > start && j < newick.size() && newick[j] == ']') {
				// It's [digits] — convert to {digits}
				result += '{';
				result.append(newick, start, j - start);
				result += '}';
				i = j + 1;
			} else {
				// Not a pure integer — keep as-is (actual comment)
				result += newick[i];
				++i;
			}
		} else {
			result += newick[i];
			++i;
		}
	}

	return result;
}

// Extract the "tree" string value from jplace JSON content.
// The jplace format (Matsen et al., 2012) has a top-level "tree" key whose value
// is always a Newick string. We find this key and extract the JSON string value.
std::string ReadJplaceNewickTableFunction::ExtractTreeFromJplaceContent(const std::string &content,
                                                                        const std::string &path) {
	auto pos = content.find("\"tree\"");
	if (pos == std::string::npos) {
		throw IOException("Jplace file '%s' does not contain a 'tree' field", path);
	}

	// Skip past "tree", whitespace, and colon
	pos += 6;
	while (pos < content.size() &&
	       (content[pos] == ' ' || content[pos] == '\t' || content[pos] == '\n' || content[pos] == '\r')) {
		pos++;
	}
	if (pos >= content.size() || content[pos] != ':') {
		throw IOException("Jplace file '%s': malformed 'tree' field (expected ':')", path);
	}
	pos++;
	while (pos < content.size() &&
	       (content[pos] == ' ' || content[pos] == '\t' || content[pos] == '\n' || content[pos] == '\r')) {
		pos++;
	}
	if (pos >= content.size() || content[pos] != '"') {
		throw IOException("Jplace file '%s': 'tree' field value is not a string", path);
	}

	// Quick-and-dirty JSON string extraction: only handles \\ and \" escapes.
	// Does NOT decode \n, \t, \uXXXX etc. — sufficient for Newick tree strings
	// which contain no such escapes in practice.
	pos++;
	std::string result;
	while (pos < content.size() && content[pos] != '"') {
		if (content[pos] == '\\' && pos + 1 < content.size()) {
			pos++;
			result += content[pos];
		} else {
			result += content[pos];
		}
		pos++;
	}

	if (result.empty()) {
		throw IOException("Jplace file '%s': 'tree' field is empty", path);
	}

	// Convert [integer] edge IDs to {integer} syntax.
	// Many jplace tools (pplacer, SEPP, gappa) use square brackets for edge numbering,
	// but the Newick parser treats [...] as comments. Only pure integers are converted.
	return ConvertBracketEdgeIds(result);
}

ReadJplaceNewickTableFunction::Data::Data(const std::vector<std::string> &paths, bool include_fp)
    : file_paths(paths), include_filepath(include_fp) {
	ReadNewickTableFunction::GetSchema(names, types, include_filepath);
}

unique_ptr<FunctionData> ReadJplaceNewickTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> file_paths;

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_jplace_newick: first argument cannot be NULL");
	}
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		file_paths = std::move(result.paths);
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_jplace_newick: file path list cannot contain NULL");
			}
			file_paths.push_back(child.ToString());
		}
		if (file_paths.empty()) {
			throw InvalidInputException("read_jplace_newick: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_jplace_newick: first argument must be VARCHAR or VARCHAR[]");
	}

	// Validate files exist (skip remote paths)
	for (const auto &path : file_paths) {
		if (!miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	bool include_filepath = ParseIncludeFilepathParameter(input.named_parameters);

	auto data = duckdb::make_uniq<Data>(file_paths, include_filepath);

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadJplaceNewickTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &fs = FileSystem::GetFileSystem(context);
	return duckdb::make_uniq<GlobalState>(data.file_paths, fs);
}

unique_ptr<LocalTableFunctionState> ReadJplaceNewickTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadJplaceNewickTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	idx_t output_idx = 0;

	while (output_idx < STANDARD_VECTOR_SIZE) {
		// Emit remaining rows from current file
		if (local_state.current_row_idx < local_state.current_rows.size()) {
			size_t rows_to_output = std::min<size_t>(STANDARD_VECTOR_SIZE - output_idx,
			                                         local_state.current_rows.size() - local_state.current_row_idx);

			ReadNewickTableFunction::EmitNodeRows(local_state.current_rows, local_state.current_row_idx, rows_to_output,
			                                      output, output_idx, bind_data.include_filepath,
			                                      local_state.current_filepath);

			output_idx += rows_to_output;
			local_state.current_row_idx += rows_to_output;
			continue;
		}

		// Claim next file
		std::string path;
		{
			std::lock_guard<std::mutex> guard(global_state.lock);
			if (global_state.next_file_idx >= global_state.file_paths.size()) {
				break;
			}
			path = global_state.file_paths[global_state.next_file_idx];
			global_state.next_file_idx++;
		}

		local_state.current_filepath = path;

		// Read jplace file as text, extract tree string, parse newick
		try {
			std::string content = ReadNewickTableFunction::ReadNewickFile(global_state.fs, path);
			std::string newick_str = ExtractTreeFromJplaceContent(content, path);
			auto tree = miint::NewickTree::parse(newick_str);
			local_state.current_rows = ReadNewickTableFunction::TreeToRows(tree);
			local_state.current_row_idx = 0;
		} catch (const IOException &) {
			throw;
		} catch (const std::exception &e) {
			throw IOException("Error parsing tree from jplace file '%s': %s", path, e.what());
		}
	}

	output.SetCardinality(output_idx);
}

TableFunction ReadJplaceNewickTableFunction::GetFunction() {
	auto tf = TableFunction("read_jplace_newick", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadJplaceNewickTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
