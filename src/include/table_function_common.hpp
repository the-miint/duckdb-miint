#pragma once
#include <string>
#include <vector>
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/named_parameter_map.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/client_context.hpp"

namespace duckdb {

// Parse file paths parameter that can be VARCHAR or VARCHAR[]
// Validates that at least one path is provided
std::vector<std::string> ParseFilePathsParameter(const Value &input, const std::string &function_name);

// Parse include_filepath named parameter (optional BOOLEAN, default false)
bool ParseIncludeFilepathParameter(const named_parameter_map_t &named_parameters);

// Expand a glob pattern into a sorted list of file paths
// - If pattern contains glob characters (*, ?, []), expands and sorts alphabetically
// - If pattern is a literal path, returns it as-is
// - Throws IOException if glob matches zero files
// - Throws InvalidInputException if pattern could match stdin paths
std::vector<std::string> ExpandGlobPattern(FileSystem &fs, ClientContext &context, const std::string &pattern);

// Result of glob expansion with metadata
struct GlobExpansionResult {
	std::vector<std::string> paths;
	bool is_glob; // Whether the original pattern was a glob
};

// Expand glob pattern and return whether it was a glob pattern
// Useful for paired-end read validation where both must be globs or both literals
GlobExpansionResult ExpandGlobPatternWithInfo(FileSystem &fs, ClientContext &context, const std::string &pattern);

} // namespace duckdb
