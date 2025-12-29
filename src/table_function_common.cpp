#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/types/value.hpp"
#include <algorithm>

namespace duckdb {

// Check if a path refers to stdin
static bool IsStdinPath(const std::string &path) {
	return path == "-" || path == "/dev/stdin" || path == "/dev/fd/0" || path == "/proc/self/fd/0";
}

// Check if any part of the pattern could match stdin paths
static bool GlobCouldMatchStdin(const std::string &pattern) {
	// Check for patterns that could match stdin-like paths
	// This is conservative - we check if the pattern starts with these prefixes or contains wildcards that could
	if (pattern.find("/dev/std") != std::string::npos || pattern.find("/dev/fd/") != std::string::npos ||
	    pattern.find("/proc/self/fd/") != std::string::npos) {
		return true;
	}
	// Check for dash patterns like "-*" or patterns starting with "-"
	if (!pattern.empty() && pattern[0] == '-') {
		return true;
	}
	return false;
}

std::vector<std::string> ParseFilePathsParameter(const Value &input, const std::string &function_name) {
	std::vector<std::string> paths;

	if (input.type().id() == LogicalTypeId::VARCHAR) {
		paths.push_back(input.ToString());
	} else if (input.type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input);
		for (const auto &child : list_children) {
			paths.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException(function_name + ": first argument must be VARCHAR or VARCHAR[]");
	}

	if (paths.empty()) {
		throw InvalidInputException(function_name + ": at least one file path must be provided");
	}

	return paths;
}

bool ParseIncludeFilepathParameter(const named_parameter_map_t &named_parameters) {
	auto fp_param = named_parameters.find("include_filepath");
	if (fp_param != named_parameters.end() && !fp_param->second.IsNull()) {
		return fp_param->second.GetValue<bool>();
	}
	return false;
}

std::vector<std::string> ExpandGlobPattern(FileSystem &fs, ClientContext &context, const std::string &pattern) {
	// Check if this is a glob pattern
	if (!FileSystem::HasGlob(pattern)) {
		// Not a glob pattern, return as-is (validation happens elsewhere)
		return {pattern};
	}

	// Only check for stdin paths when it's actually a glob pattern
	if (GlobCouldMatchStdin(pattern)) {
		throw InvalidInputException("Glob patterns cannot include stdin paths");
	}

	// Expand the glob pattern
	auto files = fs.GlobFiles(pattern, context, FileGlobOptions::ALLOW_EMPTY);

	if (files.empty()) {
		throw IOException("No files matched pattern: " + pattern);
	}

	// Extract paths and sort alphabetically
	std::vector<std::string> paths;
	paths.reserve(files.size());
	for (const auto &file : files) {
		paths.push_back(file.path);
	}
	std::sort(paths.begin(), paths.end());

	return paths;
}

GlobExpansionResult ExpandGlobPatternWithInfo(FileSystem &fs, ClientContext &context, const std::string &pattern) {
	GlobExpansionResult result;
	result.is_glob = FileSystem::HasGlob(pattern);

	// Check if pattern could match stdin paths
	if (result.is_glob && GlobCouldMatchStdin(pattern)) {
		throw InvalidInputException("Glob patterns cannot include stdin paths");
	}

	if (!result.is_glob) {
		// Not a glob pattern, return as-is
		result.paths = {pattern};
		return result;
	}

	// Expand the glob pattern
	auto files = fs.GlobFiles(pattern, context, FileGlobOptions::ALLOW_EMPTY);

	if (files.empty()) {
		throw IOException("No files matched pattern: " + pattern);
	}

	// Extract paths and sort alphabetically
	result.paths.reserve(files.size());
	for (const auto &file : files) {
		result.paths.push_back(file.path);
	}
	std::sort(result.paths.begin(), result.paths.end());

	return result;
}

} // namespace duckdb
