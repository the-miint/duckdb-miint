#pragma once
#include <string>
#include <vector>
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/named_parameter_map.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/client_context.hpp"
#include "QualScore.hpp"

namespace duckdb {

// Check if a path refers to stdin (-, /dev/stdin, /dev/fd/0, /proc/self/fd/0)
bool IsStdinPath(const std::string &path);

// Check if a path has .gz extension (for gzip-compressed files)
bool IsGzipped(const std::string &path);

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

// --- Result vector helpers ---
// Shared utilities for populating DuckDB result vectors from batch data.
// Used by read_fastx, read_alignments, read_sequences_sam, and other table functions.

void SetResultVectorNull(Vector &result_vector);
void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values);
void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values);
void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath);
void SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values);
void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values);
void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values);
void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values,
                                  const std::vector<bool> &valid);
void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values, uint8_t qual_offset);

} // namespace duckdb
