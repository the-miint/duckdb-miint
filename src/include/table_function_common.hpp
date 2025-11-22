#pragma once
#include <string>
#include <vector>
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/named_parameter_map.hpp"

namespace duckdb {

// Parse file paths parameter that can be VARCHAR or VARCHAR[]
// Validates that at least one path is provided
std::vector<std::string> ParseFilePathsParameter(const Value &input, const std::string &function_name);

// Parse include_filepath named parameter (optional BOOLEAN, default false)
bool ParseIncludeFilepathParameter(const named_parameter_map_t &named_parameters);

} // namespace duckdb
