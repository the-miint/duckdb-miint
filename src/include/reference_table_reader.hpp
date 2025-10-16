#pragma once

#include "duckdb/main/client_context.hpp"
#include <string>
#include <unordered_map>

namespace duckdb {

std::unordered_map<std::string, uint64_t> ReadReferenceTable(ClientContext &context, const std::string &table_name);

} // namespace duckdb
