#pragma once

#include "massql_parser.hpp"

#include <string>

namespace miint {

class MassQLTranspiler {
public:
	static std::string to_sql(const MassQLQuery &query, const std::string &source);
	static bool is_file_path(const std::string &source);
};

} // namespace miint
