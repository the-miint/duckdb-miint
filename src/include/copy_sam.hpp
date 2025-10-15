#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"

namespace duckdb {

class CopySAMFunction {
public:
	static CopyFunction GetFunction();
};

} // namespace duckdb
