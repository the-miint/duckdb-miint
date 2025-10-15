#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include <zlib.h>

namespace duckdb {

class CopyFastqFunction {
public:
	static CopyFunction GetFunction();
};

} // namespace duckdb
