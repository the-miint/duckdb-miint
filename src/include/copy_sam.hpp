#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class CopySAMFunction {
public:
	static CopyFunction GetFunction();
	static CopyFunction GetBAMFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
