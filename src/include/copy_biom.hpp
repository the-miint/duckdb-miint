#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class CopyBiomFunction {
public:
	static CopyFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
