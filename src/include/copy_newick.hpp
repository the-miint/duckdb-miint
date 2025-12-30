#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class CopyNewickFunction {
public:
	static CopyFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
