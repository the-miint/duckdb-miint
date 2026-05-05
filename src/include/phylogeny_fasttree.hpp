#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class PhylogenyFastTreeTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

class PhylogenyFastTreeAvailableScalar {
public:
	static ScalarFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
