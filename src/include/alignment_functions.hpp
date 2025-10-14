#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignmentSeqIdentityFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

} // namespace duckdb
