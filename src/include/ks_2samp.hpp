#pragma once

#include "KsTwoSample.hpp"
#include "duckdb.hpp"
#include "duckdb/function/scalar_function.hpp"

namespace duckdb {

class KsTwoSampleFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
