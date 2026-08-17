#pragma once

#include "CumulativeCoverageAccumulator.hpp"
#include "duckdb.hpp"
#include "duckdb/function/aggregate_function.hpp"

namespace duckdb {

class CumulativeCoverageFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
