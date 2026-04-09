#pragma once

#include "CoverageDepthCalculator.hpp"
#include "duckdb.hpp"
#include "duckdb/function/aggregate_function.hpp"

namespace duckdb {

class ComputeCoverageDepthFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
