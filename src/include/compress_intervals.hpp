#pragma once

#include "IntervalCompressor.hpp"
#include "duckdb.hpp"
#include "duckdb/function/aggregate_function.hpp"

namespace duckdb {

class CompressIntervalsFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
