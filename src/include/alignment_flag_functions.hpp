#pragma once
#include "duckdb.hpp"

namespace duckdb {

class ExtensionLoader;

class AlignmentFlagFunctions {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
