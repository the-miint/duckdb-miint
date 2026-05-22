#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ComputePileupTableFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
