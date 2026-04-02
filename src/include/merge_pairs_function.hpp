#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class MergePairsFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
