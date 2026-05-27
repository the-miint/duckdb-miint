#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ComputeMsaConsensusFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
