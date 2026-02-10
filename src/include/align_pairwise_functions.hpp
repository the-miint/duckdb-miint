#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignPairwiseScoreFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseCigarFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseFullFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
