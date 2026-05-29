#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignPairwiseWfa2ScoreFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseWfa2CigarFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseWfa2FullFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
