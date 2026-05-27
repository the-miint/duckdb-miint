#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignPairwiseKsw2SpliceScoreFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseKsw2SpliceCigarFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseKsw2SpliceFullFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
