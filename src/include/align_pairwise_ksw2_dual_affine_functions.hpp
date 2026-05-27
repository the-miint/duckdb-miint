#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignPairwiseKsw2DualAffineScoreFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseKsw2DualAffineCigarFunction {
public:
	static void Register(ExtensionLoader &loader);
};

class AlignPairwiseKsw2DualAffineFullFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
