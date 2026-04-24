#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class AlignmentSeqIdentityFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

class CigarSequenceIdentityFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

class CigarQueryLengthFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

class AlignmentQueryCoverageFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

} // namespace duckdb
