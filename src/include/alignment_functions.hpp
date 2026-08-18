#pragma once

#include "duckdb/function/aggregate_function.hpp"
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

class CigarQueryCoverageFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

class CigarQueryIntervalsFunction {
public:
	static void Register(ExtensionLoader &loader);
	static ScalarFunction GetFunction();
};

// Aggregate rather than scalar: identity for a read the aligner split into several records
// has to be computed from all of them at once.
class CigarPooledIdentityFunction {
public:
	static void Register(ExtensionLoader &loader);
	static AggregateFunction GetFunction();
};

} // namespace duckdb
