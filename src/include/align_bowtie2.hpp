#pragma once

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// align_bowtie2 routes through the gpl-boundary daemon's bowtie2-align /
// bowtie2-build tools (schema_version 2 / 1 respectively). Subjects are
// indexed via bowtie2-build at InitGlobal time, queries are streamed per
// DuckDB chunk through bowtie2-align in Execute. Index files live under a
// per-query temp dir and are unlinked when the GlobalState destructor runs.
class AlignBowtie2TableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

// Scalar function: returns true when the daemon's tool registry advertises
// both bowtie2-align and bowtie2-build. Implementation in align_bowtie2.cpp.
void RegisterBowtie2AvailableFunction(ExtensionLoader &loader);

} // namespace duckdb
