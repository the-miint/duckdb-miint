#pragma once

#include "duckdb/common/vector.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/parsed_data/create_function_info.hpp"

#include <initializer_list>
#include <string>
#include <vector>

// Helpers that wrap DuckDB's standard CreateXxxFunctionInfo registration
// path so descriptions, parameter names, examples, aliases, and categories
// flow into FunctionDescription and out through duckdb_functions(). The
// docs site (site/) introspects duckdb_functions() to generate reference
// pages, so registrations that go through these helpers appear in the docs
// automatically with no sidecar files.
//
// alias_of, when non-empty, links a registration to its canonical name
// (e.g. is_paired -> alignment_is_paired). The catalog surfaces it in the
// duckdb_functions().alias_of column.
//
// DuckDB also ships FunctionList::RegisterExtensionFunctions for bulk
// registration of static metadata tables (used by core_functions). It is
// scalar/aggregate-only and assumes a static StaticFunctionDefinition[]
// driven by \1-delimited variant strings — wrong shape for one-off
// extension code, hence these thin wrappers instead.
//
// Bodies live in src/documented_function.cpp; keeping them out of the
// header avoids vague-linkage emission of duckdb static-const class
// members in every TU that includes this file.
namespace duckdb {

void RegisterDocumentedScalar(ExtensionLoader &loader, ScalarFunction function, const std::string &description,
                              std::initializer_list<const char *> parameter_names,
                              const std::vector<std::string> &examples, const std::string &alias_of = "",
                              std::initializer_list<const char *> categories = {});

// Per-overload metadata for ScalarFunctionSet registrations. Pairs each
// overload's positional parameter names with their LogicalTypeId. Using
// LogicalTypeId (the plain enum) instead of LogicalType (the static-const
// class members) avoids ODR-using the duckdb static members in every TU.
struct DocumentedOverload {
	std::initializer_list<const char *> parameter_names;
	std::initializer_list<LogicalTypeId> parameter_types;
};

// Register a ScalarFunctionSet where every overload shares the same
// description, examples, and categories — only parameter names/types
// differ. Builds one FunctionDescription per overload internally.
void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set, const std::string &description,
                                 std::initializer_list<DocumentedOverload> overloads,
                                 const std::vector<std::string> &examples, const std::string &alias_of = "",
                                 std::initializer_list<const char *> categories = {});

// Lower-level ScalarFunctionSet registration for when overloads need
// genuinely different descriptions or examples. Caller builds each
// FunctionDescription by hand.
void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set,
                                 std::initializer_list<FunctionDescription> per_overload_descriptions,
                                 const std::string &alias_of = "");

// Register a table function. positional_parameter_names names the
// positional arguments; the named-parameter map on TableFunction is
// surfaced separately by the catalog.
//
// IMPORTANT: named-parameter defaults are not catalog-accessible via any
// DuckDB API. They live inside each table function's Bind callback.
// Document them in the description prose (the field accepts full markdown,
// so a `### Named parameters` table renders fine) until the proposed
// miint_function_metadata() debug function lands.
void RegisterDocumentedTableFunction(ExtensionLoader &loader, TableFunction function, const std::string &description,
                                     std::initializer_list<const char *> positional_parameter_names,
                                     const std::vector<std::string> &examples, const std::string &alias_of = "",
                                     std::initializer_list<const char *> categories = {});

} // namespace duckdb
