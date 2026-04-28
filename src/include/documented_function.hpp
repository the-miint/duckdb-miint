#pragma once

#include "duckdb/common/vector.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/parsed_data/create_function_info.hpp"

#include <initializer_list>
#include <string>
#include <vector>

namespace duckdb {

// Register a scalar function whose description, parameter names, and examples
// are exposed via duckdb_functions(). The doc-site generator queries
// duckdb_functions() to build the reference pages, so registrations that go
// through this helper appear in the docs automatically and stay in sync with
// the C++ as a side effect of building the extension.
//
// alias_of, when non-empty, links this registration to its canonical name
// (e.g. is_paired -> alignment_is_paired). The catalog surfaces it in the
// duckdb_functions().alias_of column.
//
// Bodies live in src/documented_function.cpp; keeping them out of the header
// avoids vague-linkage emission of duckdb static-const class members in every
// translation unit that includes this file.
void RegisterDocumentedScalar(ExtensionLoader &loader, ScalarFunction function, const std::string &description,
                              std::initializer_list<const char *> parameter_names,
                              const std::vector<std::string> &examples, const std::string &alias_of = "",
                              std::initializer_list<const char *> categories = {});

// Variant for ScalarFunctionSet (multiple overloads under one name). Each
// FunctionDescription in the descriptions vector matches one overload by
// parameter_types; duckdb_functions() routes the right description to the
// right variant via CalcDescriptionSpecificity.
void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set,
                                 std::initializer_list<FunctionDescription> per_overload_descriptions,
                                 const std::string &alias_of = "");

// Register a table function with documentation routed through duckdb_functions()
// the same way scalar functions are. positional_parameter_names names the
// positional arguments (the named-parameter map on the TableFunction is
// surfaced separately by the catalog and is not covered by this list).
//
// IMPORTANT: named-parameter defaults are not catalog-accessible via any
// DuckDB API. They live inside each table function's Bind callback. Document
// them in the description prose (e.g. "**Default:** `false`") until the
// proposed miint_function_metadata() debug function lands.
void RegisterDocumentedTableFunction(ExtensionLoader &loader, TableFunction function, const std::string &description,
                                     std::initializer_list<const char *> positional_parameter_names,
                                     const std::vector<std::string> &examples, const std::string &alias_of = "",
                                     std::initializer_list<const char *> categories = {});

} // namespace duckdb
