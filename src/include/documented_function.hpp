#pragma once

#include "duckdb/common/vector.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

#include <initializer_list>
#include <string>
#include <utility>
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
inline void RegisterDocumentedScalar(ExtensionLoader &loader, ScalarFunction function, const std::string &description,
                                     std::initializer_list<const char *> parameter_names,
                                     const std::vector<std::string> &examples, const std::string &alias_of = "",
                                     std::initializer_list<const char *> categories = {}) {
	auto parameter_types = function.arguments;
	CreateScalarFunctionInfo info(std::move(function));
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	FunctionDescription fd;
	fd.description = description;
	for (const auto *name : parameter_names) {
		fd.parameter_names.emplace_back(name);
	}
	fd.parameter_types = std::move(parameter_types);
	for (const auto &ex : examples) {
		fd.examples.emplace_back(ex);
	}
	for (const auto *c : categories) {
		fd.categories.emplace_back(c);
	}
	info.descriptions.push_back(std::move(fd));
	loader.RegisterFunction(std::move(info));
}

// Variant for ScalarFunctionSet (multiple overloads under one name). Each
// FunctionDescription in the descriptions vector matches one overload by
// parameter_types; duckdb_functions() routes the right description to the
// right variant via CalcDescriptionSpecificity.
inline void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set,
                                        std::initializer_list<FunctionDescription> per_overload_descriptions,
                                        const std::string &alias_of = "") {
	CreateScalarFunctionInfo info(std::move(set));
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	for (const auto &d : per_overload_descriptions) {
		info.descriptions.push_back(d);
	}
	loader.RegisterFunction(std::move(info));
}

// Register a table function with documentation routed through duckdb_functions()
// the same way scalar functions are. positional_parameter_names names the
// positional arguments (the named-parameter map on the TableFunction is
// surfaced separately by the catalog and is not covered by this list).
//
// IMPORTANT: named-parameter defaults are not catalog-accessible via any
// DuckDB API. They live inside each table function's Bind callback. Document
// them in the description prose (e.g. "**Default:** `false`") until the
// proposed miint_function_metadata() debug function lands.
inline void RegisterDocumentedTableFunction(ExtensionLoader &loader, TableFunction function,
                                            const std::string &description,
                                            std::initializer_list<const char *> positional_parameter_names,
                                            const std::vector<std::string> &examples,
                                            const std::string &alias_of = "",
                                            std::initializer_list<const char *> categories = {}) {
	auto positional_types = function.arguments;
	CreateTableFunctionInfo info(std::move(function));
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	FunctionDescription fd;
	fd.description = description;
	for (const auto *name : positional_parameter_names) {
		fd.parameter_names.emplace_back(name);
	}
	fd.parameter_types = std::move(positional_types);
	for (const auto &ex : examples) {
		fd.examples.emplace_back(ex);
	}
	for (const auto *c : categories) {
		fd.categories.emplace_back(c);
	}
	info.descriptions.push_back(std::move(fd));
	loader.RegisterFunction(std::move(info));
}

} // namespace duckdb
