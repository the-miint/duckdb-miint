#include "documented_function.hpp"

#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

#include <utility>

namespace duckdb {

void RegisterDocumentedScalar(ExtensionLoader &loader, ScalarFunction function, const std::string &description,
                              std::initializer_list<const char *> parameter_names,
                              const std::vector<std::string> &examples, const std::string &alias_of,
                              std::initializer_list<const char *> categories) {
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

void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set,
                                 std::initializer_list<FunctionDescription> per_overload_descriptions,
                                 const std::string &alias_of) {
	CreateScalarFunctionInfo info(std::move(set));
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	for (const auto &d : per_overload_descriptions) {
		info.descriptions.push_back(d);
	}
	loader.RegisterFunction(std::move(info));
}

void RegisterDocumentedTableFunction(ExtensionLoader &loader, TableFunction function, const std::string &description,
                                     std::initializer_list<const char *> positional_parameter_names,
                                     const std::vector<std::string> &examples, const std::string &alias_of,
                                     std::initializer_list<const char *> categories) {
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
