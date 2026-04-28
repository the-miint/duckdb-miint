#include "documented_function.hpp"

#include "duckdb/parser/parsed_data/create_scalar_function_info.hpp"
#include "duckdb/parser/parsed_data/create_table_function_info.hpp"

#include <utility>

namespace duckdb {

namespace {

// Build a FunctionDescription from the user-supplied prose and the function's
// own argument types. Shared by the scalar and table-function helpers below;
// they only differ in which CreateXxxFunctionInfo wraps the result.
FunctionDescription BuildDescription(vector<LogicalType> parameter_types, const std::string &description,
                                     std::initializer_list<const char *> parameter_names,
                                     const std::vector<std::string> &examples,
                                     std::initializer_list<const char *> categories) {
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
	return fd;
}

template <class InfoT>
void RegisterWithInfo(ExtensionLoader &loader, InfoT &&info, FunctionDescription &&fd, const std::string &alias_of) {
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	info.descriptions.push_back(std::move(fd));
	loader.RegisterFunction(std::move(info));
}

} // namespace

void RegisterDocumentedScalar(ExtensionLoader &loader, ScalarFunction function, const std::string &description,
                              std::initializer_list<const char *> parameter_names,
                              const std::vector<std::string> &examples, const std::string &alias_of,
                              std::initializer_list<const char *> categories) {
	auto fd = BuildDescription(function.arguments, description, parameter_names, examples, categories);
	RegisterWithInfo(loader, CreateScalarFunctionInfo(std::move(function)), std::move(fd), alias_of);
}

void RegisterDocumentedScalarSet(ExtensionLoader &loader, ScalarFunctionSet set, const std::string &description,
                                 std::initializer_list<DocumentedOverload> overloads,
                                 const std::vector<std::string> &examples, const std::string &alias_of,
                                 std::initializer_list<const char *> categories) {
	CreateScalarFunctionInfo info(std::move(set));
	if (!alias_of.empty()) {
		info.alias_of = alias_of;
	}
	for (const auto &ov : overloads) {
		vector<LogicalType> types;
		for (auto id : ov.parameter_types) {
			types.emplace_back(id);
		}
		info.descriptions.push_back(
		    BuildDescription(std::move(types), description, ov.parameter_names, examples, categories));
	}
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
	auto fd = BuildDescription(function.arguments, description, positional_parameter_names, examples, categories);
	RegisterWithInfo(loader, CreateTableFunctionInfo(std::move(function)), std::move(fd), alias_of);
}

} // namespace duckdb
