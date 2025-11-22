#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"

namespace duckdb {

std::vector<std::string> ParseFilePathsParameter(const Value &input, const std::string &function_name) {
	std::vector<std::string> paths;

	if (input.type().id() == LogicalTypeId::VARCHAR) {
		paths.push_back(input.ToString());
	} else if (input.type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input);
		for (const auto &child : list_children) {
			paths.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException(function_name + ": first argument must be VARCHAR or VARCHAR[]");
	}

	if (paths.empty()) {
		throw InvalidInputException(function_name + ": at least one file path must be provided");
	}

	return paths;
}

bool ParseIncludeFilepathParameter(const named_parameter_map_t &named_parameters) {
	auto fp_param = named_parameters.find("include_filepath");
	if (fp_param != named_parameters.end() && !fp_param->second.IsNull()) {
		return fp_param->second.GetValue<bool>();
	}
	return false;
}

} // namespace duckdb
