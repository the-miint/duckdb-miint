#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include <algorithm>
#include <cstring>

namespace duckdb {

// Check if a path refers to stdin
bool IsStdinPath(const std::string &path) {
	return path == "-" || path == "/dev/stdin" || path == "/dev/fd/0" || path == "/proc/self/fd/0";
}

// Check if a path has .gz extension
bool IsGzipped(const std::string &path) {
	return path.size() >= 3 && path.substr(path.size() - 3) == ".gz";
}

// Check if any part of the pattern could match stdin paths
static bool GlobCouldMatchStdin(const std::string &pattern) {
	// Check for patterns that could match stdin-like paths
	// This is conservative - we check if the pattern starts with these prefixes or contains wildcards that could
	if (pattern.find("/dev/std") != std::string::npos || pattern.find("/dev/fd/") != std::string::npos ||
	    pattern.find("/proc/self/fd/") != std::string::npos) {
		return true;
	}
	// Check for dash patterns like "-*" or patterns starting with "-"
	if (!pattern.empty() && pattern[0] == '-') {
		return true;
	}
	return false;
}

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

std::vector<std::string> ExpandGlobPattern(FileSystem &fs, ClientContext &context, const std::string &pattern) {
	// Check if this is a glob pattern
	if (!FileSystem::HasGlob(pattern)) {
		// Not a glob pattern, return as-is (validation happens elsewhere)
		return {pattern};
	}

	// Only check for stdin paths when it's actually a glob pattern
	if (GlobCouldMatchStdin(pattern)) {
		throw InvalidInputException("Glob patterns cannot include stdin paths");
	}

	// Expand the glob pattern
	auto files = fs.GlobFiles(pattern, context, FileGlobOptions::ALLOW_EMPTY);

	if (files.empty()) {
		throw IOException("No files matched pattern: " + pattern);
	}

	// Extract paths and sort alphabetically
	std::vector<std::string> paths;
	paths.reserve(files.size());
	for (const auto &file : files) {
		paths.push_back(file.path);
	}
	std::sort(paths.begin(), paths.end());

	return paths;
}

GlobExpansionResult ExpandGlobPatternWithInfo(FileSystem &fs, ClientContext &context, const std::string &pattern) {
	GlobExpansionResult result;
	result.is_glob = FileSystem::HasGlob(pattern);

	// Check if pattern could match stdin paths
	if (result.is_glob && GlobCouldMatchStdin(pattern)) {
		throw InvalidInputException("Glob patterns cannot include stdin paths");
	}

	if (!result.is_glob) {
		// Not a glob pattern, return as-is
		result.paths = {pattern};
		return result;
	}

	// Expand the glob pattern
	auto files = fs.GlobFiles(pattern, context, FileGlobOptions::ALLOW_EMPTY);

	if (files.empty()) {
		throw IOException("No files matched pattern: " + pattern);
	}

	// Extract paths and sort alphabetically
	result.paths.reserve(files.size());
	for (const auto &file : files) {
		result.paths.push_back(file.path);
	}
	std::sort(result.paths.begin(), result.paths.end());

	return result;
}

// --- Result vector helpers ---

void SetResultVectorNull(Vector &result_vector) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	result_vector.SetValue(0, Value());
}

void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
	}
}

// Empty strings are treated as NULL. This convention is used for fields like
// spectrum_type, polarity, and activation_method where an empty string indicates
// the value was absent in the source data.
void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
		if (!values[j].empty()) {
			validity.SetValid(j);
		}
	}
}

void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values) {
	auto result_data = FlatVector::GetData<uint8_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values) {
	auto result_data = FlatVector::GetData<uint16_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values,
                                  const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
		if (valid[j]) {
			validity.SetValid(j);
		}
	}
}

void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values, uint8_t qual_offset) {
	idx_t total_child_elements = 0;
	for (auto &qual : values) {
		total_child_elements += qual.size();
	}

	ListVector::Reserve(result_vector, total_child_elements);
	ListVector::SetListSize(result_vector, total_child_elements);

	auto &child_vector = ListVector::GetEntry(result_vector);
	auto child_data = FlatVector::GetData<uint8_t>(child_vector);
	auto list_entries = FlatVector::GetData<list_entry_t>(result_vector);

	const auto output_count = values.size();
	idx_t value_offset = 0;
	for (idx_t row_offset = 0; row_offset < output_count; row_offset++) {
		auto len = values[row_offset].size();
		list_entries[row_offset].offset = value_offset;
		list_entries[row_offset].length = len;

		values[row_offset].write_decoded(child_data + value_offset, qual_offset);
		value_offset += len;
	}

	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllValid(output_count);

	auto &child_validity = FlatVector::Validity(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

void SetResultVectorInt32(Vector &result_vector, const std::vector<int32_t> &values) {
	auto result_data = FlatVector::GetData<int32_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt32Nullable(Vector &result_vector, const std::vector<int32_t> &values,
                                  const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetData<int32_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
		if (valid[j]) {
			validity.SetValid(j);
		}
	}
}

void SetResultVectorDouble(Vector &result_vector, const std::vector<double> &values) {
	auto result_data = FlatVector::GetData<double>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorDoubleNullable(Vector &result_vector, const std::vector<double> &values,
                                   const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetData<double>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
		if (valid[j]) {
			validity.SetValid(j);
		}
	}
}

void SetResultVectorListDouble(Vector &result_vector, const std::vector<std::vector<double>> &values) {
	idx_t total_child_elements = 0;
	for (const auto &vec : values) {
		total_child_elements += vec.size();
	}

	ListVector::Reserve(result_vector, total_child_elements);
	ListVector::SetListSize(result_vector, total_child_elements);

	auto &child_vector = ListVector::GetEntry(result_vector);
	auto child_data = FlatVector::GetData<double>(child_vector);
	auto list_entries = FlatVector::GetData<list_entry_t>(result_vector);

	idx_t value_offset = 0;
	for (idx_t row_offset = 0; row_offset < values.size(); row_offset++) {
		auto len = values[row_offset].size();
		list_entries[row_offset].offset = value_offset;
		list_entries[row_offset].length = len;

		if (len > 0) {
			std::memcpy(child_data + value_offset, values[row_offset].data(), len * sizeof(double));
		}
		value_offset += len;
	}

	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllValid(values.size());

	auto &child_validity = FlatVector::Validity(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

} // namespace duckdb
