#include "table_function_common.hpp"
#include "remote_file_helper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/common/types/vector.hpp"
#include "duckdb/common/vector/constant_vector.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include <algorithm>
#include <cstring>
#include <filesystem>

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

bool ParseIncludeFilepathParameter(const named_argument_map_t &named_parameters) {
	auto fp_param = named_parameters.find("include_filepath");
	if (fp_param != named_parameters.end() && !fp_param->second.IsNull()) {
		return fp_param->second.GetValue<bool>();
	}
	return false;
}

unique_ptr<NodeStatistics> EstimateFileScanCardinality(const std::vector<std::string> &paths, idx_t bytes_per_row) {
	constexpr idx_t GZIP_EXPANSION_FACTOR = 4;
	constexpr idx_t UNKNOWN_SIZE_ESTIMATE = 10'000'000; // err large

	idx_t total_bytes = 0;
	for (const auto &path : paths) {
		if (IsStdinPath(path) || miint::RemoteFileHelper::IsRemotePath(path)) {
			return make_uniq<NodeStatistics>(UNKNOWN_SIZE_ESTIMATE);
		}
		std::error_code ec;
		auto size = std::filesystem::file_size(path, ec);
		if (ec) {
			return make_uniq<NodeStatistics>(UNKNOWN_SIZE_ESTIMATE);
		}
		total_bytes += IsGzipped(path) ? size * GZIP_EXPANSION_FACTOR : size;
	}
	return make_uniq<NodeStatistics>(std::max<idx_t>(1, total_bytes / bytes_per_row));
}

std::vector<std::string> ExpandGlobPattern(FileSystem &fs, ClientContext &context, const std::string &pattern) {
	return ExpandGlobPatternWithInfo(fs, context, pattern).paths;
}

std::vector<std::string> ExpandPathList(FileSystem &fs, const Value &list, const std::string &function_name,
                                        const std::string &list_name) {
	std::vector<std::string> paths;
	for (const auto &child : ListValue::GetChildren(list)) {
		if (child.IsNull()) {
			throw InvalidInputException(function_name + ": " + list_name + " cannot contain NULL");
		}
		paths.push_back(fs.ExpandPath(child.ToString()));
	}
	if (paths.empty()) {
		throw InvalidInputException(function_name + ": at least one file path must be provided");
	}
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
		// Not a glob pattern: expanded as DuckDB's glob expands one
		result.paths = {fs.ExpandPath(pattern)};
		return result;
	}

	// Expand the glob pattern
	auto files = fs.GlobFiles(pattern, FileGlobOptions::ALLOW_EMPTY);

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

GlobExpansionResult ExpandPathArgument(FileSystem &fs, ClientContext &context, const Value &arg,
                                       const std::string &function_name) {
	if (arg.IsNull()) {
		throw InvalidInputException(function_name + ": first argument cannot be NULL");
	}
	if (arg.type().id() == LogicalTypeId::VARCHAR) {
		return ExpandGlobPatternWithInfo(fs, context, arg.ToString());
	}
	if (arg.type().id() == LogicalTypeId::LIST) {
		GlobExpansionResult result;
		result.paths = ExpandPathList(fs, arg, function_name);
		result.is_glob = false;
		return result;
	}
	throw InvalidInputException(function_name + ": first argument must be VARCHAR or VARCHAR[]");
}

void RequireLocalFilesExist(FileSystem &fs, const std::vector<std::string> &paths) {
	for (const auto &path : paths) {
		if (!IsStdinPath(path) && !miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}
}

std::vector<std::string> SequenceColumnNames() {
	return {"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"};
}

std::vector<LogicalType> SequenceColumnTypes() {
	return {LogicalType::BIGINT,
	        LogicalType::VARCHAR,
	        LogicalType::VARCHAR,
	        LogicalType::VARCHAR,
	        LogicalType::VARCHAR,
	        LogicalType::LIST(LogicalType::UTINYINT),
	        LogicalType::LIST(LogicalType::UTINYINT)};
}

idx_t GetNonNegativeOption(ClientContext &context, const char *name, idx_t if_null, const char *zero_means) {
	Value v;
	// Registered at extension load (LoadInternal in src/miint_extension.cpp).
	// Missing means the option was never registered, which is a load-order
	// programmer error rather than a user mistake.
	if (!context.TryGetCurrentSetting(name, v)) {
		throw InternalException("miint extension option '%s' is not registered — was LoadInternal completed?", name);
	}
	if (v.IsNull()) {
		return if_null;
	}
	const auto value = v.GetValue<int64_t>();
	if (value < 0) {
		throw InvalidInputException("%s must be >= 0 (%s)", name, zero_means);
	}
	return NumericCast<idx_t>(value);
}

// --- Result vector helpers ---

void SetResultVectorNull(Vector &result_vector) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	result_vector.SetValue(0, Value());
}

void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetDataMutable<string_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
	}
}

// Empty strings are treated as NULL. This convention is used for fields like
// spectrum_type, polarity, and activation_method where an empty string indicates
// the value was absent in the source data.
void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetDataMutable<string_t>(result_vector);
	auto &validity = FlatVector::ValidityMutable(result_vector);
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
	auto result_data = FlatVector::GetDataMutable<uint8_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values) {
	auto result_data = FlatVector::GetDataMutable<uint16_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values) {
	auto result_data = FlatVector::GetDataMutable<int64_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values,
                                  const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetDataMutable<int64_t>(result_vector);
	auto &validity = FlatVector::ValidityMutable(result_vector);
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

	auto &child_vector = ListVector::GetChildMutable(result_vector);
	auto child_data = FlatVector::GetDataMutable<uint8_t>(child_vector);
	auto list_entries = FlatVector::GetDataMutable<list_entry_t>(result_vector);

	const auto output_count = values.size();
	idx_t value_offset = 0;
	for (idx_t row_offset = 0; row_offset < output_count; row_offset++) {
		auto len = values[row_offset].size();
		list_entries[row_offset].offset = value_offset;
		list_entries[row_offset].length = len;

		values[row_offset].write_decoded(child_data + value_offset, qual_offset);
		value_offset += len;
	}

	// A mask with no buffer already means all-valid, and giving it one sends every score down the
	// storage append's row-by-row path (docs/internals/duckdb-engine-notes.md), so only a buffer
	// left over from earlier use is rewritten.
	auto &validity = FlatVector::ValidityMutable(result_vector);
	if (validity.CanHaveNull()) {
		validity.SetAllValid(output_count);
	}

	auto &child_validity = FlatVector::ValidityMutable(child_vector);
	if (child_validity.CanHaveNull()) {
		child_validity.SetAllValid(total_child_elements);
	}
}

void SetResultVectorInt32(Vector &result_vector, const std::vector<int32_t> &values) {
	auto result_data = FlatVector::GetDataMutable<int32_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorInt32Nullable(Vector &result_vector, const std::vector<int32_t> &values,
                                  const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetDataMutable<int32_t>(result_vector);
	auto &validity = FlatVector::ValidityMutable(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
		if (valid[j]) {
			validity.SetValid(j);
		}
	}
}

void SetResultVectorDouble(Vector &result_vector, const std::vector<double> &values) {
	auto result_data = FlatVector::GetDataMutable<double>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void SetResultVectorDoubleNullable(Vector &result_vector, const std::vector<double> &values,
                                   const std::vector<bool> &valid) {
	auto result_data = FlatVector::GetDataMutable<double>(result_vector);
	auto &validity = FlatVector::ValidityMutable(result_vector);
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

	auto &child_vector = ListVector::GetChildMutable(result_vector);
	auto child_data = FlatVector::GetDataMutable<double>(child_vector);
	auto list_entries = FlatVector::GetDataMutable<list_entry_t>(result_vector);

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

	auto &validity = FlatVector::ValidityMutable(result_vector);
	validity.SetAllValid(values.size());

	auto &child_validity = FlatVector::ValidityMutable(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

void PopulateSpectrumBatchOutput(DataChunk &output, const miint::MzMLSpectrumBatch &batch, bool include_filepath,
                                 const std::string &filepath) {
	size_t col = 0;
	SetResultVectorInt32(output.data[col++], batch.spectrum_index);
	SetResultVectorString(output.data[col++], batch.spectrum_id);
	SetResultVectorInt32Nullable(output.data[col++], batch.scan_number, batch.scan_number_valid);
	SetResultVectorInt32(output.data[col++], batch.ms_level);
	SetResultVectorDoubleNullable(output.data[col++], batch.retention_time, batch.retention_time_valid);
	SetResultVectorStringNullable(output.data[col++], batch.spectrum_type);
	SetResultVectorStringNullable(output.data[col++], batch.polarity);
	SetResultVectorDoubleNullable(output.data[col++], batch.base_peak_mz, batch.base_peak_mz_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.base_peak_intensity, batch.base_peak_intensity_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.total_ion_current, batch.total_ion_current_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.lowest_mz, batch.lowest_mz_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.highest_mz, batch.highest_mz_valid);
	SetResultVectorInt32(output.data[col++], batch.default_array_length);
	SetResultVectorDoubleNullable(output.data[col++], batch.precursor_mz, batch.precursor_mz_valid);
	SetResultVectorInt32Nullable(output.data[col++], batch.precursor_charge, batch.precursor_charge_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.precursor_intensity, batch.precursor_intensity_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_target,
	                              batch.isolation_window_target_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_lower, batch.isolation_window_lower_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_upper, batch.isolation_window_upper_valid);
	SetResultVectorStringNullable(output.data[col++], batch.activation_method);
	SetResultVectorDoubleNullable(output.data[col++], batch.collision_energy, batch.collision_energy_valid);
	SetResultVectorListDouble(output.data[col++], batch.mz_array);
	SetResultVectorListDouble(output.data[col++], batch.intensity_array);
	SetResultVectorStringNullable(output.data[col++], batch.filter_string);
	SetResultVectorDoubleNullable(output.data[col++], batch.scan_window_lower, batch.scan_window_lower_valid);
	SetResultVectorDoubleNullable(output.data[col++], batch.scan_window_upper, batch.scan_window_upper_valid);
	SetResultVectorInt32Nullable(output.data[col++], batch.ms1_scan_index, batch.ms1_scan_index_valid);
	if (include_filepath) {
		SetResultVectorFilepath(output.data[col++], filepath);
	}
	output.SetChildCardinality(batch.size());
}

void GetListUInt8Slice(Vector &list_vec, UnifiedVectorFormat &list_data, idx_t row_idx, const uint8_t *&out_data,
                       idx_t &out_length) {
	auto list_entries = UnifiedVectorFormat::GetData<list_entry_t>(list_data);
	auto mapped_idx = list_data.sel->get_index(row_idx);
	auto &entry = list_entries[mapped_idx];

	auto &child = ListVector::GetChildMutable(list_vec);
	auto child_data = FlatVector::GetData<uint8_t>(child);

	out_data = child_data + entry.offset;
	out_length = entry.length;
}

void SetStructRowNull(Vector &result, idx_t row_idx) {
	FlatVector::SetNull(result, row_idx, true);
}

} // namespace duckdb
