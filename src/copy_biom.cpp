#include "copy_biom.hpp"
#include "BIOMTable.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include "H5Cpp.h"
#include <chrono>
#include <iomanip>
#include <sstream>
#include <limits>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct BiomCopyBindData : public FunctionData {
	string file_path;
	string id = "No Table ID";
	string generated_by = "miint";
	bool use_compression = true; // Default to gzip compression
	vector<string> names;
	idx_t feature_id_idx = DConstants::INVALID_INDEX;
	idx_t sample_id_idx = DConstants::INVALID_INDEX;
	idx_t value_idx = DConstants::INVALID_INDEX;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<BiomCopyBindData>();
		result->file_path = file_path;
		result->id = id;
		result->generated_by = generated_by;
		result->use_compression = use_compression;
		result->names = names;
		result->feature_id_idx = feature_id_idx;
		result->sample_id_idx = sample_id_idx;
		result->value_idx = value_idx;
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<BiomCopyBindData>();
		return file_path == other.file_path && id == other.id && generated_by == other.generated_by &&
		       use_compression == other.use_compression && names == other.names &&
		       feature_id_idx == other.feature_id_idx && sample_id_idx == other.sample_id_idx &&
		       value_idx == other.value_idx;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> BiomCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                             const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<BiomCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;

	// Find required columns
	for (idx_t i = 0; i < names.size(); i++) {
		if (names[i] == "feature_id") {
			result->feature_id_idx = i;
		} else if (names[i] == "sample_id") {
			result->sample_id_idx = i;
		} else if (names[i] == "value") {
			result->value_idx = i;
		}
	}

	// Validate required columns exist
	if (result->feature_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT BIOM requires 'feature_id' column");
	}
	if (result->sample_id_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT BIOM requires 'sample_id' column");
	}
	if (result->value_idx == DConstants::INVALID_INDEX) {
		throw BinderException("COPY FORMAT BIOM requires 'value' column");
	}

	// Validate column types
	if (sql_types[result->feature_id_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("COPY FORMAT BIOM: 'feature_id' must be VARCHAR");
	}
	if (sql_types[result->sample_id_idx].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("COPY FORMAT BIOM: 'sample_id' must be VARCHAR");
	}
	if (sql_types[result->value_idx].id() != LogicalTypeId::DOUBLE) {
		throw BinderException("COPY FORMAT BIOM: 'value' must be DOUBLE");
	}

	// Parse optional parameters
	for (auto &option : input.info.options) {
		auto &key = option.first;
		auto &values = option.second;

		if (StringUtil::CIEquals(key, "compression")) {
			string comp_str = StringUtil::Lower(values[0].ToString());
			if (comp_str == "gzip" || comp_str == "gz") {
				result->use_compression = true;
			} else if (comp_str == "none") {
				result->use_compression = false;
			} else {
				throw InvalidInputException("COPY FORMAT BIOM: compression must be 'gzip', 'gz', or 'none'");
			}
		} else if (StringUtil::CIEquals(key, "id")) {
			result->id = values[0].ToString();
		} else if (StringUtil::CIEquals(key, "generated_by")) {
			result->generated_by = values[0].ToString();
		} else {
			throw BinderException("Unknown option for COPY FORMAT BIOM: %s", key);
		}
	}

	return result;
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct BiomCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	// Incremental dictionary building for optimal performance
	unordered_map<string, int32_t> feature_map;
	unordered_map<string, int32_t> sample_map;
	vector<string> feature_ids_ordered;
	vector<string> sample_ids_ordered;
	// Store integer indices instead of strings
	vector<size_t> feature_indices;
	vector<size_t> sample_indices;
	vector<double> values;
};

static unique_ptr<GlobalFunctionData> BiomCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                               const string &file_path) {
	auto gstate = make_uniq<BiomCopyGlobalState>();
	return gstate;
}

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct BiomCopyLocalState : public LocalFunctionData {
	// Store strings locally - dictionary building happens in Combine
	vector<string> local_feature_ids;
	vector<string> local_sample_ids;
	vector<double> local_values;
};

static unique_ptr<LocalFunctionData> BiomCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	auto lstate = make_uniq<BiomCopyLocalState>();
	return lstate;
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void BiomCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                         LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<BiomCopyBindData>();
	auto &gstate = gstate_p.Cast<BiomCopyGlobalState>();
	auto &lstate = lstate_p.Cast<BiomCopyLocalState>();

	// Process each row
	UnifiedVectorFormat feature_data, sample_data, value_data;

	input.data[fdata.feature_id_idx].ToUnifiedFormat(input.size(), feature_data);
	input.data[fdata.sample_id_idx].ToUnifiedFormat(input.size(), sample_data);
	input.data[fdata.value_idx].ToUnifiedFormat(input.size(), value_data);

	auto feature_strings = UnifiedVectorFormat::GetData<string_t>(feature_data);
	auto sample_strings = UnifiedVectorFormat::GetData<string_t>(sample_data);
	auto value_doubles = UnifiedVectorFormat::GetData<double>(value_data);

	for (idx_t row = 0; row < input.size(); row++) {
		auto feature_idx = feature_data.sel->get_index(row);
		auto sample_idx = sample_data.sel->get_index(row);
		auto value_idx = value_data.sel->get_index(row);

		// Check for NULLs
		if (!feature_data.validity.RowIsValid(feature_idx)) {
			throw InvalidInputException("COPY FORMAT BIOM: NULL values not allowed in feature_id column");
		}
		if (!sample_data.validity.RowIsValid(sample_idx)) {
			throw InvalidInputException("COPY FORMAT BIOM: NULL values not allowed in sample_id column");
		}
		if (!value_data.validity.RowIsValid(value_idx)) {
			throw InvalidInputException("COPY FORMAT BIOM: NULL values not allowed in value column");
		}

		// Get values
		string feature_id = feature_strings[feature_idx].GetString();
		string sample_id = sample_strings[sample_idx].GetString();
		double value = value_doubles[value_idx];

		// Validate feature_id and sample_id are not empty
		if (feature_id.empty()) {
			throw InvalidInputException("COPY FORMAT BIOM: empty feature_id not allowed");
		}
		if (sample_id.empty()) {
			throw InvalidInputException("COPY FORMAT BIOM: empty sample_id not allowed");
		}

		// Store in local state (no locking needed - each thread has its own local state)
		lstate.local_feature_ids.push_back(std::move(feature_id));
		lstate.local_sample_ids.push_back(std::move(sample_id));
		lstate.local_values.push_back(value);
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void BiomCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                            LocalFunctionData &lstate_p) {
	auto &gstate = gstate_p.Cast<BiomCopyGlobalState>();
	auto &lstate = lstate_p.Cast<BiomCopyLocalState>();

	// Convert local strings to indices using global dictionaries
	// Lock only during dictionary updates and global state append
	lock_guard<mutex> glock(gstate.lock);

	idx_t n = lstate.local_feature_ids.size();
	vector<size_t> local_feature_indices;
	vector<size_t> local_sample_indices;
	local_feature_indices.reserve(n);
	local_sample_indices.reserve(n);

	for (idx_t i = 0; i < n; i++) {
		const auto &feature_id = lstate.local_feature_ids[i];
		const auto &sample_id = lstate.local_sample_ids[i];

		// Look up or insert feature_id
		auto feature_it = gstate.feature_map.find(feature_id);
		int32_t feature_index;
		if (feature_it == gstate.feature_map.end()) {
			feature_index = static_cast<int32_t>(gstate.feature_ids_ordered.size());
			gstate.feature_map[feature_id] = feature_index;
			gstate.feature_ids_ordered.push_back(feature_id);
		} else {
			feature_index = feature_it->second;
		}

		// Look up or insert sample_id
		auto sample_it = gstate.sample_map.find(sample_id);
		int32_t sample_index;
		if (sample_it == gstate.sample_map.end()) {
			sample_index = static_cast<int32_t>(gstate.sample_ids_ordered.size());
			gstate.sample_map[sample_id] = sample_index;
			gstate.sample_ids_ordered.push_back(sample_id);
		} else {
			sample_index = sample_it->second;
		}

		local_feature_indices.push_back(static_cast<size_t>(feature_index));
		local_sample_indices.push_back(static_cast<size_t>(sample_index));
	}

	// Append local indices to global
	gstate.feature_indices.insert(gstate.feature_indices.end(), local_feature_indices.begin(),
	                              local_feature_indices.end());
	gstate.sample_indices.insert(gstate.sample_indices.end(), local_sample_indices.begin(),
	                             local_sample_indices.end());
	gstate.values.insert(gstate.values.end(), lstate.local_values.begin(), lstate.local_values.end());
}

//===--------------------------------------------------------------------===//
// HDF5 Writing Helpers
//===--------------------------------------------------------------------===//
static string GetCurrentTimestamp() {
	auto now = std::chrono::system_clock::now();
	auto time_t_now = std::chrono::system_clock::to_time_t(now);
	auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()) % 1000;

	std::ostringstream oss;
	oss << std::put_time(std::gmtime(&time_t_now), "%Y-%m-%dT%H:%M:%S");
	oss << "." << std::setfill('0') << std::setw(6) << (ms.count() * 1000); // Convert ms to microseconds
	return oss.str();
}

static void WriteRootAttributes(H5::H5File &file, size_t nnz, size_t n_features, size_t n_samples, const string &id,
                                const string &generated_by) {
	H5::Group root = file.openGroup("/");

	// creation-date
	string creation_date = GetCurrentTimestamp();
	H5::StrType str_type(H5::PredType::C_S1, creation_date.size());
	H5::DataSpace scalar_space(H5S_SCALAR);
	H5::Attribute attr_date = root.createAttribute("creation-date", str_type, scalar_space);
	attr_date.write(str_type, creation_date);

	// format-url
	string format_url = "http://biom-format.org";
	H5::StrType url_type(H5::PredType::C_S1, format_url.size());
	H5::Attribute attr_url = root.createAttribute("format-url", url_type, scalar_space);
	attr_url.write(url_type, format_url);

	// format-version
	hsize_t dims[1] = {2};
	H5::DataSpace version_space(1, dims);
	H5::Attribute attr_version = root.createAttribute("format-version", H5::PredType::NATIVE_INT, version_space);
	int version[2] = {2, 1};
	attr_version.write(H5::PredType::NATIVE_INT, version);

	// generated-by
	H5::StrType gen_type(H5::PredType::C_S1, generated_by.size());
	H5::Attribute attr_gen = root.createAttribute("generated-by", gen_type, scalar_space);
	attr_gen.write(gen_type, generated_by);

	// id
	H5::StrType id_type(H5::PredType::C_S1, id.size());
	H5::Attribute attr_id = root.createAttribute("id", id_type, scalar_space);
	attr_id.write(id_type, id);

	// nnz
	H5::Attribute attr_nnz = root.createAttribute("nnz", H5::PredType::NATIVE_INT64, scalar_space);
	int64_t nnz_val = static_cast<int64_t>(nnz);
	attr_nnz.write(H5::PredType::NATIVE_INT64, &nnz_val);

	// shape
	H5::DataSpace shape_space(1, dims);
	H5::Attribute attr_shape = root.createAttribute("shape", H5::PredType::NATIVE_INT64, shape_space);
	int64_t shape[2] = {static_cast<int64_t>(n_features), static_cast<int64_t>(n_samples)};
	attr_shape.write(H5::PredType::NATIVE_INT64, shape);

	// type (empty string)
	string type_str = "";
	H5::StrType type_type(H5::PredType::C_S1, 1);
	H5::Attribute attr_type = root.createAttribute("type", type_type, scalar_space);
	attr_type.write(type_type, type_str);
}

static void WriteStringDataset(H5::Group &group, const string &name, const std::vector<std::string> &data,
                               bool use_compression) {
	if (data.empty()) {
		// Write empty dataset
		hsize_t dims[1] = {0};
		H5::DataSpace dataspace(1, dims);
		H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);
		H5::DataSet dataset = group.createDataSet(name, datatype, dataspace);
		return;
	}

	// Create variable-length string type
	H5::StrType datatype(H5::PredType::C_S1, H5T_VARIABLE);

	// Create dataspace
	hsize_t dims[1] = {data.size()};
	H5::DataSpace dataspace(1, dims);

	// Create dataset with optional compression
	H5::DSetCreatPropList plist;
	if (use_compression && data.size() > 0) {
		hsize_t chunk_dims[1] = {std::min(data.size(), size_t(65536))};
		plist.setChunk(1, chunk_dims);
		plist.setDeflate(4); // gzip compression level 4 (h5py/BIOM default)
	}

	H5::DataSet dataset = group.createDataSet(name, datatype, dataspace, plist);

	// Convert to char* array for HDF5
	std::vector<const char *> c_strs(data.size());
	for (size_t i = 0; i < data.size(); i++) {
		c_strs[i] = data[i].c_str();
	}

	dataset.write(c_strs.data(), datatype);
}

template <typename T>
static void WriteNumericDataset(H5::Group &group, const string &name, const std::vector<T> &data,
                                const H5::DataType &dtype, bool use_compression) {
	if (data.empty()) {
		// Write empty dataset
		hsize_t dims[1] = {0};
		H5::DataSpace dataspace(1, dims);
		H5::DataSet dataset = group.createDataSet(name, dtype, dataspace);
		return;
	}

	hsize_t dims[1] = {data.size()};
	H5::DataSpace dataspace(1, dims);

	// Create dataset with optional compression
	H5::DSetCreatPropList plist;
	if (use_compression && data.size() > 0) {
		hsize_t chunk_dims[1] = {std::min(data.size(), size_t(65536))};
		plist.setChunk(1, chunk_dims);
		plist.setDeflate(4); // gzip compression level 4 (h5py/BIOM default)
	}

	H5::DataSet dataset = group.createDataSet(name, dtype, dataspace, plist);
	dataset.write(data.data(), dtype);
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void BiomCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &fdata = bind_data.Cast<BiomCopyBindData>();
	auto &gstate = gstate_p.Cast<BiomCopyGlobalState>();

	// Create BIOMTable using pre-computed integer indices and dictionaries
	// This skips the expensive string hashing that the old constructor did
	miint::BIOMTable table(gstate.feature_indices, gstate.sample_indices, gstate.values,
	                       std::move(gstate.feature_ids_ordered), std::move(gstate.sample_ids_ordered));

	auto n_features = table.NumFeatures();
	auto n_samples = table.NumSamples();
	auto nnz = table.nnz();

	// Validate int32 range
	if (n_features > static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
		throw InvalidInputException("COPY FORMAT BIOM: Too many features (>2.1B), exceeds int32_t limit");
	}
	if (n_samples > static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
		throw InvalidInputException("COPY FORMAT BIOM: Too many samples (>2.1B), exceeds int32_t limit");
	}

	// Convert to CSR and CSC
	auto csr = table.ToCSR(); // observation/matrix
	auto csc = table.ToCSC(); // sample/matrix

	// Create HDF5 file
	H5::H5File file(fdata.file_path, H5F_ACC_TRUNC);

	// Write root attributes
	WriteRootAttributes(file, nnz, n_features, n_samples, fdata.id, fdata.generated_by);

	// Create groups
	H5::Group obs_group = file.createGroup("/observation");
	H5::Group obs_matrix_group = file.createGroup("/observation/matrix");
	H5::Group obs_metadata_group = file.createGroup("/observation/metadata");
	H5::Group obs_group_metadata_group = file.createGroup("/observation/group-metadata");

	H5::Group samp_group = file.createGroup("/sample");
	H5::Group samp_matrix_group = file.createGroup("/sample/matrix");
	H5::Group samp_metadata_group = file.createGroup("/sample/metadata");
	H5::Group samp_group_metadata_group = file.createGroup("/sample/group-metadata");

	// Write observation (feature) IDs
	WriteStringDataset(obs_group, "ids", table.FeatureIDs(), fdata.use_compression);

	// Write observation/matrix (CSR)
	WriteNumericDataset(obs_matrix_group, "data", csr.data, H5::PredType::NATIVE_DOUBLE, fdata.use_compression);
	WriteNumericDataset(obs_matrix_group, "indices", csr.indices, H5::PredType::NATIVE_INT32, fdata.use_compression);
	WriteNumericDataset(obs_matrix_group, "indptr", csr.indptr, H5::PredType::NATIVE_INT32, fdata.use_compression);

	// Write sample IDs
	WriteStringDataset(samp_group, "ids", table.SampleIDs(), fdata.use_compression);

	// Write sample/matrix (CSC)
	WriteNumericDataset(samp_matrix_group, "data", csc.data, H5::PredType::NATIVE_DOUBLE, fdata.use_compression);
	WriteNumericDataset(samp_matrix_group, "indices", csc.indices, H5::PredType::NATIVE_INT32, fdata.use_compression);
	WriteNumericDataset(samp_matrix_group, "indptr", csc.indptr, H5::PredType::NATIVE_INT32, fdata.use_compression);

	// Close file
	file.close();
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
CopyFunction CopyBiomFunction::GetFunction() {
	CopyFunction func("biom");
	func.copy_to_bind = BiomCopyBind;
	func.copy_to_initialize_local = BiomCopyInitializeLocal;
	func.copy_to_initialize_global = BiomCopyInitializeGlobal;
	func.copy_to_sink = BiomCopySink;
	func.copy_to_combine = BiomCopyCombine;
	func.copy_to_finalize = BiomCopyFinalize;
	return func;
}

void CopyBiomFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
