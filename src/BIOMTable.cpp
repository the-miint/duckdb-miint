#include "BIOMTable.hpp"
#include <vector>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <unordered_set>

namespace miint {

BIOMTable::BIOMTable()
    : coo_values({}), feature_ids_ordered({}), sample_ids_ordered({}), coo_feature_indices({}), coo_sample_indices({}) {
}

BIOMTable::BIOMTable(const H5::DataSet &ds_indices, const H5::DataSet &ds_indptr, const H5::DataSet &ds_data,
                     const H5::DataSet &ds_obs_ids, const H5::DataSet &ds_samp_ids) {
	// adapted from
	// https://github.com/biocore/unifrac-binaries/blob/main/src/biom.cpp
	feature_ids_ordered = load_dataset_1D_str(ds_obs_ids);
	sample_ids_ordered = load_dataset_1D_str(ds_samp_ids);
	std::vector<int32_t> indptr = load_dataset_1D<int32_t>(ds_indptr, H5::PredType::NATIVE_INT32);
	std::vector<int32_t> indices = load_dataset_1D<int32_t>(ds_indices, H5::PredType::NATIVE_INT32);
	coo_values = load_dataset_1D<double>(ds_data, H5::PredType::NATIVE_DOUBLE);

	InitCOOFromCSC(indptr, indices);
}
BIOMTable::BIOMTable(const std::vector<std::string> &feature_ids, const std::vector<std::string> &sample_ids,
                     const std::vector<double> &values)
    : coo_values(values) {
	InitCOOFromCOO(feature_ids, sample_ids);
}

void BIOMTable::InitCOOFromCOO(const std::vector<std::string> &feature_ids,
                               const std::vector<std::string> &sample_ids) {
	feature_ids_ordered = unique_ids_in_order(feature_ids);
	sample_ids_ordered = unique_ids_in_order(sample_ids);
	coo_feature_indices = ids_to_indices(feature_ids, feature_ids_ordered);
	coo_sample_indices = ids_to_indices(sample_ids, sample_ids_ordered);
	compress_coo();
	coo_sample_indices_as_ids = indices_to_ids(coo_sample_indices, sample_ids_ordered);
	coo_feature_indices_as_ids = indices_to_ids(coo_feature_indices, feature_ids_ordered);
}

std::vector<size_t> ids_to_indices(const std::vector<std::string> &ids, const std::vector<std::string> &ordered) {
	if (ids.empty()) {
		return {}; // No IDs to convert
	}

	if (ordered.empty()) {
		throw std::invalid_argument("Cannot convert IDs: ordered is empty");
	}

	std::unordered_map<std::string, size_t> lookup;
	lookup.reserve(ordered.size());
	for (size_t i = 0; i < ordered.size(); i++) {
		lookup[ordered[i]] = i;
	}

	std::vector<size_t> coo_indices(ids.size());
	for (size_t i = 0; i < ids.size(); i++) {
		coo_indices[i] = lookup.at(ids[i]);
	}

	return coo_indices;
}

std::vector<std::string> unique_ids_in_order(const std::vector<std::string> &ids) {
	if (ids.empty()) {
		return {};
	}

	std::vector<std::string> unique_ordered_ids;
	std::unordered_set<std::string> seen;

	// we don't know these actual sizes, so let's reserve some space up front
	// which should cover many cases and avoid many of the lower end of
	// reallocations
	unique_ordered_ids.reserve(std::min(ids.size(), size_t(10000)));
	seen.reserve(std::min(ids.size(), size_t(10000)));
	for (auto &id : ids) {
		if (!seen.contains(id)) {
			unique_ordered_ids.push_back(id);
			seen.insert(id);
		}
	}

	return unique_ordered_ids;
}

void BIOMTable::InitCOOFromCSC(const std::vector<int32_t> &indptr, const std::vector<int32_t> &indices) {
	auto nnz = coo_values.size();
	auto n_cols = indptr.size() - 1;

	coo_feature_indices.reserve(nnz);
	coo_sample_indices.reserve(nnz);

	for (auto col = 0; col < n_cols; col++) {
		auto start = indptr[col];
		auto end = indptr[col + 1];

		for (auto offset = start; offset < end; offset++) {
			auto index = indices[offset];

			coo_sample_indices.push_back(col);
			coo_feature_indices.push_back(index);
		}
	}

	// Skip compress_coo() since BIOM files are already in canonical form
	// (sorted by row then column, no duplicates, no zeros)
	coo_sample_indices_as_ids = indices_to_ids(coo_sample_indices, sample_ids_ordered);
	coo_feature_indices_as_ids = indices_to_ids(coo_feature_indices, feature_ids_ordered);
}

std::vector<std::string> BIOMTable::indices_to_ids(const std::vector<size_t> &indices,
                                                   const std::vector<std::string> &names) {
	std::vector<std::string> result;
	result.reserve(indices.size());
	for (auto index : indices) {
		result.push_back(names[index]);
	}

	return result;
}

const std::vector<std::string> &BIOMTable::COOFeatures() const {
	return coo_feature_indices_as_ids;
}
const std::vector<size_t> &BIOMTable::COOFeatureIndices() const {
	return coo_feature_indices;
}

const std::vector<std::string> &BIOMTable::COOSamples() const {
	return coo_sample_indices_as_ids;
}

const std::vector<size_t> &BIOMTable::COOSampleIndices() const {
	return coo_sample_indices;
}
const std::vector<double> &BIOMTable::COOValues() const {
	return coo_values;
}

std::vector<std::string> BIOMTable::load_dataset_1D_str(const H5::DataSet &ds_ids) {
	// adapted from
	// https://github.com/biocore/unifrac-binaries/blob/main/src/biom.cpp
	H5::DataSpace dataspace;
	try {
		dataspace = ds_ids.getSpace();
	} catch (H5::DataSpaceIException &e) {
		throw std::runtime_error("Failed to access dataspace: " + std::string(e.getDetailMsg()));
	}

	hsize_t dims[1];
	dataspace.getSimpleExtentDims(dims, nullptr);

	if (dims[0] == 0) {
		return std::vector<std::string>();
	}

	// Get and check datatype
	H5::DataType dtype = ds_ids.getDataType();
	H5T_class_t type_class = dtype.getClass();

	// Check if it's a string type
	if (type_class != H5T_STRING) {
		throw std::runtime_error("Dataset is not a string type");
	}

	// Check if it's variable length
	H5::StrType strtype = ds_ids.getStrType();
	if (!strtype.isVariableStr()) {
		throw std::runtime_error("Dataset is not variable-length string type");
	}

	/* the IDs are a dataset of variable length strings */
	std::vector<char *> rdata(dims[0]);
	H5::StrType datatype = ds_ids.getStrType();
	ds_ids.read((void *)rdata.data(), datatype);

	// Convert to std::vector<std::string>
	std::vector<std::string> ids;
	ids.reserve(dims[0]);
	for (hsize_t i = 0; i < dims[0]; i++) {
		ids.emplace_back(rdata[i]);
	}

	// Free memory allocated by HDF5
	ds_ids.vlenReclaim((void *)rdata.data(), datatype, dataspace);

	for (auto &id : ids) {
	}
	return ids;
}

// Helper to get the appropriate HDF5 datatype
template <typename T>
H5::DataType BIOMTable::get_hdf5_type() {
	if constexpr (std::is_same_v<T, double>) {
		return H5::PredType::NATIVE_DOUBLE;
	} else if constexpr (std::is_same_v<T, int32_t>) {
		return H5::PredType::NATIVE_INT32;
	} else {
		throw std::runtime_error("Unsupported type");
	}
}

template <typename T>
std::vector<T> BIOMTable::load_dataset_1D(const H5::DataSet &ds, const H5::PredType &exp_dtype) {
	H5::DataSpace dataspace;
	auto dtype = ds.getDataType();
	auto name = ds.getObjName();

	if (dtype != exp_dtype) {
		throw std::runtime_error("Invalid data type on BIOM load from dataset: " + name);
	}

	if (dtype != get_hdf5_type<T>()) {
		throw std::runtime_error("Type mismatch");
	}

	try {
		dataspace = ds.getSpace();
	} catch (H5::DataSpaceIException &e) {
		throw std::runtime_error("Failed to access dataspace: " + std::string(e.getDetailMsg()));
	}

	hsize_t dims[1];
	dataspace.getSimpleExtentDims(dims, nullptr);

	std::vector<T> data(dims[0]);
	ds.read(data.data(), get_hdf5_type<T>());

	return data;
}

uint32_t BIOMTable::nnz() {
	auto values_size = coo_values.size();
	auto sids_size = coo_sample_indices.size();
	auto fids_size = coo_feature_indices.size();

	if (values_size != sids_size || values_size != fids_size) {
		throw std::runtime_error("Invalid data size on BIOM load");
	}

	return values_size;
}

void sort_by_row_then_column(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values) {
	std::vector<size_t> indices(rows.size());
	std::iota(indices.begin(), indices.end(), 0);

	std::ranges::sort(indices,
	                  [&](size_t i, size_t j) { return std::tie(rows[i], cols[i]) < std::tie(rows[j], cols[j]); });

	// naive approach, but which requires memory allocation
	// std::vector<size_t> sorted_rows(rows.size());
	// std::vector<size_t> sorted_cols(cols.size());
	// std::vector<double> sorted_values(values.size());

	// for (size_t i = 0; i < indices.size(); ++i) {
	//	sorted_rows[i] = rows[indices[i]];
	//	sorted_cols[i] = cols[indices[i]];
	//	sorted_values[i] = values[indices[i]];
	// }

	// rows = std::move(sorted_rows);
	// cols = std::move(sorted_cols);
	// values = std::move(sorted_values);

	// so instead, let's swap elements
	apply_permutation(rows, cols, values, indices);
}

void apply_permutation(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values,
                       const std::vector<size_t> &indices) {
	// swap values in place being careful about cycle detection
	std::vector<bool> visited(indices.size(), false);
	for (size_t start = 0; start < indices.size(); ++start) {
		if (visited[start] || indices[start] == start) {
			continue;
		}

		// Save the starting values
		size_t start_row = rows[start];
		size_t start_col = cols[start];
		double start_val = values[start];

		size_t j = start;
		while (indices[j] != start) {
			visited[j] = true;
			size_t next = indices[j];
			rows[j] = rows[next];
			cols[j] = cols[next];
			values[j] = values[next];
			j = next;
		}
		// Close the cycle
		visited[j] = true;
		rows[j] = start_row;
		cols[j] = start_col;
		values[j] = start_val;
	}
}

void BIOMTable::compress_coo() {
	auto n = coo_sample_indices.size();
	if (n == 0) {
		return;
	}

	sort_by_row_then_column(coo_sample_indices, coo_feature_indices, coo_values);

	std::vector<size_t> coo_sample_indices_compressed;
	std::vector<size_t> coo_feature_indices_compressed;
	std::vector<double> coo_values_compressed;

	coo_sample_indices_compressed.reserve(n);
	coo_feature_indices_compressed.reserve(n);
	coo_values_compressed.reserve(n);

	size_t last_row = coo_sample_indices[0];
	size_t last_col = coo_feature_indices[0];
	double accum_value = coo_values[0];

	for (size_t i = 1; i < n; i++) {
		size_t current_row = coo_sample_indices[i];
		size_t current_col = coo_feature_indices[i];
		double current_value = coo_values[i];

		if (current_row == last_row && current_col == last_col) {
			accum_value += current_value;
		} else {
			if (accum_value > 0.) {
				coo_sample_indices_compressed.push_back(last_row);
				coo_feature_indices_compressed.push_back(last_col);
				coo_values_compressed.push_back(accum_value);
			}

			last_row = current_row;
			last_col = current_col;
			accum_value = current_value;
		}
	}

	if (accum_value > 0.) {
		coo_sample_indices_compressed.push_back(last_row);
		coo_feature_indices_compressed.push_back(last_col);
		coo_values_compressed.push_back(accum_value);
	}

	coo_sample_indices = std::move(coo_sample_indices_compressed);
	coo_feature_indices = std::move(coo_feature_indices_compressed);
	coo_values = std::move(coo_values_compressed);
}

SparseMatrix BIOMTable::ToCSR() const {
	// CSR: features are rows (major axis), samples are columns (minor axis)
	// Sort by (feature, sample) and build row pointers
	return ConvertCOOToCompressed(true);
}

SparseMatrix BIOMTable::ToCSC() const {
	// CSC: samples are columns (major axis), features are rows (minor axis)
	// Sort by (sample, feature) and build column pointers
	// Note: COO is already sorted by (sample, feature) after compress_coo()
	return ConvertCOOToCompressed(false);
}

SparseMatrix BIOMTable::ConvertCOOToCompressed(bool by_row) const {
	// Convert COO to compressed format (CSR or CSC)
	// by_row=true:  CSR (features=rows=major, samples=cols=minor)
	// by_row=false: CSC (samples=cols=major, features=rows=minor)

	auto n_major = by_row ? NumFeatures() : NumSamples();
	auto n_values = coo_values.size();

	if (n_values == 0) {
		return SparseMatrix({}, {}, {0});
	}

	// Create index array for sorting
	std::vector<size_t> sort_indices(n_values);
	std::iota(sort_indices.begin(), sort_indices.end(), 0);

	// Sort by (major_axis, minor_axis)
	if (by_row) {
		// CSR: sort by (feature, sample)
		std::ranges::sort(sort_indices, [&](size_t i, size_t j) {
			return std::tie(coo_feature_indices[i], coo_sample_indices[i]) <
			       std::tie(coo_feature_indices[j], coo_sample_indices[j]);
		});
	} else {
		// CSC: sort by (sample, feature)
		// COO is already sorted this way after compress_coo(), but be explicit
		std::ranges::sort(sort_indices, [&](size_t i, size_t j) {
			return std::tie(coo_sample_indices[i], coo_feature_indices[i]) <
			       std::tie(coo_sample_indices[j], coo_feature_indices[j]);
		});
	}

	// Build compressed arrays
	std::vector<double> data;
	std::vector<int32_t> minor_indices; // column indices for CSR, row indices for CSC
	std::vector<int32_t> indptr;        // row pointers for CSR, column pointers for CSC

	data.reserve(n_values);
	minor_indices.reserve(n_values);
	indptr.reserve(n_major + 1);

	indptr.push_back(0);

	size_t current_major = 0;
	for (size_t i = 0; i < n_values; i++) {
		auto idx = sort_indices[i];
		auto major_axis = by_row ? coo_feature_indices[idx] : coo_sample_indices[idx];
		auto minor_axis = by_row ? coo_sample_indices[idx] : coo_feature_indices[idx];
		auto value = coo_values[idx];

		// Add indptr entries for any empty major axis positions
		while (current_major < major_axis) {
			indptr.push_back(static_cast<int32_t>(data.size()));
			current_major++;
		}

		data.push_back(value);
		minor_indices.push_back(static_cast<int32_t>(minor_axis));
	}

	// Add final indptr entries for any trailing empty major axis positions
	while (current_major < n_major) {
		indptr.push_back(static_cast<int32_t>(data.size()));
		current_major++;
	}

	return SparseMatrix(std::move(data), std::move(minor_indices), std::move(indptr));
}
} // namespace miint
