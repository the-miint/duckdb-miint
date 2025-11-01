#include "BIOMReader.hpp"
#include <iostream>
#include <sys/stat.h>
#include <cstring>

namespace miint {

BIOMReader::BIOMReader(const std::string &path1) {
	// Disable HDF5's automatic error printing
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

	// Check if file exists and get size
	struct stat st;
	if (stat(path1.c_str(), &st) != 0) {
		throw std::runtime_error("File does not exist or cannot be accessed: " + path1);
	}

	// Use C API for maximum control and error reporting
	herr_t status;
	hid_t file_id = -1;
	hid_t ds_indices_id = -1, ds_indptr_id = -1, ds_data_id = -1;
	hid_t ds_samp_ids_id = -1, ds_obs_ids_id = -1;

	try {
		// Open file with C API
		file_id = H5Fopen(path1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		if (file_id < 0) {
			throw std::runtime_error("H5Fopen failed with return value: " + std::to_string(file_id));
		}

		// Check if file is actually HDF5
		htri_t is_hdf5 = H5Fis_hdf5(path1.c_str());

		// Get file intent
		unsigned intent;
		status = H5Fget_intent(file_id, &intent);

		// Try to open first dataset
		ds_indices_id = H5Dopen2(file_id, SAMPLE_INDICES, H5P_DEFAULT);

		if (ds_indices_id < 0) {
			// Get more detailed error info

			// Check if the dataset path exists by trying to open as object
			htri_t exists = H5Lexists(file_id, SAMPLE_INDICES, H5P_DEFAULT);

			if (exists > 0) {
				// Path exists but can't open as dataset - check what type it is
				H5O_info2_t obj_info;
				status = H5Oget_info_by_name3(file_id, SAMPLE_INDICES, &obj_info, H5O_INFO_BASIC, H5P_DEFAULT);
			}

			H5Fclose(file_id);
			throw std::runtime_error("Failed to open dataset " + std::string(SAMPLE_INDICES));
		}

		ds_indptr_id = H5Dopen2(file_id, SAMPLE_INDPTR, H5P_DEFAULT);
		if (ds_indptr_id < 0)
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_INDPTR));

		ds_data_id = H5Dopen2(file_id, SAMPLE_DATA, H5P_DEFAULT);
		if (ds_data_id < 0)
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_DATA));

		ds_samp_ids_id = H5Dopen2(file_id, SAMPLE_IDS, H5P_DEFAULT);
		if (ds_samp_ids_id < 0)
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_IDS));

		ds_obs_ids_id = H5Dopen2(file_id, OBS_IDS, H5P_DEFAULT);
		if (ds_obs_ids_id < 0)
			throw std::runtime_error("Failed to open " + std::string(OBS_IDS));

		// Store the C handles directly
		file_handle = file_id;
		ds_indices = ds_indices_id;
		ds_indptr = ds_indptr_id;
		ds_data = ds_data_id;
		ds_samp_ids = ds_samp_ids_id;
		ds_obs_ids = ds_obs_ids_id;

	} catch (const std::exception &e) {
		// Clean up any open handles
		if (ds_obs_ids_id >= 0)
			H5Dclose(ds_obs_ids_id);
		if (ds_samp_ids_id >= 0)
			H5Dclose(ds_samp_ids_id);
		if (ds_data_id >= 0)
			H5Dclose(ds_data_id);
		if (ds_indptr_id >= 0)
			H5Dclose(ds_indptr_id);
		if (ds_indices_id >= 0)
			H5Dclose(ds_indices_id);
		if (file_id >= 0)
			H5Fclose(file_id);

		throw;
	}
}

BIOMReader::~BIOMReader() {
	// Explicitly close datasets and file to ensure locks are released
	// Check if each handle is valid before attempting to close
	if (ds_indices >= 0) {
		H5Dclose(ds_indices);
	}
	if (ds_indptr >= 0) {
		H5Dclose(ds_indptr);
	}
	if (ds_data >= 0) {
		H5Dclose(ds_data);
	}
	if (ds_samp_ids >= 0) {
		H5Dclose(ds_samp_ids);
	}
	if (ds_obs_ids >= 0) {
		H5Dclose(ds_obs_ids);
	}
	if (file_handle >= 0) {
		H5Fclose(file_handle);
	}
}

BIOMTable BIOMReader::read() const {
	return BIOMTable(ds_indices, ds_indptr, ds_data, ds_obs_ids, ds_samp_ids);
}

bool BIOMReader::IsBIOM(const std::string &path) {
	// Disable HDF5's automatic error printing
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);

	const char *target = "format-version";
	bool valid = false;

	hid_t file_id = -1;
	hid_t root_id = -1;
	hid_t attr_id = -1;

	// Open file
	file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		return false;
	}

	// Open root group
	root_id = H5Gopen(file_id, "/", H5P_DEFAULT);
	if (root_id < 0) {
		H5Fclose(file_id);
		return false;
	}

	// Get number of attributes
	H5O_info2_t obj_info;
	if (H5Oget_info3(root_id, &obj_info, H5O_INFO_NUM_ATTRS) < 0) {
		H5Gclose(root_id);
		H5Fclose(file_id);
		return false;
	}

	// Check each attribute
	for (hsize_t i = 0; i < obj_info.num_attrs; i++) {
		attr_id = H5Aopen_by_idx(root_id, ".", H5_INDEX_NAME, H5_ITER_NATIVE, i, H5P_DEFAULT, H5P_DEFAULT);
		if (attr_id < 0) {
			continue;
		}

		// Get attribute name
		char attr_name[256];
		ssize_t name_size = H5Aget_name(attr_id, sizeof(attr_name), attr_name);

		if (name_size > 0 && strcmp(attr_name, target) == 0) {
			// Found the format-version attribute, read it
			int values[2];
			if (H5Aread(attr_id, H5T_NATIVE_INT, values) >= 0) {
				if (values[0] == 2) {
					valid = true;
				}
			}
			H5Aclose(attr_id);
			break;
		}

		H5Aclose(attr_id);
	}

	// Clean up
	H5Gclose(root_id);
	H5Fclose(file_id);

	return valid;
}
} // namespace miint
