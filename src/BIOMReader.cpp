#include "BIOMReader.hpp"
#include <iostream>

namespace miint {

BIOMReader::BIOMReader(const std::string &path1) {
	try {
		file_handle = H5::H5File(path1, H5F_ACC_RDONLY);

		try {
			ds_indices = file_handle.openDataSet(SAMPLE_INDICES);
			ds_indptr = file_handle.openDataSet(SAMPLE_INDPTR);
			ds_data = file_handle.openDataSet(SAMPLE_DATA);
			ds_samp_ids = file_handle.openDataSet(SAMPLE_IDS);
			ds_obs_ids = file_handle.openDataSet(OBS_IDS);
		} catch (...) {
			// If dataset opening fails, close file before re-throwing
			if (file_handle.getId() >= 0) {
				file_handle.close();
			}
			throw;
		}
	} catch (H5::FileIException &e) {
		throw std::runtime_error("Failed to open HDF5 file: " + std::string(e.getDetailMsg()));
	} catch (H5::DataSetIException &e) {
		throw std::runtime_error("Failed to access dataset: " + std::string(e.getDetailMsg()));
	}
}

BIOMReader::~BIOMReader() {
	// Explicitly close datasets and file to ensure locks are released
	// Check if each object is valid before attempting to close
	if (ds_indices.getId() >= 0) {
		ds_indices.close();
	}
	if (ds_indptr.getId() >= 0) {
		ds_indptr.close();
	}
	if (ds_data.getId() >= 0) {
		ds_data.close();
	}
	if (ds_samp_ids.getId() >= 0) {
		ds_samp_ids.close();
	}
	if (ds_obs_ids.getId() >= 0) {
		ds_obs_ids.close();
	}
	if (file_handle.getId() >= 0) {
		file_handle.close();
	}
}

BIOMTable BIOMReader::read() const {
	return BIOMTable(ds_indices, ds_indptr, ds_data, ds_obs_ids, ds_samp_ids);
}

bool BIOMReader::IsBIOM(const std::string &path) {
	std::string target = "format-version";
	bool valid = false;
	H5::Exception::dontPrint();

	try {
		H5::H5File file(path, H5F_ACC_RDONLY);

		try {
			H5::Group root = file.openGroup("/");
			int num_attrs = root.getNumAttrs();

			for (int i = 0; i < num_attrs; i++) {
				H5::Attribute attr = root.openAttribute(i);
				std::string attr_name = attr.getName();
				if (attr_name == target) {
					int values[2];
					attr.read(H5::PredType::NATIVE_INT, &values);

					if (values[0] == 2) {
						valid = true;
						break;
					}
				}
			}
			// Explicitly close file before returning - ensures lock is released
			if (file.getId() >= 0) {
				file.close();
			}
		} catch (H5::Exception &e) {
			// Ensure file is closed even on error
			if (file.getId() >= 0) {
				file.close();
			}
			throw;
		}
	} catch (H5::Exception &e) {
		return false;
	}
	return valid;
}
} // namespace miint
