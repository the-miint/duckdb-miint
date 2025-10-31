#include "BIOMReader.hpp"
#include <iostream>
#include <sys/stat.h>

namespace miint {

BIOMReader::BIOMReader(const std::string &path1) {
	// Log HDF5 version info
	unsigned majnum, minnum, relnum;
	H5get_libversion(&majnum, &minnum, &relnum);
	std::cerr << "DEBUG: HDF5 version " << majnum << "." << minnum << "." << relnum << std::endl;
	std::cerr << "DEBUG: Opening file: " << path1 << std::endl;

	// Check if file exists and get size
	struct stat st;
	if (stat(path1.c_str(), &st) != 0) {
		throw std::runtime_error("File does not exist or cannot be accessed: " + path1);
	}
	std::cerr << "DEBUG: File exists, size: " << st.st_size << " bytes" << std::endl;

	// Use C API for maximum control and error reporting
	herr_t status;
	hid_t file_id = -1;
	hid_t ds_indices_id = -1, ds_indptr_id = -1, ds_data_id = -1;
	hid_t ds_samp_ids_id = -1, ds_obs_ids_id = -1;

	try {
		// Open file with C API
		std::cerr << "DEBUG: Calling H5Fopen..." << std::endl;
		file_id = H5Fopen(path1.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
		std::cerr << "DEBUG: H5Fopen returned: " << file_id << std::endl;

		if (file_id < 0) {
			throw std::runtime_error("H5Fopen failed with return value: " + std::to_string(file_id));
		}

		// Check if file is actually HDF5
		htri_t is_hdf5 = H5Fis_hdf5(path1.c_str());
		std::cerr << "DEBUG: H5Fis_hdf5 returned: " << is_hdf5 << " (1=yes, 0=no, <0=error)" << std::endl;

		// Get file intent
		unsigned intent;
		status = H5Fget_intent(file_id, &intent);
		std::cerr << "DEBUG: H5Fget_intent returned status=" << status << ", intent=" << intent << std::endl;

		// Try to open first dataset
		std::cerr << "DEBUG: Attempting H5Dopen2 for: " << SAMPLE_INDICES << std::endl;
		ds_indices_id = H5Dopen2(file_id, SAMPLE_INDICES, H5P_DEFAULT);
		std::cerr << "DEBUG: H5Dopen2 returned: " << ds_indices_id << std::endl;

		if (ds_indices_id < 0) {
			// Get more detailed error info
			std::cerr << "ERROR: H5Dopen2 failed for " << SAMPLE_INDICES << std::endl;
			std::cerr << "ERROR: Checking if path exists..." << std::endl;

			// Check if the dataset path exists by trying to open as object
			htri_t exists = H5Lexists(file_id, SAMPLE_INDICES, H5P_DEFAULT);
			std::cerr << "ERROR: H5Lexists returned: " << exists << " (1=exists, 0=no, <0=error)" << std::endl;

			if (exists > 0) {
				// Path exists but can't open as dataset - check what type it is
				H5O_info2_t obj_info;
				status = H5Oget_info_by_name3(file_id, SAMPLE_INDICES, &obj_info, H5O_INFO_BASIC, H5P_DEFAULT);
				std::cerr << "ERROR: H5Oget_info_by_name3 status=" << status << std::endl;
				if (status >= 0) {
					std::cerr << "ERROR: Object type: " << (int)obj_info.type
					          << " (0=group, 1=dataset, 2=datatype)" << std::endl;
				}
			}

			H5Fclose(file_id);
			throw std::runtime_error("Failed to open dataset " + std::string(SAMPLE_INDICES));
		}

		std::cerr << "DEBUG: Opening remaining datasets..." << std::endl;
		ds_indptr_id = H5Dopen2(file_id, SAMPLE_INDPTR, H5P_DEFAULT);
		if (ds_indptr_id < 0) throw std::runtime_error("Failed to open " + std::string(SAMPLE_INDPTR));

		ds_data_id = H5Dopen2(file_id, SAMPLE_DATA, H5P_DEFAULT);
		if (ds_data_id < 0) throw std::runtime_error("Failed to open " + std::string(SAMPLE_DATA));

		ds_samp_ids_id = H5Dopen2(file_id, SAMPLE_IDS, H5P_DEFAULT);
		if (ds_samp_ids_id < 0) throw std::runtime_error("Failed to open " + std::string(SAMPLE_IDS));

		ds_obs_ids_id = H5Dopen2(file_id, OBS_IDS, H5P_DEFAULT);
		if (ds_obs_ids_id < 0) throw std::runtime_error("Failed to open " + std::string(OBS_IDS));

		std::cerr << "DEBUG: All C API calls succeeded, wrapping in C++ objects..." << std::endl;

		// Wrap C handles in C++ objects
		file_handle = H5::H5File();
		file_handle.setId(file_id);

		ds_indices = H5::DataSet();
		ds_indices.setId(ds_indices_id);

		ds_indptr = H5::DataSet();
		ds_indptr.setId(ds_indptr_id);

		ds_data = H5::DataSet();
		ds_data.setId(ds_data_id);

		ds_samp_ids = H5::DataSet();
		ds_samp_ids.setId(ds_samp_ids_id);

		ds_obs_ids = H5::DataSet();
		ds_obs_ids.setId(ds_obs_ids_id);

		std::cerr << "DEBUG: BIOMReader constructor completed successfully" << std::endl;

	} catch (const std::exception &e) {
		std::cerr << "ERROR: Exception in BIOMReader constructor: " << e.what() << std::endl;

		// Clean up any open handles
		if (ds_obs_ids_id >= 0) H5Dclose(ds_obs_ids_id);
		if (ds_samp_ids_id >= 0) H5Dclose(ds_samp_ids_id);
		if (ds_data_id >= 0) H5Dclose(ds_data_id);
		if (ds_indptr_id >= 0) H5Dclose(ds_indptr_id);
		if (ds_indices_id >= 0) H5Dclose(ds_indices_id);
		if (file_id >= 0) H5Fclose(file_id);

		throw;
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
		// Open file without file locking property list
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
