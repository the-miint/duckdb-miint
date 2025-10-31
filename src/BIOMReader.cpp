#include "BIOMReader.hpp"
#include <iostream>
#include <sys/stat.h>

namespace miint {

BIOMReader::BIOMReader(const std::string &path1) {
	// Enable detailed HDF5 error reporting
	H5::Exception::printErrorStack();

	try {
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

		// Try opening WITHOUT the file locking property list first
		std::cerr << "DEBUG: Attempting to open HDF5 file (without file locking property)..." << std::endl;
		file_handle = H5::H5File(path1, H5F_ACC_RDONLY);
		std::cerr << "DEBUG: File opened successfully, handle ID: " << file_handle.getId() << std::endl;

		try {
			std::cerr << "DEBUG: Opening dataset: " << SAMPLE_INDICES << std::endl;
			ds_indices = file_handle.openDataSet(SAMPLE_INDICES);
			std::cerr << "DEBUG: Opening dataset: " << SAMPLE_INDPTR << std::endl;
			ds_indptr = file_handle.openDataSet(SAMPLE_INDPTR);
			std::cerr << "DEBUG: Opening dataset: " << SAMPLE_DATA << std::endl;
			ds_data = file_handle.openDataSet(SAMPLE_DATA);
			std::cerr << "DEBUG: Opening dataset: " << SAMPLE_IDS << std::endl;
			ds_samp_ids = file_handle.openDataSet(SAMPLE_IDS);
			std::cerr << "DEBUG: Opening dataset: " << OBS_IDS << std::endl;
			ds_obs_ids = file_handle.openDataSet(OBS_IDS);
			std::cerr << "DEBUG: All datasets opened successfully" << std::endl;
		} catch (...) {
			std::cerr << "ERROR: Failed to open one or more datasets" << std::endl;
			if (file_handle.getId() >= 0) {
				file_handle.close();
			}
			throw;
		}
	} catch (H5::FileIException &e) {
		std::cerr << "ERROR: H5::FileIException - " << e.getDetailMsg() << std::endl;
		std::cerr << "ERROR: Function: " << e.getFuncName() << std::endl;
		throw std::runtime_error("Failed to open HDF5 file: " + std::string(e.getDetailMsg()));
	} catch (H5::DataSetIException &e) {
		std::cerr << "ERROR: H5::DataSetIException - " << e.getDetailMsg() << std::endl;
		std::cerr << "ERROR: Function: " << e.getFuncName() << std::endl;
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
