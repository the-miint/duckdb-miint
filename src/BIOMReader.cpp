#include "BIOMReader.hpp"

namespace miint {

BIOMReader::BIOMReader(const std::string &path1) {
	try {
		file_handle = H5::H5File(path1, H5F_ACC_RDONLY);

		ds_indices = file_handle.openDataSet(SAMPLE_INDICES);
		ds_indptr = file_handle.openDataSet(SAMPLE_INDPTR);
		ds_data = file_handle.openDataSet(SAMPLE_DATA);
		ds_samp_ids = file_handle.openDataSet(SAMPLE_IDS);
		ds_obs_ids = file_handle.openDataSet(OBS_IDS);

		// load the datasets etc
	} catch (H5::FileIException &e) {
		throw std::runtime_error("Failed to open HDF5 file: " + std::string(e.getDetailMsg()));
	} catch (H5::DataSetIException &e) {
		throw std::runtime_error("Failed to access dataset: " + std::string(e.getDetailMsg()));
	}
}

BIOMTable BIOMReader::read() const {
	return BIOMTable(ds_indices, ds_indptr, ds_data, ds_samp_ids, ds_obs_ids);
}
} // namespace miint
