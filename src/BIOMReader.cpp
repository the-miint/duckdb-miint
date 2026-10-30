#include "BIOMReader.hpp"
#include <sys/stat.h>
#include <cerrno>
#include <cstring>

namespace miint {

// Storage for the process-wide HDF5 mutex declared in BIOMReader.hpp.
// See the header for the contract and rationale.
std::mutex g_hdf5_mutex;

static herr_t KeepInnermostHdf5Error(unsigned n, const H5E_error2_t *err, void *desc) {
	if (n == 0) {
		*static_cast<const char **>(desc) = err->desc;
	}
	return 0;
}

// The innermost entry of an HDF5 error stack, cut to the OS error when the entry quotes one (a failed open's
// "error message = '...'"); empty for an empty stack
static std::string Hdf5ErrorCause(hid_t error_stack) {
	const char *innermost = nullptr;
	H5Ewalk2(error_stack, H5E_WALK_UPWARD, KeepInnermostHdf5Error, &innermost);
	if (!innermost) {
		return {};
	}
	std::string cause = innermost;
	static const std::string os_error_tag = "error message = '";
	const auto tag = cause.find(os_error_tag);
	if (tag != std::string::npos) {
		const auto start = tag + os_error_tag.size();
		const auto end = cause.find('\'', start);
		if (end != std::string::npos) {
			cause = cause.substr(start, end - start);
		}
	}
	return cause;
}

Hdf5FailureRecorder::Hdf5FailureRecorder() : hdf5_guard(g_hdf5_mutex) {
	H5Eget_auto2(H5E_DEFAULT, &saved_func, &saved_data);
	H5Eset_auto2(H5E_DEFAULT, Record, this);
}

Hdf5FailureRecorder::~Hdf5FailureRecorder() {
	H5Eset_auto2(H5E_DEFAULT, saved_func, saved_data);
}

// HDF5 calls this as a failing call returns, with that call's error stack. HDF5 is C, so no exception may
// leave it: a cause that cannot be kept is left empty.
herr_t Hdf5FailureRecorder::Record(hid_t error_stack, void *recorder) {
	auto &self = *static_cast<Hdf5FailureRecorder *>(recorder);
	if (!self.failed) {
		self.failed = true;
		try {
			self.cause = Hdf5ErrorCause(error_stack);
		} catch (...) { // NOLINT: must not unwind through HDF5
			self.cause.clear();
		}
	}
	return 0;
}

BIOMReader::BIOMReader(const std::string &path) {
	// Check if file exists and get size
	struct stat st;
	if (stat(path.c_str(), &st) != 0) {
		const int stat_errno = errno;
		throw std::runtime_error(std::string("cannot stat the file: ") + std::strerror(stat_errno));
	}

	htri_t is_hdf5 = -1;
	hid_t file_id = -1;
	hid_t ds_indices_id = -1, ds_indptr_id = -1, ds_data_id = -1;
	hid_t ds_samp_ids_id = -1, ds_obs_ids_id = -1;

	try {
		// Check if file is actually HDF5
		is_hdf5 = H5Fis_hdf5(path.c_str());
		if (is_hdf5 < 0) {
			throw std::runtime_error("H5Fis_hdf5 failed");
		}
		if (is_hdf5 == 0) {
			throw std::runtime_error("not an HDF5 file");
		}

		// Open file with C API
		file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

		if (file_id < 0) {
			throw std::runtime_error("H5Fopen failed");
		}

		ds_indices_id = H5Dopen2(file_id, SAMPLE_INDICES, H5P_DEFAULT);
		if (ds_indices_id < 0) {
			throw std::runtime_error("Failed to open dataset " + std::string(SAMPLE_INDICES));
		}

		ds_indptr_id = H5Dopen2(file_id, SAMPLE_INDPTR, H5P_DEFAULT);
		if (ds_indptr_id < 0) {
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_INDPTR));
		}

		ds_data_id = H5Dopen2(file_id, SAMPLE_DATA, H5P_DEFAULT);
		if (ds_data_id < 0) {
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_DATA));
		}

		ds_samp_ids_id = H5Dopen2(file_id, SAMPLE_IDS, H5P_DEFAULT);
		if (ds_samp_ids_id < 0) {
			throw std::runtime_error("Failed to open " + std::string(SAMPLE_IDS));
		}

		ds_obs_ids_id = H5Dopen2(file_id, OBS_IDS, H5P_DEFAULT);
		if (ds_obs_ids_id < 0) {
			throw std::runtime_error("Failed to open " + std::string(OBS_IDS));
		}

		// Store the C handles directly
		file_handle = file_id;
		ds_indices = ds_indices_id;
		ds_indptr = ds_indptr_id;
		ds_data = ds_data_id;
		ds_samp_ids = ds_samp_ids_id;
		ds_obs_ids = ds_obs_ids_id;

	} catch (const std::exception &e) {
		// Clean up any open handles
		if (ds_obs_ids_id >= 0) {
			H5Dclose(ds_obs_ids_id);
		}
		if (ds_samp_ids_id >= 0) {
			H5Dclose(ds_samp_ids_id);
		}
		if (ds_data_id >= 0) {
			H5Dclose(ds_data_id);
		}
		if (ds_indptr_id >= 0) {
			H5Dclose(ds_indptr_id);
		}
		if (ds_indices_id >= 0) {
			H5Dclose(ds_indices_id);
		}
		if (file_id >= 0) {
			H5Fclose(file_id);
		}
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
	const char *target = "format-version";

	htri_t is_hdf5 = -1;
	hid_t file_id = -1;
	hid_t attr_id = -1;

	struct stat st;
	if (stat(path.c_str(), &st) != 0) {
		return false;
	}

	// A file HDF5 cannot open is an error, not a file that is not BIOM, so its reason reaches the user
	is_hdf5 = H5Fis_hdf5(path.c_str());
	if (is_hdf5 < 0) {
		throw std::runtime_error("H5Fis_hdf5 failed");
	}
	if (is_hdf5 == 0) {
		return false;
	}

	file_id = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if (file_id < 0) {
		throw std::runtime_error("H5Fopen failed");
	}

	// From here every HDF5 failure is the file's, except a format-version that is not two integers, which
	// makes the file not BIOM. Looking the attribute up by name, rather than walking the root's attributes,
	// leaves no other failure for a Hdf5FailureRecorder to mistake for the cause of a later one.
	bool valid = false;
	const char *failed_call = nullptr;
	htri_t has_version = H5Aexists_by_name(file_id, "/", target, H5P_DEFAULT);
	if (has_version < 0) {
		failed_call = "H5Aexists_by_name failed";
	} else if (has_version > 0) {
		attr_id = H5Aopen_by_name(file_id, "/", target, H5P_DEFAULT, H5P_DEFAULT);
		if (attr_id < 0) {
			failed_call = "H5Aopen_by_name failed";
		} else {
			hid_t space_id = H5Aget_space(attr_id);
			if (space_id < 0) {
				failed_call = "H5Aget_space failed";
			} else {
				hssize_t n_values = H5Sget_simple_extent_npoints(space_id);
				H5Sclose(space_id);
				// Read only an attribute of exactly two values: H5Aread writes all of them into `values`
				int values[2];
				valid = n_values == 2 && H5Aread(attr_id, H5T_NATIVE_INT, values) >= 0 && values[0] == 2;
			}
			H5Aclose(attr_id);
		}
	}

	// Clean up
	H5Fclose(file_id);
	if (failed_call) {
		throw std::runtime_error(failed_call);
	}

	return valid;
}
} // namespace miint
