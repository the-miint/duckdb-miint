#pragma once
#include <optional>
#include <vector>
#include <string>
#include <memory>
#include <mutex>
#include "BIOMTable.hpp"
#include <hdf5.h>

namespace miint {

//! Process-wide mutex serializing every HDF5 call from miint.
//!
//! HDF5 (vcpkg build: features cpp/zlib/szip, no `threadsafe`) has global
//! library state — initialization, type system, error stack, dataset caches
//! — that races under concurrent access. Two read_biom() table-function
//! instances each have their own GlobalState, so a per-state mutex does
//! NOT serialize cross-instance calls. Mirrors `g_mafft_mutex` in
//! MafftAligner.cpp.
//!
//! Acquire around every HDF5 entry point: H5Fopen, H5Dopen2, H5Dread,
//! H5Fclose, H5Dclose, plus the BIOMReader / BIOMTable scope that drives
//! them — including the BIOMReader destructor so file/dataset handles are
//! closed under the lock.
extern std::mutex g_hdf5_mutex;


/* datasets defined by the BIOM 2.x spec */
static constexpr const char *OBS_INDPTR = "/observation/matrix/indptr";
static constexpr const char *OBS_INDICES = "/observation/matrix/indices";
static constexpr const char *OBS_DATA = "/observation/matrix/data";
static constexpr const char *OBS_IDS = "/observation/ids";

static constexpr const char *SAMPLE_INDPTR = "/sample/matrix/indptr";
static constexpr const char *SAMPLE_INDICES = "/sample/matrix/indices";
static constexpr const char *SAMPLE_DATA = "/sample/matrix/data";
static constexpr const char *SAMPLE_IDS = "/sample/ids";

class BIOMReader {
private:
	hid_t file_handle;
	hid_t ds_indices;
	hid_t ds_indptr;
	hid_t ds_data;
	hid_t ds_samp_ids;
	hid_t ds_obs_ids;

public:
	explicit BIOMReader(const std::string &path1);
	~BIOMReader();
	BIOMTable read() const;
	static bool IsBIOM(const std::string &path);
};
} // namespace miint
