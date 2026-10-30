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
//! Acquire around every HDF5 entry point, including library-global
//! configuration at extension load (H5Eset_auto), plus the BIOMReader /
//! BIOMTable scope that drives them — including the BIOMReader destructor so
//! file/dataset handles are closed under the lock. COPY FORMAT BIOM
//! (copy_biom.cpp) holds it for its whole HDF5 write, the closing destructors
//! included.
//!
//! HDF5 tears down from atexit handlers outside this lock, so a process must
//! not exit while any HDF5 call is in progress.
extern std::mutex g_hdf5_mutex;

//! Why an HDF5 call failed: the innermost entry of an HDF5 error stack, cut to the OS error when the entry
//! quotes one (a failed open's "error message = '...'"); empty for an empty stack.
std::string Hdf5ErrorCause(hid_t error_stack = H5E_DEFAULT);

//! While alive, keeps the cause (Hdf5ErrorCause) of the first HDF5 call that fails, taken as that call fails.
//! Every HDF5 call clears the error stack, so a cause read later is lost to the calls made in between, such
//! as the closes in the destructors an exception unwinds through. The stack and the hook this sets are the
//! process's (this HDF5 is not thread-safe): hold g_hdf5_mutex for the recorder's whole life.
class Hdf5FailureRecorder {
public:
	Hdf5FailureRecorder();
	~Hdf5FailureRecorder();
	Hdf5FailureRecorder(const Hdf5FailureRecorder &) = delete;
	Hdf5FailureRecorder &operator=(const Hdf5FailureRecorder &) = delete;

	//! Empty if no call has failed, or the failure's stack named no cause
	const std::string &Cause() const {
		return cause;
	}

private:
	static herr_t Record(hid_t error_stack, void *recorder);

	H5E_auto2_t saved_func = nullptr;
	void *saved_data = nullptr;
	bool failed = false;
	std::string cause;
};

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
