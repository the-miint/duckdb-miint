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
//! file/dataset handles are closed under the lock. read_biom and COPY FORMAT
//! BIOM hold it through a Hdf5FailureRecorder for all their HDF5 work, the
//! closing destructors included.
//!
//! HDF5 tears down from atexit handlers outside this lock, so a process must
//! not exit while any HDF5 call is in progress.
extern std::mutex g_hdf5_mutex;

//! While alive, keeps the cause of the first HDF5 call that fails, taken as that call fails: the innermost
//! entry of its error stack, cut to the OS error when the entry quotes one ("File exists", for example).
//! Every HDF5 call clears the error stack, so a cause read later is lost to the calls made in between, such
//! as the closes in the destructors an exception unwinds through. Use it around code that stops at its first
//! failure: a failure the code does not stop at is kept as the cause of a later one. The stack and the hook
//! this sets are the process's (this HDF5 is not thread-safe), so the recorder holds g_hdf5_mutex for its
//! whole life, and code it covers must not take that mutex again.
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
	//! `what` failed, and why: "<what>: <cause>", or `what` alone without a cause
	std::string Explain(const std::string &what) const {
		return cause.empty() ? what : what + ": " + cause;
	}

private:
	static herr_t Record(hid_t error_stack, void *recorder);

	// First, so the mutex is taken before the hook is set and released after it is restored
	std::lock_guard<std::mutex> hdf5_guard;
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
	//! False for a missing file, a file that is not HDF5, and an HDF5 file whose format-version is not two
	//! values HDF5 reads as integers, the first 2. Throws std::runtime_error naming the HDF5 call that failed
	//! to open the file or to find its format-version; a Hdf5FailureRecorder covering the call has the reason.
	static bool IsBIOM(const std::string &path);
};
} // namespace miint
