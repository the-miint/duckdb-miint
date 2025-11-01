#pragma once
#include <optional>
#include <vector>
#include <string>
#include <memory>
#include "BIOMTable.hpp"
#include <hdf5.h>

namespace miint {

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
