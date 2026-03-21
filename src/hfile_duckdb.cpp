#include "hfile_duckdb.hpp"

#ifdef MIINT_STATIC_BUILD

#include <cerrno>
#include <memory>
#include <htslib-1.22.1/hfile_internal.h>
#include "duckdb/common/file_system.hpp"

using duckdb::FileHandle;
using duckdb::FileOpenFlags;
using duckdb::FileSystem;
using std::string;
using std::unique_ptr;

namespace miint {

struct hFILE_duckdb {
	hFILE base; // must be first member
	unique_ptr<FileHandle> handle;
};

static ssize_t duckdb_hfile_read(hFILE *fpv, void *buffer, size_t nbytes) {
	auto *fp = reinterpret_cast<hFILE_duckdb *>(fpv);
	try {
		auto n = fp->handle->Read(buffer, nbytes);
		return static_cast<ssize_t>(n);
	} catch (...) {
		errno = EIO;
		return -1;
	}
}

static off_t duckdb_hfile_seek(hFILE *fpv, off_t offset, int whence) {
	auto *fp = reinterpret_cast<hFILE_duckdb *>(fpv);
	try {
		if (whence == SEEK_SET) {
			fp->handle->Seek(offset);
		} else if (whence == SEEK_CUR) {
			auto cur = fp->handle->SeekPosition();
			fp->handle->Seek(cur + offset);
		} else if (whence == SEEK_END) {
			auto size = fp->handle->GetFileSize();
			fp->handle->Seek(size + offset);
		} else {
			errno = EINVAL;
			return -1;
		}
		return static_cast<off_t>(fp->handle->SeekPosition());
	} catch (...) {
		// Seek not supported (e.g., some streaming HTTP connections)
		errno = ESPIPE;
		return -1;
	}
}

static int duckdb_hfile_close(hFILE *fpv) {
	auto *fp = reinterpret_cast<hFILE_duckdb *>(fpv);
	// Release the FileHandle — unique_ptr destructor calls Close()
	fp->handle.reset();
	return 0;
}

static const struct hFILE_backend duckdb_backend = {duckdb_hfile_read, nullptr, duckdb_hfile_seek, nullptr,
                                                    duckdb_hfile_close};

hFILE *hfile_duckdb_open(FileSystem &fs, const string &path) {
	try {
		auto handle = fs.OpenFile(path, FileOpenFlags::FILE_FLAGS_READ);
		if (!handle) {
			errno = ENOENT;
			return nullptr;
		}

		auto *fp = reinterpret_cast<hFILE_duckdb *>(hfile_init(sizeof(hFILE_duckdb), "r", 0));
		if (!fp) {
			return nullptr;
		}

		// Placement-new the unique_ptr into the already-allocated struct
		new (&fp->handle) unique_ptr<FileHandle>(std::move(handle));
		fp->base.backend = &duckdb_backend;
		return &fp->base;
	} catch (...) {
		errno = EIO;
		return nullptr;
	}
}

} // namespace miint

#endif
