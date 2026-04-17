#include "remote_file_helper.hpp"

#include <chrono>
#include <cstdio>
#include <cstring>
#include <random>
#include <thread>
#ifdef _WIN32
#include <windows.h>
#endif

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/storage/buffer_manager.hpp"
#endif

namespace miint {

// ===== ResolvedFileSet =====

ResolvedFileSet::~ResolvedFileSet() {
	for (auto &f : files_) {
		if (f.is_temp) {
			try {
				std::remove(f.local_path.c_str());
			} catch (...) {
				// Exception-safe: best-effort cleanup
			}
		}
	}
}

void ResolvedFileSet::Add(ResolvedFile file) {
	files_.push_back(std::move(file));
}

const std::vector<ResolvedFile> &ResolvedFileSet::Files() const {
	return files_;
}

// ===== RemoteFileHelper =====

static bool HasSchemePrefix(const std::string &path, const char *scheme) {
	size_t len = std::strlen(scheme);
	return path.size() > len && path.compare(0, len, scheme) == 0;
}

bool RemoteFileHelper::IsRemotePath(const std::string &path) {
	return HasSchemePrefix(path, "https://") || HasSchemePrefix(path, "http://") || HasSchemePrefix(path, "s3://") ||
	       HasSchemePrefix(path, "s3a://") || HasSchemePrefix(path, "s3n://") || HasSchemePrefix(path, "gcs://") ||
	       HasSchemePrefix(path, "gs://") || HasSchemePrefix(path, "r2://") || HasSchemePrefix(path, "hf://");
}

ResolvedFile RemoteFileHelper::ResolveToLocal(const std::string &path) {
	if (IsRemotePath(path)) {
		throw std::runtime_error(
		    "Remote path requires FileSystem context. Use the overload with FileSystem & ClientContext.");
	}
	return {path, false};
}

#ifdef MIINT_STATIC_BUILD

static constexpr idx_t DOWNLOAD_BUFFER_SIZE = 1048576; // 1MB — matches DuckDB httpfs chunk size

std::string RemoteFileHelper::GetTempDirectory(duckdb::ClientContext &context) {
	auto &buffer_manager = duckdb::BufferManager::GetBufferManager(context);
	const auto &temp_dir = buffer_manager.GetTemporaryDirectory();
	if (!temp_dir.empty()) {
		return temp_dir;
	}
#ifdef EMSCRIPTEN
	return "/tmp";
#elif defined(_WIN32)
	char buf[260];
	DWORD len = GetTempPathA(260, buf);
	if (len > 0 && len < 260) {
		return std::string(buf, len);
	}
	return "C:\\Temp";
#else
	return "/tmp";
#endif
}

static std::string GenerateTempFileName() {
	// Seed from real entropy to avoid collisions across concurrent threads
	auto seed = static_cast<uint64_t>(std::chrono::steady_clock::now().time_since_epoch().count()) ^
	            (std::hash<std::thread::id> {}(std::this_thread::get_id()) * 2654435761ULL);
	std::mt19937_64 rng(seed);
	auto rand_val = rng();
	return duckdb::StringUtil::Format("miint_%016llx.tmp", static_cast<unsigned long long>(rand_val));
}

ResolvedFile RemoteFileHelper::ResolveToLocal(duckdb::FileSystem &fs, duckdb::ClientContext &context,
                                              const std::string &path) {
	if (!IsRemotePath(path)) {
		return {path, false};
	}

	using duckdb::FileOpenFlags;
	auto source = fs.OpenFile(path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));

	auto temp_dir = GetTempDirectory(context);
	// DuckDB creates the temp directory lazily when the buffer manager spills;
	// a pure-read workload may never trigger that path, so ensure it exists
	// before we try to write. Recursive form handles user-configured nested
	// paths (e.g. SET temp_directory = '/var/run/myapp/nested/tmp').
	fs.CreateDirectoriesRecursive(temp_dir);
	auto temp_path = fs.JoinPath(temp_dir, GenerateTempFileName());

	try {
		auto dest = fs.OpenFile(
		    temp_path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_WRITE | FileOpenFlags::FILE_FLAGS_FILE_CREATE_NEW));

		auto buffer = duckdb::unique_ptr<char[]>(new char[DOWNLOAD_BUFFER_SIZE]);
		while (true) {
			auto bytes_read = source->Read(buffer.get(), DOWNLOAD_BUFFER_SIZE);
			if (bytes_read == 0) {
				break;
			}
			dest->Write(buffer.get(), bytes_read);
		}
		dest->Close();
	} catch (...) {
		std::remove(temp_path.c_str());
		throw;
	}

	return {temp_path, true};
}

ResolvedFileSet RemoteFileHelper::ResolveAllToLocal(duckdb::FileSystem &fs, duckdb::ClientContext &context,
                                                    const std::vector<std::string> &paths) {
	ResolvedFileSet result;
	for (const auto &path : paths) {
		result.Add(ResolveToLocal(fs, context, path));
	}
	return result;
}

#endif // MIINT_STATIC_BUILD

} // namespace miint
