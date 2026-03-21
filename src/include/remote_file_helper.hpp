#pragma once
#include <string>
#include <vector>

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/file_system.hpp"
#include "duckdb/main/client_context.hpp"
#endif

namespace miint {

// Result of resolving a path that may be remote
struct ResolvedFile {
	std::string local_path;
	bool is_temp;

	ResolvedFile(std::string path, bool temp) : local_path(std::move(path)), is_temp(temp) {
	}
};

// RAII wrapper — deletes temp files on destruction
class ResolvedFileSet {
public:
	~ResolvedFileSet();

	// Non-copyable (destructor deletes temp files; copy would double-delete)
	ResolvedFileSet(const ResolvedFileSet &) = delete;
	ResolvedFileSet &operator=(const ResolvedFileSet &) = delete;
	ResolvedFileSet(ResolvedFileSet &&) = default;
	ResolvedFileSet &operator=(ResolvedFileSet &&) = default;
	ResolvedFileSet() = default;

	void Add(ResolvedFile file);
	const std::vector<ResolvedFile> &Files() const;

private:
	std::vector<ResolvedFile> files_;
};

class RemoteFileHelper {
public:
	// Detect remote paths (https://, http://, s3://, etc.)
	static bool IsRemotePath(const std::string &path);

	// For local files only — returns path as-is with is_temp=false.
	// Use the overload with FileSystem for remote paths.
	static ResolvedFile ResolveToLocal(const std::string &path);

#ifdef MIINT_STATIC_BUILD
	// Download remote file to temp, return local path. Local paths returned as-is.
	static ResolvedFile ResolveToLocal(duckdb::FileSystem &fs, duckdb::ClientContext &context, const std::string &path);

	// Resolve multiple paths
	static ResolvedFileSet ResolveAllToLocal(duckdb::FileSystem &fs, duckdb::ClientContext &context,
	                                         const std::vector<std::string> &paths);

	// Get cross-platform temp directory
	static std::string GetTempDirectory(duckdb::ClientContext &context);
#endif
};

} // namespace miint
