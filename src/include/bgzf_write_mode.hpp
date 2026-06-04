#pragma once

// Open-mode selection for the local bgzf gzip writer (COPY ... COMPRESSION gzip).
//
// Returns the htslib bgzf_open mode string. Regular-file (or not-yet-existing) targets get "wx6"
// -- write, exclusive-create (O_EXCL), gzip level 6 -- which preserves the FILE_CREATE_NEW "fail if
// the file already exists" contract atomically (no TOCTOU). Targets that ALWAYS already exist --
// stdout/stderr, /dev/fd aliases, and pre-created pipes/devices -- must NOT get O_EXCL, or the open
// fails with EEXIST and breaks gzip streaming to a pipe or /dev/stdout (a regression vs the
// uncompressed BufferedFileWriter path, which streams to them fine); those get plain "w6".
//
// Kept DuckDB-free (only <string> + <sys/stat.h>) so it is unit-testable in the standalone catch2
// test binary, which deliberately does not link libduckdb.

#include <string>
#include <sys/stat.h>

namespace miint {

inline const char *BgzfWriteMode(const std::string &path) {
	auto starts_with = [&](const char *prefix) {
		return path.rfind(prefix, 0) == 0;
	};
	// Named stdout/stderr/fd aliases: never exclusive. Matched by name, not by stat target --
	// stdout may be redirected to a regular file, which would also trip O_EXCL.
	if (path == "-" || path == "/dev/stdout" || path == "/dev/stderr" || starts_with("/dev/fd/") ||
	    starts_with("/proc/self/fd/")) {
		return "w6";
	}
	// An existing non-regular file (FIFO / character / block / socket) is a deliberate streaming
	// sink; only a regular file (or a not-yet-existing path) keeps the no-clobber/no-TOCTOU guard.
	struct stat st;
	if (stat(path.c_str(), &st) == 0 && !S_ISREG(st.st_mode)) {
		return "w6";
	}
	return "wx6";
}

} // namespace miint
