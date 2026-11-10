#pragma once

#include <string_view>

namespace miint {

// Whether `path` is an http(s) URL, rather than a local path or another scheme (s3://,
// ...). Apart from DuckDB, so the DuckDB-free file rules can use it too.
inline bool IsHttpUrl(std::string_view path) {
	return path.substr(0, 7) == "http://" || path.substr(0, 8) == "https://";
}

} // namespace miint
