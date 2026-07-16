#pragma once

#include <string>

namespace miint {

// The miint on-disk cache directory for a given subdirectory. For example
// MiintCacheDir("taxonomy") returns "${XDG_CACHE_HOME:-$HOME/.cache}/miint/taxonomy".
// Returns an empty string if neither XDG_CACHE_HOME nor HOME is set.
std::string MiintCacheDir(const std::string &subdir);

} // namespace miint
