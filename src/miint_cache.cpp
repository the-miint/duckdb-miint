#include "miint_cache.hpp"

#include <cstdlib>

namespace miint {

std::string MiintCacheDir(const std::string &subdir) {
	if (const char *xdg = ::getenv("XDG_CACHE_HOME"); xdg && xdg[0] != '\0') {
		return std::string(xdg) + "/miint/" + subdir;
	}
	if (const char *home = ::getenv("HOME"); home && home[0] != '\0') {
		return std::string(home) + "/.cache/miint/" + subdir;
	}
	return {};
}

} // namespace miint
