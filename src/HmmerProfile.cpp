#include "HmmerProfile.hpp"

#include <cstdarg>
#include <mutex>

extern "C" {
#include "easel.h"
}

namespace miint {
namespace hmmer {

namespace {
void IgnoreException(int, int, char *, int, char *, va_list) {
}
} // namespace

void Initialize() {
	static std::once_flag initialized;
	std::call_once(initialized, [] { esl_exception_SetHandler(IgnoreException); });
}

std::string TakeExceptionMessage() {
	return "";
}

ProfileHeader ParseProfile(std::string_view) {
	throw std::logic_error("ParseProfile: not implemented");
}

namespace detail {
void ThrowIfReadFailed(int, bool, const char *, int64_t) {
	throw std::logic_error("ThrowIfReadFailed: not implemented");
}
} // namespace detail

} // namespace hmmer
} // namespace miint
