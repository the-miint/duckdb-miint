// SPDX-License-Identifier: MIT
//
// Portable env-var manipulation for tests. MinGW's libc doesn't declare
// POSIX setenv/unsetenv, so test files that touch env vars must route
// through SetEnvPortable (or the EnvGuard RAII wrapper) instead of calling
// ::setenv / ::unsetenv directly — otherwise the windows_amd64_mingw build
// breaks at compile time.

#pragma once

#include <cstdlib>
#include <string>
#include <utility>

namespace miint_test {

// _putenv_s with "" unsets on Windows.
inline void SetEnvPortable(const char *name, const char *value) {
#ifdef _WIN32
	_putenv_s(name, value ? value : "");
#else
	if (value) {
		::setenv(name, value, 1);
	} else {
		::unsetenv(name);
	}
#endif
}

// RAII: set on construct, restore prior value (or unset) on destruct.
class EnvGuard {
public:
	EnvGuard(std::string name, const char *value) : name_(std::move(name)) {
		const char *prior = std::getenv(name_.c_str());
		if (prior) {
			prior_ = prior;
			had_prior_ = true;
		}
		SetEnvPortable(name_.c_str(), value);
	}
	~EnvGuard() {
		SetEnvPortable(name_.c_str(), had_prior_ ? prior_.c_str() : nullptr);
	}

	EnvGuard(const EnvGuard &) = delete;
	EnvGuard &operator=(const EnvGuard &) = delete;

private:
	std::string name_;
	std::string prior_;
	bool had_prior_ {false};
};

} // namespace miint_test
