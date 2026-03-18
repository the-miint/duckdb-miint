#pragma once

#include <climits>
#include <cstdint>
#include <cstdlib>
#include <cerrno>
#include <string>

namespace miint {

// Safe numeric conversion: returns false on failure instead of throwing.
// Used by MzMLReader and MzXMLReader for parsing XML attribute values.

[[nodiscard]] inline bool safe_stoi(const std::string &s, int32_t &out) {
	if (s.empty()) {
		return false;
	}
	try {
		size_t pos = 0;
		long val = std::stol(s, &pos);
		if (pos != s.size() || val < INT32_MIN || val > INT32_MAX) {
			return false;
		}
		out = static_cast<int32_t>(val);
		return true;
	} catch (...) {
		return false;
	}
}

[[nodiscard]] inline bool safe_stod(const std::string &s, double &out) {
	if (s.empty()) {
		return false;
	}
	try {
		size_t pos = 0;
		out = std::stod(s, &pos);
		return pos == s.size();
	} catch (...) {
		return false;
	}
}

// C-string overloads use strtol/strtod directly to avoid std::string allocation
[[nodiscard]] inline bool safe_stoi(const char *s, int32_t &out) {
	if (!s || !*s) {
		return false;
	}
	char *end = nullptr;
	errno = 0;
	long val = std::strtol(s, &end, 10);
	if (end == s || *end != '\0' || errno == ERANGE || val < INT32_MIN || val > INT32_MAX) {
		return false;
	}
	out = static_cast<int32_t>(val);
	return true;
}

[[nodiscard]] inline bool safe_stod(const char *s, double &out) {
	if (!s || !*s) {
		return false;
	}
	char *end = nullptr;
	errno = 0;
	double val = std::strtod(s, &end);
	if (end == s || *end != '\0' || errno == ERANGE) {
		return false;
	}
	out = val;
	return true;
}

} // namespace miint
