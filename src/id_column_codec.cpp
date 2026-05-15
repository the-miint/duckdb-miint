#include "id_column_codec.hpp"

#include <charconv>
#include <stdexcept>
#include <string>
#include <system_error>

namespace miint {

std::optional<int64_t> ParseIdAsInt64(const std::string &s) {
	if (s.empty() || s == "*") {
		return std::nullopt;
	}
	// std::from_chars is locale-independent and strict: it rejects leading
	// whitespace, '+' sign, and float-like input. Together with the trailing
	// check below this gives us exactly the round-trip contract we want.
	const char *begin = s.data();
	const char *end = begin + s.size();
	int64_t v = 0;
	auto [ptr, ec] = std::from_chars(begin, end, v);
	if (ec == std::errc::invalid_argument) {
		throw std::invalid_argument("ParseIdAsInt64: not a decimal integer: '" + s + "'");
	}
	if (ec == std::errc::result_out_of_range) {
		throw std::out_of_range("ParseIdAsInt64: out of int64 range: '" + s + "'");
	}
	if (ptr != end) {
		throw std::invalid_argument("ParseIdAsInt64: trailing characters in '" + s + "'");
	}
	return v;
}

std::string FormatIdFromInt64(int64_t v) {
	return std::to_string(v);
}

} // namespace miint
