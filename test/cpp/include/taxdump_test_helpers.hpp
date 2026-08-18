#pragma once

#include <stdexcept>
#include <string>
#include <vector>

namespace miint_test {

// True if `fn` throws a std::runtime_error whose message contains every fragment.
//
// The taxdump layer reports malformed input as std::runtime_error (it is deliberately
// DuckDB-free, so it cannot use DuckDB's exception types). Its messages have to be
// actionable on their own -- a consumer reading one should learn which member and which
// line moved -- so the fragments are asserted, not merely the throw. Catch2's
// CHECK_THROWS_WITH matches the whole string, which would pin the exact wording; this
// pins only the facts that must be present.
template <typename Fn>
bool ThrowsMentioning(Fn &&fn, const std::vector<std::string> &fragments) {
	try {
		fn();
	} catch (const std::runtime_error &e) {
		std::string msg = e.what();
		for (const auto &frag : fragments) {
			if (msg.find(frag) == std::string::npos) {
				return false;
			}
		}
		return true;
	} catch (...) {
		return false;
	}
	return false;
}

} // namespace miint_test
