#pragma once

#include "interruptible_sleep.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/client_context.hpp"

#include <chrono>
#include <string>

namespace miint {

// The URL a setting names, without a trailing slash, since the clients append paths to
// it. The setting is registered at load (src/miint_extension.cpp).
inline std::string UrlSetting(duckdb::ClientContext &context, const char *name) {
	duckdb::Value value;
	if (!context.TryGetCurrentSetting(name, value)) {
		throw duckdb::InternalException("miint extension option '%s' is not registered — was LoadInternal completed?",
		                                name);
	}
	std::string url = value.IsNull() ? std::string() : value.ToString();
	duckdb::StringUtil::RTrim(url, "/");
	if (url.empty()) {
		throw duckdb::InvalidInputException("%s cannot be empty; RESET it to use the default", name);
	}
	return url;
}

// Sleeps `duration`, unless Ctrl-C or max_execution_time ends the query first. Every
// wait in the HTTP clients goes through here. Not DuckDB's ThreadUtil::SleepMs, which
// slices the same way but checks only IsInterrupted, which does not see the deadline.
inline void SleepUnlessInterrupted(duckdb::ClientContext &context, std::chrono::milliseconds duration) {
	InterruptibleSleep(duration, [&context] { context.BlockingInterruptCheck(); });
}

} // namespace miint
