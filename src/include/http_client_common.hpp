#pragma once

#include "interruptible_sleep.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/main/client_context.hpp"

#include <chrono>
#include <cstddef>
#include <functional>
#include <string>

namespace duckdb {
struct HTTPResponse;
} // namespace duckdb

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

// Whether `path` is an http(s) URL, which the clients fetch through DuckDB's HTTP client.
// Other remote schemes (s3://, ...) go through DuckDB's file system instead.
inline bool IsHttpUrl(const std::string &path) {
	return duckdb::StringUtil::StartsWith(path, "http://") || duckdb::StringUtil::StartsWith(path, "https://");
}

// Sleeps `duration`, unless Ctrl-C or max_execution_time ends the query first. Every
// wait in the HTTP clients goes through here. Not DuckDB's ThreadUtil::SleepMs, which
// slices the same way but checks only IsInterrupted, which does not see the deadline.
inline void SleepUnlessInterrupted(duckdb::ClientContext &context, std::chrono::milliseconds duration) {
	InterruptibleSleep(duration, [&context] { context.BlockingInterruptCheck(); });
}

// GETs `url` and hands its body to `write` as it arrives, rather than holding it, so a
// file of any size can be saved, or only its start read: `write` returns false to stop
// the transfer once it has what it needs. `begin` is called, with the response, before
// each attempt's body, which starts over from its first byte: whatever an earlier
// attempt delivered is to be dropped. Ctrl-C and max_execution_time stop the transfer.
//
// Transient failures (5xx, 429, 408) are tried again after 1, 2 and 4 s, waits that
// Ctrl-C and max_execution_time end. DuckDB's own retries (http_retries) run inside each
// attempt. A failure throws IOException, or HTTPException with the status, naming
// `fn_name` and the URL; an exception from `begin` or `write` is rethrown as it is.
void HttpGetStreaming(duckdb::ClientContext &context, const std::string &url, const char *fn_name,
                      const std::function<void(const duckdb::HTTPResponse &)> &begin,
                      const std::function<bool(const char *, size_t)> &write);

} // namespace miint
