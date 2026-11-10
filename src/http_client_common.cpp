#include "http_client_common.hpp"

#include "ensure_httpfs.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/main/database.hpp"

#include <exception>

namespace miint {

void HttpGetStreaming(duckdb::ClientContext &context, const std::string &url, const char *fn_name,
                      const std::function<void(const duckdb::HTTPResponse &)> &begin,
                      const std::function<void(const char *, size_t)> &write) {
	EnsureHttpfsLoaded(context);
	duckdb::DatabaseInstance &db = *context.db;
	duckdb::HTTPHeaders headers(db);
	auto &http_util = duckdb::HTTPUtil::Get(db);
	auto params = http_util.InitializeParameters(db, url);

	// How httpfs's clients call the handlers: docs/internals/duckdb-engine-notes.md, "A GET
	// that streams its body". An exception must not unwind through a client, which may not
	// catch it, so it is held and rethrown once the request returns. Only a body that
	// `begin` has started reaches `write`. httpfs's clients never hand an error page to the
	// content handler; duckdb-wasm's hands it every body, error pages included, and never
	// calls the response handler.
	std::exception_ptr error;
	bool body_started = false;
	duckdb::GetRequestInfo get_request(
	    url, headers, *params,
	    // Called once per attempt with its response's status: a 2xx's before its body,
	    // which then starts over.
	    [&](const duckdb::HTTPResponse &response) {
		    body_started = false;
		    if (error) {
			    return false;
		    }
		    const auto status = static_cast<int>(response.status);
		    if (status < 200 || status >= 300) {
			    return true;
		    }
		    try {
			    begin(response);
		    } catch (...) {
			    error = std::current_exception();
			    return false;
		    }
		    body_started = true;
		    return true;
	    },
	    [&](duckdb::const_data_ptr_t data, duckdb::idx_t size) {
		    if (error) {
			    return false;
		    }
		    if (!body_started) {
			    return true;
		    }
		    try {
			    context.BlockingInterruptCheck();
			    write(duckdb::const_char_ptr_cast(data), size);
		    } catch (...) {
			    error = std::current_exception();
			    return false;
		    }
		    return true;
	    });
	get_request.try_request = true;

	constexpr int MAX_RETRIES = 3;
	int retry_delay_ms = 1000;
	for (int attempt = 0; attempt <= MAX_RETRIES; ++attempt) {
		error = nullptr;
		auto response = http_util.Request(get_request);
		if (error) {
			std::rethrow_exception(error);
		}
		if (response->Success()) {
			// A client that never calls the response handler (duckdb-wasm's) leaves the body
			// in response->body, and what it handed the content handler was dropped above.
			if (!body_started && !response->body.empty()) {
				begin(*response);
				write(response->body.data(), response->body.size());
			}
			return;
		}
		const auto status = static_cast<int>(response->status);
		const bool retryable = !response->HasRequestError() && (status >= 500 || status == 429 || status == 408);
		if (attempt < MAX_RETRIES && retryable) {
			SleepUnlessInterrupted(context, std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2;
			continue;
		}
		if (response->HasRequestError()) {
			throw duckdb::IOException("%s: download failed: %s (URL: %s)", fn_name, response->GetRequestError(), url);
		}
		throw duckdb::HTTPException(*response, "%s: download failed with HTTP %d (URL: %s)", fn_name, status, url);
	}
	throw duckdb::IOException("%s: download failed after retries (URL: %s)", fn_name, url);
}

} // namespace miint
