// SPDX-License-Identifier: MIT
//
// See curl_send.hpp.

#include "curl_send.hpp"

#include <curl/curl.h>

#include <atomic>
#include <cstddef>
#include <cstring>
#include <exception>
#include <mutex>
#include <stdexcept>
#include <string>

namespace miint {

namespace {

// One-time global libcurl init. `curl_global_init` is the documented entry
// point; calling it more than once is benign but we gate it via std::call_once
// to keep the setup cost out of the hot path. Cleanup is intentionally not
// registered — process teardown handles it, and DuckDB extension unload paths
// are too varied to safely call curl_global_cleanup in.
void EnsureCurlGlobalInit() {
	static std::once_flag flag;
	std::call_once(flag, [] {
		CURLcode rc = curl_global_init(CURL_GLOBAL_DEFAULT);
		if (rc != CURLE_OK) {
			throw std::runtime_error(std::string("libcurl: curl_global_init failed: ") + curl_easy_strerror(rc));
		}
	});
}

// Thread-through state for the C-style read callback. libcurl invokes the
// callback synchronously from the same thread as curl_easy_perform; we use
// `captured` to defer producer-thrown exceptions until after curl returns
// control to us.
struct ReadContext {
	const BodyProducer *producer;
	std::exception_ptr captured;
};

std::size_t ReadCallback(char *buffer, std::size_t size, std::size_t nmemb, void *userdata) {
	auto *ctx = static_cast<ReadContext *>(userdata);
	const std::size_t want = size * nmemb;
	try {
		return (*ctx->producer)(buffer, want);
	} catch (...) {
		ctx->captured = std::current_exception();
		// CURL_READFUNC_ABORT signals libcurl to abort the transfer cleanly.
		// libcurl will return CURLE_ABORTED_BY_CALLBACK; we suppress that
		// error and rethrow the producer's exception in RunCurlUpload.
		return CURL_READFUNC_ABORT;
	}
}

// Capture the server's response body so HTTP 4xx/5xx + FTP failure messages
// surface in `error_message`. Capped at 4 KB — error bodies (XML / JSON /
// plain text) virtually never exceed this; the cap protects against a
// runaway server response.
struct WriteContext {
	std::string body;
	static constexpr std::size_t MAX_CAPTURED = 4096;
};

std::size_t WriteCallback(char *data, std::size_t size, std::size_t nmemb, void *userdata) {
	auto *ctx = static_cast<WriteContext *>(userdata);
	const std::size_t total = size * nmemb;
	if (ctx->body.size() < WriteContext::MAX_CAPTURED) {
		const std::size_t to_take = std::min(total, WriteContext::MAX_CAPTURED - ctx->body.size());
		ctx->body.append(data, to_take);
	}
	// Always claim full consumption — libcurl will retry the write
	// otherwise, which we don't want.
	return total;
}

} // namespace

std::string GetCurlVersion() {
	// curl_version() returns a static string; it does NOT require
	// curl_global_init. Skipping the init keeps `miint_versions()` cheap.
	const char *v = curl_version();
	return v ? std::string(v) : std::string();
}

CurlUploadResult RunCurlUpload(const CurlUploadOptions &opts, const BodyProducer &producer) {
	EnsureCurlGlobalInit();

	CurlUploadResult out;

	if (opts.url.empty()) {
		out.error_message = "RunCurlUpload: url is empty";
		return out;
	}
	// USERPWD splits on the FIRST colon — a colon in the userid would
	// silently let a server authenticate the wrong identity (RFC 7617 §2.1
	// forbids it for HTTP Basic; FTP has the same parsing). Reject early
	// rather than corrupt the auth handshake. Mirrors the validation in
	// `BuildBasicAuthHeader` (src/include/ena_client.hpp).
	if (opts.user.find(':') != std::string::npos) {
		out.error_message = "RunCurlUpload: user must not contain ':' (RFC 7617 §2.1 / FTP USER)";
		return out;
	}

	CURL *curl = curl_easy_init();
	if (!curl) {
		throw std::runtime_error("libcurl: curl_easy_init returned nullptr");
	}

	ReadContext rctx;
	rctx.producer = &producer;
	WriteContext wctx;

	char errbuf[CURL_ERROR_SIZE];
	errbuf[0] = '\0';

	curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, errbuf);
	curl_easy_setopt(curl, CURLOPT_URL, opts.url.c_str());
	curl_easy_setopt(curl, CURLOPT_UPLOAD, 1L);
	curl_easy_setopt(curl, CURLOPT_READFUNCTION, &ReadCallback);
	curl_easy_setopt(curl, CURLOPT_READDATA, &rctx);
	// Capture the server's response body. Without this, libcurl writes it
	// to stdout. Combined with checking http_code post-perform, this gives
	// us ENA's actual XML/JSON error diagnostic (e.g. "Sample alias already
	// exists") instead of a generic "HTTP 4xx returned" surface.
	curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, &WriteCallback);
	curl_easy_setopt(curl, CURLOPT_WRITEDATA, &wctx);
	curl_easy_setopt(curl, CURLOPT_NOSIGNAL, 1L); // don't install SIGPIPE/SIGALRM handlers
	if (opts.timeout_seconds > 0) {
		curl_easy_setopt(curl, CURLOPT_TIMEOUT, opts.timeout_seconds);
	}
	if (opts.create_dirs) {
		// FTP-only option; libcurl ignores it for HTTP. Tells the server-
		// side handling to MKD intermediate directories.
		curl_easy_setopt(curl, CURLOPT_FTP_CREATE_MISSING_DIRS, 1L);
	}
	if (opts.expected_size > 0) {
		curl_easy_setopt(curl, CURLOPT_INFILESIZE_LARGE, static_cast<curl_off_t>(opts.expected_size));
	}
	if (!opts.user.empty() || !opts.password.empty()) {
		// USERPWD format: "user:password". libcurl copies the string into
		// its internal options buffer (since 7.17), so the temporary
		// lifetime is fine.
		const std::string userpwd = opts.user + ":" + opts.password;
		curl_easy_setopt(curl, CURLOPT_USERPWD, userpwd.c_str());
	}

	const CURLcode rc = curl_easy_perform(curl);

	long http_code = 0;
	curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &http_code);
	out.http_code = http_code;

	curl_easy_cleanup(curl);

	if (rctx.captured) {
		// Producer threw — re-raise. The libcurl error message is
		// secondary in this case (we caused the abort).
		std::rethrow_exception(rctx.captured);
	}

	if (rc != CURLE_OK) {
		// Prefer the per-handle error buffer if libcurl filled it; fall
		// back to the generic strerror.
		const std::size_t buf_len = std::strlen(errbuf);
		out.error_message = buf_len > 0 ? std::string(errbuf, buf_len) : std::string(curl_easy_strerror(rc));
		// Append the captured body if any — it's where ENA puts the
		// human-readable failure reason for HTTP errors.
		if (!wctx.body.empty()) {
			out.error_message += " | response body: " + wctx.body;
		}
		return out;
	}
	// Success at the libcurl level. For HTTP, treat 4xx/5xx as a
	// failure and surface the body. (Without CURLOPT_FAILONERROR we have
	// to do this ourselves — the upside is the body is now in our hands.)
	if (http_code >= 400) {
		out.error_message = "HTTP " + std::to_string(http_code);
		if (!wctx.body.empty()) {
			out.error_message += " response body: " + wctx.body;
		}
	}
	return out;
}

} // namespace miint
