#pragma once

#include "ena_parser.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/main/database.hpp"
#include <chrono>
#include <mutex>
#include <stdexcept>
#include <string>
#include <vector>

namespace miint {

// Standard ("plus") Base64 encoder. Self-contained so the unit-test target
// can call it without linking against duckdb. Output is canonical (uses '+' /
// '/', pads to a multiple of 4 with '='). Input is treated as raw bytes;
// no UTF-8 / signed-char surprise.
inline std::string Base64Encode(const std::string &input) {
	static const char ALPHABET[] = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
	const auto n = input.size();
	std::string out;
	out.reserve(((n + 2) / 3) * 4);

	std::size_t i = 0;
	while (i + 3 <= n) {
		const auto a = static_cast<unsigned char>(input[i]);
		const auto b = static_cast<unsigned char>(input[i + 1]);
		const auto c = static_cast<unsigned char>(input[i + 2]);
		out.push_back(ALPHABET[(a >> 2) & 0x3F]);
		out.push_back(ALPHABET[((a << 4) | (b >> 4)) & 0x3F]);
		out.push_back(ALPHABET[((b << 2) | (c >> 6)) & 0x3F]);
		out.push_back(ALPHABET[c & 0x3F]);
		i += 3;
	}
	const auto remainder = n - i;
	if (remainder == 1) {
		const auto a = static_cast<unsigned char>(input[i]);
		out.push_back(ALPHABET[(a >> 2) & 0x3F]);
		out.push_back(ALPHABET[(a << 4) & 0x3F]);
		out.push_back('=');
		out.push_back('=');
	} else if (remainder == 2) {
		const auto a = static_cast<unsigned char>(input[i]);
		const auto b = static_cast<unsigned char>(input[i + 1]);
		out.push_back(ALPHABET[(a >> 2) & 0x3F]);
		out.push_back(ALPHABET[((a << 4) | (b >> 4)) & 0x3F]);
		out.push_back(ALPHABET[(b << 2) & 0x3F]);
		out.push_back('=');
	}
	return out;
}

// Build an RFC 7617 Basic auth header value: "Basic " + Base64(user:password).
// Throws std::runtime_error if either user or password is empty, or if user
// contains a colon (RFC 7617 §2.1 forbids that — the server would split the
// userid on the first colon and silently authenticate as someone else).
//
// Note: credentials live in heap-allocated std::string copies in the caller,
// here, and inside Base64Encode. miint runs as a single-tenant CLI/embedded
// extension so we don't go to the trouble of wiping memory; if you are using
// this in a multitenant context, scrub explicitly.
inline std::string BuildBasicAuthHeader(const std::string &user, const std::string &password) {
	if (user.empty()) {
		throw std::runtime_error("Basic auth: user must be non-empty");
	}
	if (password.empty()) {
		throw std::runtime_error("Basic auth: password must be non-empty");
	}
	if (user.find(':') != std::string::npos) {
		throw std::runtime_error("Basic auth: user must not contain ':' (RFC 7617 §2.1)");
	}
	std::string credentials = user + ":" + password;
	return "Basic " + Base64Encode(credentials);
}

// HTTP client for EBI/ENA APIs.
// For parsing utilities, use ENAParser directly (ena_parser.hpp).
//
// Rate Limiting: ~3 requests/second (no API key required)
// Retry: HTTP 429, 500, 502, 503 with exponential backoff
class ENAClient {
public:
	static constexpr double RATE_LIMIT = 3.0;
	static constexpr int MAX_RETRIES = 3;
	static constexpr int INITIAL_RETRY_DELAY_MS = 1000;

	// Inline so unit tests can call it without linking ena_client.cpp
	// (which depends on duckdb symbols not present in the test target).
	// 504 (Gateway Timeout) is included alongside 502/503 because ENA's
	// frontend occasionally returns it for slow validation passes.
	static inline bool IsRetryableStatus(int status) {
		return status == 429 || status == 500 || status == 502 || status == 503 || status == 504;
	}

	explicit ENAClient(duckdb::DatabaseInstance &db);

	// HTTP methods — build URL via ENAParser, then fetch
	std::string Search(const std::string &accession, const std::string &result_type, const std::string &fields);
	std::string FetchXML(const std::vector<std::string> &accessions);

	// Fetch an arbitrary URL with the same rate-limit + retry semantics as Search/FetchXML.
	// Intended for callers that build URLs externally (e.g., ResolveRunsBatch with compound
	// "accession IN (...)" queries).
	std::string FetchURL(const std::string &url);

	// POST a JSON or XML body with HTTP Basic authentication. Returns the
	// response body on success. Retries 429/5xx with exponential backoff
	// matching the GET path. Throws on auth failure (401/403), other 4xx,
	// or transport error. The Webin V2 endpoints all use Basic auth (see
	// localdocs/ena-research-webin-v2-deep.md §3).
	//
	// `accept_type` is the response format we ask for. Webin V2 honours the
	// Accept header to select XML vs JSON response shape regardless of the
	// request body format, so callers can mix (e.g. POST JSON body, get XML
	// receipt). PostJSON/PostXML set Accept = Content-Type;
	// PostJSONReceiveXML is the V2-Webin pairing (JSON envelope, XML receipt).
	std::string PostJSON(const std::string &url, const std::string &body, const std::string &user,
	                     const std::string &password);

	std::string PostXML(const std::string &url, const std::string &body, const std::string &user,
	                    const std::string &password);

	// Webin V2 envelopes are JSON but the receipt parser only consumes the
	// canonical XSD-governed XML form (JSON receipts intentionally not
	// supported — see Phase 3 design). Use this for the INSERT operators.
	std::string PostJSONReceiveXML(const std::string &url, const std::string &body, const std::string &user,
	                               const std::string &password);

	// Authenticated GET. Used by the pre-INSERT alias collision check to
	// search the user's submission account on the ENA portal API
	// (anonymous portal search would not see HOLD/private records owned by
	// the authenticated submission account). Same retry/backoff semantics
	// as the unauthenticated GET path.
	std::string AuthenticatedGet(const std::string &url, const std::string &user, const std::string &password);

private:
	duckdb::DatabaseInstance &db;
	std::chrono::steady_clock::time_point last_request_time;
	mutable std::mutex rate_limit_mutex;

	std::string MakeRequest(const std::string &url);
	// Shared GET retry loop. `headers` is built by the caller (empty for
	// anonymous, Authorization-set for authenticated) so we share the
	// rate-limit + retry + backoff semantics across both paths.
	std::string ExecuteGet(const std::string &url, duckdb::HTTPHeaders &headers, const char *failure_label);
	// content_type and accept_type are independent (mismatching is allowed,
	// e.g. POST XML, ask for JSON receipt — Webin V2 supports this).
	std::string PostBody(const std::string &url, const std::string &body, const std::string &content_type,
	                     const std::string &accept_type, const std::string &user, const std::string &password);
	void RespectRateLimit();
};

} // namespace miint
