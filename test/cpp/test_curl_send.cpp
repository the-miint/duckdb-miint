// SPDX-License-Identifier: MIT
//
// Phase 6.5 unit tests for the libcurl streaming-upload helper. These tests
// run offline — no network, no test server required.
//
// Built only when `MIINT_ENABLE_CURL=ON` (default everywhere except macOS,
// where vsearch's vendored MD5/SHA1 collide with libcurl→libcrypto and
// Apple's ld lacks `--allow-multiple-definition`). On builds without
// libcurl, this file is excluded by CMakeLists.txt; the guard below is a
// belt-and-braces backup if a future refactor unconditionally adds it.
//   - GetCurlVersion returns a non-empty string starting with "libcurl"
//   - RunCurlUpload to a non-existent local URL returns an error_message
//     (proves the wrapper actually invokes libcurl + surfaces failures)
//   - The producer is consulted before libcurl reports the connection error
//     when the body is needed, but with a non-resolving DNS we may not get
//     that far — the important invariant is that producer-thrown exceptions
//     propagate cleanly, not that libcurl always asks for body data.
//   - Producer-thrown exceptions propagate out of RunCurlUpload (caller can
//     `try/catch` around it without leaking).
//   - An empty URL is rejected with a clear error_message (does NOT throw).
//
// Live HTTPS-PUT / FTP integration is exercised via SQL when a real
// `target_url` of those schemes is wired into ena_upload_reads in a
// follow-up commit.

#ifdef MIINT_HAS_CURL

#include "curl_send.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <string>

using miint::BodyProducer;
using miint::CurlUploadOptions;
using miint::CurlUploadResult;
using miint::GetCurlVersion;
using miint::RunCurlUpload;

TEST_CASE("GetCurlVersion returns a non-empty libcurl identifier", "[curl_send]") {
	const auto v = GetCurlVersion();
	REQUIRE_FALSE(v.empty());
	REQUIRE_THAT(v, Catch::Matchers::ContainsSubstring("libcurl"));
}

TEST_CASE("RunCurlUpload: empty URL is rejected without invoking libcurl", "[curl_send]") {
	CurlUploadOptions opts;
	bool producer_called = false;
	auto producer = [&](char *, std::size_t) -> std::size_t {
		producer_called = true;
		return 0;
	};
	auto result = RunCurlUpload(opts, producer);
	REQUIRE_FALSE(result.error_message.empty());
	REQUIRE_THAT(result.error_message, Catch::Matchers::ContainsSubstring("url is empty"));
	REQUIRE_FALSE(producer_called);
}

TEST_CASE("RunCurlUpload: rejects a user containing ':' (auth ambiguity)", "[curl_send]") {
	// USERPWD parses on the FIRST colon — a colon in the userid would
	// silently authenticate as the wrong identity (RFC 7617 §2.1 and FTP
	// USER both forbid). RunCurlUpload validates this before invoking
	// libcurl rather than relying on the upstream secret type.
	CurlUploadOptions opts;
	opts.url = "http://example.org/";
	opts.user = "Webin-12345:bad";
	opts.password = "secret";
	bool producer_called = false;
	auto producer = [&](char *, std::size_t) -> std::size_t {
		producer_called = true;
		return 0;
	};
	auto result = RunCurlUpload(opts, producer);
	REQUIRE_THAT(result.error_message, Catch::Matchers::ContainsSubstring("RFC 7617"));
	REQUIRE_FALSE(producer_called);
}

TEST_CASE("RunCurlUpload: bogus host surfaces libcurl error", "[curl_send]") {
	// `127.0.0.1:1` should refuse all incoming connections (port 1 is
	// privileged + nothing listens). Some CI sandboxes block outbound
	// network entirely; either way libcurl returns CURLE_COULDNT_CONNECT
	// or similar, populating error_message.
	CurlUploadOptions opts;
	opts.url = "http://127.0.0.1:1/";
	opts.timeout_seconds = 3;
	auto producer = [](char *buf, std::size_t max_bytes) -> std::size_t {
		(void)buf;
		(void)max_bytes;
		return 0; // EOF; should not be reached if connection is refused first
	};
	auto result = RunCurlUpload(opts, producer);
	REQUIRE_FALSE(result.error_message.empty());
}

TEST_CASE("RunCurlUpload: producer-thrown exception propagates to caller", "[curl_send]") {
	CurlUploadOptions opts;
	// Use a URL libcurl will *try* to upload to; an immediate connection
	// failure could end the request before the producer is consulted.
	// `file://` upload is supported by libcurl and exercises the read
	// callback synchronously, so a deliberately-thrown exception in the
	// producer must propagate out.
	opts.url = "file:///tmp/__miint_curl_send_test_should_not_be_created";
	bool producer_called = false;
	auto producer = [&](char *, std::size_t) -> std::size_t {
		producer_called = true;
		throw std::runtime_error("ena_upload_reads sentinel: producer failure");
	};
	REQUIRE_THROWS_WITH(RunCurlUpload(opts, producer), Catch::Matchers::ContainsSubstring("ena_upload_reads sentinel"));
	REQUIRE(producer_called);
}

TEST_CASE("RunCurlUpload: file:// transport streams producer bytes through libcurl", "[curl_send]") {
	// End-to-end happy path that needs no network: upload to a file:// URL
	// in a temp location, then read the file back and assert the bytes
	// match. Proves libcurl ↔ producer integration without a server.
	const std::string tmp_path = "/tmp/__miint_curl_send_test_payload.bin";
	std::remove(tmp_path.c_str()); // clean prior runs

	CurlUploadOptions opts;
	opts.url = "file://" + tmp_path;

	const std::string payload = "miint-curl-send-streaming-test\nline2\n";
	std::size_t cursor = 0;
	auto producer = [&](char *buf, std::size_t max_bytes) -> std::size_t {
		const std::size_t remaining = payload.size() - cursor;
		if (remaining == 0) {
			return 0;
		}
		const std::size_t to_copy = std::min(remaining, max_bytes);
		std::memcpy(buf, payload.data() + cursor, to_copy);
		cursor += to_copy;
		return to_copy;
	};

	auto result = RunCurlUpload(opts, producer);
	REQUIRE(result.error_message.empty());

	// Read back and verify.
	FILE *f = std::fopen(tmp_path.c_str(), "rb");
	REQUIRE(f != nullptr);
	std::string read_back;
	char rbuf[4096];
	while (true) {
		auto n = std::fread(rbuf, 1, sizeof(rbuf), f);
		if (n == 0) {
			break;
		}
		read_back.append(rbuf, n);
	}
	std::fclose(f);
	std::remove(tmp_path.c_str());

	REQUIRE(read_back == payload);
}

#endif // MIINT_HAS_CURL
