// SPDX-License-Identifier: MIT
//
// Streaming upload via libcurl. Used by `ena_upload_reads` when the
// `target_url` scheme is `https://`, `http://`, `ftp://`, or `ftps://`.
// Unlike the Aspera path (which has to write a temp file because `ascp`
// doesn't read stdin reliably), libcurl pulls bytes from a caller-supplied
// `BodyProducer` callback on demand — no intermediate file ever exists.
//
// Threading model: the caller's thread blocks inside `RunCurlUpload` while
// libcurl drives the upload. libcurl repeatedly invokes the producer in
// the same thread until the producer signals EOF (returns 0). The producer
// owns whatever encode/gzip/MD5 state it needs to manufacture the next
// chunk on demand.

#pragma once

#include <cstddef>
#include <functional>
#include <string>

namespace miint {

struct CurlUploadOptions {
	// Target URL. Scheme decides the protocol: ftp://, ftps://, http://,
	// https://. For ftp/ftps the trailing path segment is the destination
	// filename; for http/https the URL is the PUT target.
	std::string url;

	// Optional. When both are non-empty, libcurl is told to use
	// USERPWD-style auth (Basic for HTTP, plain user/pass for FTP).
	// `user` must not contain ':' — RFC 7617 §2.1 forbids it for HTTP
	// Basic, and FTP USER has the same parsing problem. RunCurlUpload
	// validates this and surfaces an error_message rather than relying
	// on the upstream secret type to enforce it.
	std::string user;
	std::string password;

	// For FTP/FTPS: have libcurl issue MKD as needed to create the
	// destination directory tree. No-op for HTTP/HTTPS PUT.
	bool create_dirs = true;

	// 0 = no overall timeout. Positive value caps the request duration.
	long timeout_seconds = 0;

	// 0 = unknown size (chunked transfer for HTTP, no SIZE preflight for
	// FTP). Positive value lets libcurl set Content-Length / FTP SIZE.
	long long expected_size = 0;
};

struct CurlUploadResult {
	long http_code = 0;        // 0 for FTP success, 2xx for HTTP success
	std::string error_message; // empty on success
};

// CURLOPT_READFUNCTION-shaped producer. Called repeatedly by libcurl until
// it returns 0 (EOF). On each call, write up to `max_bytes` bytes into
// `buf` and return the number of bytes written. Returning 0 ends the
// upload. The producer may throw — exceptions propagate out of
// RunCurlUpload after libcurl has been told to abort.
using BodyProducer = std::function<std::size_t(char *buf, std::size_t max_bytes)>;

// Block until the upload completes (or fails). Returns the HTTP code (or
// 0 for FTP) on success; returns a non-empty `error_message` on libcurl-
// or server-level failure. Throws `duckdb::IOException` only for
// catastrophic libcurl-init failures (curl_easy_init returns nullptr).
//
// The caller is responsible for any retry policy. For ENA uploads the
// initial integration treats any non-success as a hard failure surfaced
// to SQL — Aspera's natural resume isn't available here.
CurlUploadResult RunCurlUpload(const CurlUploadOptions &opts, const BodyProducer &producer);

// Returns the runtime libcurl version string (curl_version()). Used by
// `miint_versions()` so users can audit which TLS / SSH backend is linked.
std::string GetCurlVersion();

} // namespace miint
