// SPDX-License-Identifier: MIT
//
// Webin Reports API client — alias → primary-accession lookup.
// Pure-data implementation: URL build + HTTP GET (via fetcher) + JSON parse.
// No DuckDB linkage; tests substitute a mock fetcher.

#include "ena_reports_client.hpp"

#include "ena_parser.hpp" // PercentEncodeQueryValue

#include <cstdlib>
#include <stdexcept>
#include <string>

namespace miint {

namespace {

// Strip a single trailing slash from a base URL so callers passing
// "https://.../report" and "https://.../report/" produce identical URLs.
// Mirrors `TrimTrailingSlash` in ena_alias_check.cpp; not extracted because
// the two TUs don't currently share a header for misc helpers.
std::string TrimTrailingSlash(std::string s) {
	if (!s.empty() && s.back() == '/') {
		s.pop_back();
	}
	return s;
}

// Path segment per kind. STUDY → "projects" (returns id=PRJEB, the primary
// accession lifecycle ops target). NOT "studies" (which returns id=ERP, the
// secondary). See header for the full mapping rationale.
const char *ReportsPathForKind(AliasObjectKind kind) {
	switch (kind) {
	case AliasObjectKind::STUDY:
		return "projects";
	case AliasObjectKind::SAMPLE:
		return "samples";
	case AliasObjectKind::EXPERIMENT:
		return "experiments";
	case AliasObjectKind::RUN:
		return "runs";
	}
	throw std::logic_error("unhandled AliasObjectKind in ReportsPathForKind");
}

void ValidateId(const std::string &id) {
	if (id.empty()) {
		throw std::invalid_argument("Reports lookup id must be non-empty");
	}
	for (unsigned char c : id) {
		// Header-injection guard: HTTP requests would split on a newline.
		// PercentEncodeQueryValue would tunnel these as %0A / %0D safely,
		// but a defensive caller-side reject is clearer at debugging time.
		if (c == '\n' || c == '\r') {
			throw std::invalid_argument("Reports lookup id must not contain newline characters");
		}
	}
}

// Skip ASCII whitespace (space, \t, \n, \r) starting at `pos`. Used by the
// hand-rolled JSON tolerator below.
void SkipWs(const std::string &s, size_t &pos) {
	while (pos < s.size()) {
		const auto c = static_cast<unsigned char>(s[pos]);
		if (c != ' ' && c != '\t' && c != '\n' && c != '\r') {
			break;
		}
		pos++;
	}
}

// Read a JSON string starting at `pos` (which must point at the opening `"`).
// Advances pos past the closing quote. Returns the unquoted contents.
// Doesn't decode escapes (\\, \", \n, \uXXXX) — adequate for the Reports
// payload's `id` field which is `[A-Z0-9]+` (server-assigned accession),
// and acceptable for the `key` parsing which is also escape-free in the
// pinned wire shape. Backslash-escaped quotes inside strings are honoured
// while scanning past values (so we don't mis-terminate on a JSON-escaped
// quote in some other field).
std::string ReadJsonString(const std::string &body, size_t &pos) {
	if (pos >= body.size() || body[pos] != '"') {
		throw std::runtime_error("Reports JSON: expected '\"'");
	}
	pos++; // past opening '"'
	const auto start = pos;
	while (pos < body.size() && body[pos] != '"') {
		if (body[pos] == '\\' && pos + 1 < body.size()) {
			pos += 2;
			continue;
		}
		pos++;
	}
	if (pos >= body.size()) {
		throw std::runtime_error("Reports JSON: unterminated string");
	}
	std::string out = body.substr(start, pos - start);
	pos++; // past closing '"'
	return out;
}

// Skip exactly one complete JSON value at `pos`. Handles all four value
// classes — string, object, array, bare scalar (number / true / false /
// null). Pos advances past the end of the value (NOT past any trailing
// separator). Used by ExtractFirstReportId to skip values whose key isn't
// `report` / `id`.
//
// Earlier impl was depth-tracking with a "consumed_first_token" gate that
// silently mis-skipped bare scalars at depth 0; rewrite is correct by
// construction (each branch consumes exactly one value).
void SkipJsonValue(const std::string &body, size_t &pos) {
	SkipWs(body, pos);
	if (pos >= body.size()) {
		throw std::runtime_error("Reports JSON: EOF where value expected");
	}
	const auto first = body[pos];
	if (first == '"') {
		(void)ReadJsonString(body, pos);
		return;
	}
	if (first == '{' || first == '[') {
		const auto open = first;
		const auto close = (first == '{') ? '}' : ']';
		pos++; // past opening
		int depth = 1;
		while (pos < body.size() && depth > 0) {
			const auto c = body[pos];
			if (c == '"') {
				(void)ReadJsonString(body, pos);
				continue;
			}
			if (c == open) {
				depth++;
			} else if (c == close) {
				depth--;
			}
			pos++;
		}
		if (depth != 0) {
			throw std::runtime_error("Reports JSON: unterminated object/array");
		}
		return;
	}
	// Bare scalar: number / true / false / null. Consume until we hit a
	// JSON structural character that ends the value.
	while (pos < body.size()) {
		const auto c = body[pos];
		if (c == ',' || c == '}' || c == ']' || c == ' ' || c == '\t' || c == '\n' || c == '\r') {
			break;
		}
		pos++;
	}
}

// Find a string-typed field by name inside the open object at the current
// cursor position. Returns the value when found; throws if the key isn't
// present (the caller — ExtractFirstReportId — surfaces that as a clean
// "object missing key 'id'" diagnostic). Skips non-string values via
// SkipJsonValue, so other-typed siblings don't trip parsing.
std::string FindStringValueByKey(const std::string &body, size_t &pos, const std::string &target_key) {
	if (pos >= body.size() || body[pos] != '{') {
		throw std::runtime_error("Reports JSON: expected '{'");
	}
	pos++; // past '{'
	while (true) {
		SkipWs(body, pos);
		if (pos >= body.size()) {
			throw std::runtime_error("Reports JSON: EOF while scanning object");
		}
		if (body[pos] == '}') {
			pos++; // past '}'
			throw std::runtime_error("Reports JSON: object missing key '" + target_key + "'");
		}
		const auto key = ReadJsonString(body, pos);
		SkipWs(body, pos);
		if (pos >= body.size() || body[pos] != ':') {
			throw std::runtime_error("Reports JSON: expected ':' after key");
		}
		pos++; // past ':'
		SkipWs(body, pos);
		if (key == target_key) {
			return ReadJsonString(body, pos);
		}
		SkipJsonValue(body, pos);
		SkipWs(body, pos);
		if (pos < body.size() && body[pos] == ',') {
			pos++;
		}
	}
}

// Position pos at the start of the value for `target_key` inside the open
// object at the current cursor. Used when the value is itself an object
// (e.g. `report` → `{...}`) and the caller wants to descend into it
// without copying the value first.
void DescendToObjectKey(const std::string &body, size_t &pos, const std::string &target_key) {
	if (pos >= body.size() || body[pos] != '{') {
		throw std::runtime_error("Reports JSON: expected '{'");
	}
	pos++; // past '{'
	while (true) {
		SkipWs(body, pos);
		if (pos >= body.size()) {
			throw std::runtime_error("Reports JSON: EOF while scanning object");
		}
		if (body[pos] == '}') {
			throw std::runtime_error("Reports JSON: object missing key '" + target_key + "'");
		}
		const auto key = ReadJsonString(body, pos);
		SkipWs(body, pos);
		if (pos >= body.size() || body[pos] != ':') {
			throw std::runtime_error("Reports JSON: expected ':' after key");
		}
		pos++; // past ':'
		SkipWs(body, pos);
		if (key == target_key) {
			return; // pos points at start of value
		}
		SkipJsonValue(body, pos);
		SkipWs(body, pos);
		if (pos < body.size() && body[pos] == ',') {
			pos++;
		}
	}
}

// Extract the first `id` field inside the first `report` object of the JSON
// array. The full Reports payload contains several string-typed fields (`id`,
// `alias`, `firstCreated`, `releaseStatus`, …); we only need `id`. Returns
// the id on hit, empty string on miss (response was `[]`).
//
// Pinned shapes (live-verified):
//   hit:  `[{"report":{"id":"PRJEB1","alias":"…", …}}]`
//   miss: `[]`
// Whitespace tolerated everywhere; `id` may appear at any position within
// `report` (FindStringValueByKey scans keys until it finds the target).
std::string ExtractFirstReportId(const std::string &body) {
	size_t pos = 0;
	SkipWs(body, pos);
	if (pos >= body.size() || body[pos] != '[') {
		throw std::runtime_error("Reports JSON: expected array");
	}
	pos++; // past '['
	SkipWs(body, pos);
	if (pos < body.size() && body[pos] == ']') {
		return ""; // empty array — miss
	}
	// Descend into the first array element's `report` object, then read its
	// `id` field. Both helpers throw on missing keys, which matches the
	// "hit array element missing report.id is a parse error" test expectation.
	DescendToObjectKey(body, pos, "report");
	return FindStringValueByKey(body, pos, "id");
}

} // namespace

std::string BuildReportsLookupURL(const std::string &reports_base, AliasObjectKind kind, const std::string &id) {
	ValidateId(id);
	std::string url;
	url.reserve(reports_base.size() + 64 + id.size());
	url += TrimTrailingSlash(reports_base);
	url += '/';
	url += ReportsPathForKind(kind);
	url += '/';
	// PercentEncodeQueryValue encodes everything outside the RFC 3986
	// unreserved set (alnum + `-._~`). Path segments allow a slightly larger
	// set (sub-delims + `:@`), but the unreserved-only encoding is strictly
	// safer and produces a URL the server still parses correctly — verified
	// in the test for spaces / ampersands / equals.
	url += PercentEncodeQueryValue(id);
	url += "?format=json";
	return url;
}

std::string LookupAccessionByAlias(const std::string &reports_base, AliasObjectKind kind, const std::string &alias,
                                   const URLFetcher &fetcher) {
	if (alias.empty()) {
		throw std::invalid_argument("Reports lookup alias must be non-empty");
	}
	if (!fetcher) {
		throw std::invalid_argument("Reports lookup fetcher must be set");
	}
	const auto url = BuildReportsLookupURL(reports_base, kind, alias);
	const auto body = fetcher(url);
	return ExtractFirstReportId(body);
}

std::string ResolveReportsBaseForEndpoint(const std::string &endpoint) {
	const char *override_base = std::getenv("MIINT_ENA_REPORTS_URL_BASE");
	if (override_base != nullptr && *override_base != '\0') {
		return TrimTrailingSlash(override_base);
	}
	// Mirror ResolveDefaultENAEndpointURL's mapping: only "production"
	// routes to the prod base; everything else (including "test", empty,
	// or any other label) defaults to wwwdev. The lifecycle POST goes to
	// the same endpoint by the same rule, so the Reports lookup and the
	// lifecycle action stay aligned.
	const auto base = (endpoint == "production") ? DEFAULT_REPORTS_BASE_PROD : DEFAULT_REPORTS_BASE_TEST;
	return TrimTrailingSlash(base);
}

} // namespace miint
