// SPDX-License-Identifier: MIT
//
// Pre-INSERT alias collision check (Phase 8 Step 8a).
// Pure-data implementation: URL build + chunked fetch + TSV parse. No DuckDB
// linkage, so test/cpp/test_ena_alias_check.cpp can substitute a mock fetcher.

#include "ena_alias_check.hpp"

#include <cstdlib>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace miint {

namespace {

// One portal-API query batches up to this many aliases per URL. With each
// alias <80 chars typical, the resulting URL stays well under any common
// 4KB / 8KB limit. Keep modest for clearer test failures and so a single
// stale alias doesn't fan out into a giant chunk.
constexpr size_t ALIAS_QUERY_CHUNK_SIZE = 50;

void ValidateAlias(const std::string &alias) {
	if (alias.empty()) {
		throw std::invalid_argument("alias must be non-empty");
	}
	for (unsigned char c : alias) {
		// Newline / quote would let a crafted alias break out of the IN-list
		// quoting structure even after percent encoding (the portal API
		// URL-decodes the query parameter before parsing). Reject up front
		// rather than relying on the encoder to tunnel them safely.
		if (c == '"' || c == '\n' || c == '\r') {
			throw std::invalid_argument("alias must not contain '\"' or newline characters");
		}
	}
}

// Percent-encode a value to be safe inside a URL query parameter. Encodes
// every byte except the RFC 3986 unreserved set (alnum + `-._~`). Aliases
// commonly contain `_`, `-`, `.` and digits, so the typical alias passes
// through verbatim. Anything else, including `%`, `&`, `=`, `+`, `?`, `#`,
// `/`, and any non-ASCII byte, is %HH-encoded.
std::string PercentEncodeQueryValue(const std::string &value) {
	static const char HEX[] = "0123456789ABCDEF";
	std::string out;
	out.reserve(value.size() + 8);
	for (unsigned char c : value) {
		const bool unreserved = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9') ||
		                        c == '-' || c == '_' || c == '.' || c == '~';
		if (unreserved) {
			out.push_back(static_cast<char>(c));
		} else {
			out.push_back('%');
			out.push_back(HEX[c >> 4]);
			out.push_back(HEX[c & 0x0F]);
		}
	}
	return out;
}

// URL-encoded literals used by the IN-list query, mirroring the
// `<field> IN (<csv>)` shape that ENAParser::BuildSearchURLBatch produces:
//   space = %20   ( = %28   ) = %29   , = %2C   " = %22
// Each alias value is percent-encoded so `%`, `&`, `=`, `+` and other
// reserved characters don't corrupt the query when the portal API decodes
// the URL.
std::string BuildInList(const std::vector<std::string> &aliases) {
	std::ostringstream buf;
	buf << "%20IN%20%28";
	for (size_t i = 0; i < aliases.size(); i++) {
		if (i > 0) {
			buf << "%2C";
		}
		buf << "%22" << PercentEncodeQueryValue(aliases[i]) << "%22";
	}
	buf << "%29";
	return buf.str();
}

// Inline TSV parser: first row is column names; subsequent non-empty rows are
// data. Tolerant of CRLF and trailing blank lines. Alias rows have a single
// column, so we just take field[0] of each row. Same shape as
// ENAParser::ParseTSV (src/ena_parser.cpp:154) but kept local so this file
// has no extra link dependencies.
std::vector<std::string> ParseAliasTSV(const std::string &tsv) {
	std::vector<std::string> aliases;
	if (tsv.empty()) {
		return aliases;
	}
	std::istringstream stream(tsv);
	std::string line;
	if (!std::getline(stream, line)) {
		return aliases;
	}
	// Skip the header line.
	while (std::getline(stream, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (line.empty()) {
			continue;
		}
		const auto tab = line.find('\t');
		aliases.push_back(tab == std::string::npos ? line : line.substr(0, tab));
	}
	return aliases;
}

// Strip a single trailing slash from a base URL so callers passing
// "https://.../api" and "https://.../api/" produce identical query URLs.
std::string TrimTrailingSlash(std::string s) {
	if (!s.empty() && s.back() == '/') {
		s.pop_back();
	}
	return s;
}

} // namespace

const char *AliasObjectKindName(AliasObjectKind kind) {
	switch (kind) {
	case AliasObjectKind::STUDY:
		return "study";
	case AliasObjectKind::SAMPLE:
		return "sample";
	case AliasObjectKind::EXPERIMENT:
		return "experiment";
	case AliasObjectKind::RUN:
		return "run";
	}
	throw std::logic_error("unhandled AliasObjectKind in AliasObjectKindName");
}

AliasObjectKind AliasObjectKindFromTableName(const std::string &table_name) {
	if (table_name == "projects") {
		return AliasObjectKind::STUDY;
	}
	if (table_name == "samples") {
		return AliasObjectKind::SAMPLE;
	}
	if (table_name == "experiments") {
		return AliasObjectKind::EXPERIMENT;
	}
	if (table_name == "runs") {
		return AliasObjectKind::RUN;
	}
	throw std::invalid_argument("AliasObjectKindFromTableName: unknown table '" + table_name + "'");
}

std::string BuildAliasCollisionURL(const std::string &portal_base, AliasObjectKind kind,
                                   const std::vector<std::string> &aliases) {
	if (aliases.empty()) {
		throw std::invalid_argument("BuildAliasCollisionURL: aliases vector must not be empty");
	}
	for (const auto &a : aliases) {
		ValidateAlias(a);
	}
	const std::string kind_name = AliasObjectKindName(kind);
	std::ostringstream url;
	url << TrimTrailingSlash(portal_base) << "/search?"
	    << "result=read_" << kind_name << "&query=" << kind_name << "_alias" << BuildInList(aliases)
	    << "&fields=" << kind_name << "_alias&limit=0&format=tsv";
	return url.str();
}

std::vector<std::string> CheckAliasCollisions(const std::string &portal_base, AliasObjectKind kind,
                                              const std::vector<std::string> &aliases, const URLFetcher &fetcher) {
	if (aliases.empty()) {
		return {};
	}
	if (!fetcher) {
		throw std::invalid_argument("CheckAliasCollisions: fetcher must be set");
	}

	// Stable input order, dedup, validation. We build a canonical position map
	// keyed by alias so the final hits vector preserves the position of the
	// first occurrence.
	std::vector<std::string> queue;
	std::unordered_map<std::string, size_t> first_pos;
	for (size_t i = 0; i < aliases.size(); i++) {
		const auto &a = aliases[i];
		ValidateAlias(a);
		if (first_pos.emplace(a, i).second) {
			queue.push_back(a);
		}
	}

	// Aggregate hits across chunks. Use a set for O(1) membership tests against
	// the user's input list; reject any alias the server names that we did not
	// query for.
	std::unordered_set<std::string> hit_set;
	for (size_t off = 0; off < queue.size(); off += ALIAS_QUERY_CHUNK_SIZE) {
		const size_t end = std::min(off + ALIAS_QUERY_CHUNK_SIZE, queue.size());
		std::vector<std::string> chunk(queue.begin() + static_cast<std::ptrdiff_t>(off),
		                               queue.begin() + static_cast<std::ptrdiff_t>(end));

		const auto url = BuildAliasCollisionURL(portal_base, kind, chunk);
		const auto body = fetcher(url);
		for (auto &reported : ParseAliasTSV(body)) {
			if (first_pos.find(reported) != first_pos.end()) {
				hit_set.insert(std::move(reported));
			}
			// Aliases the server names that we didn't query for are silently
			// ignored — the caller cannot act on them and a noisy warning
			// would be misleading (the user's request did not mention them).
		}
	}

	// Emit hits in caller's input order, deduped.
	std::vector<std::string> hits;
	hits.reserve(hit_set.size());
	std::unordered_set<std::string> emitted;
	for (const auto &a : aliases) {
		if (emitted.count(a)) {
			continue;
		}
		if (hit_set.count(a)) {
			hits.push_back(a);
			emitted.insert(a);
		}
	}
	return hits;
}

std::string ResolvePortalBaseFromEnv() {
	const char *override_base = std::getenv("MIINT_ENA_PORTAL_URL_BASE");
	const std::string base = (override_base != nullptr && *override_base != '\0') ? std::string(override_base)
	                                                                              : std::string(DEFAULT_PORTAL_BASE);
	return TrimTrailingSlash(base);
}

} // namespace miint
