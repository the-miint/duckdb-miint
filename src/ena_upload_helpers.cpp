// SPDX-License-Identifier: MIT
//
// Implementation of the pure-data helpers behind `ena_upload_reads`. See
// ena_upload_helpers.hpp for the contract; tests in
// test/cpp/test_ena_upload_reads.cpp pin every branch.

#include "ena_upload_helpers.hpp"

#include <cctype>
#include <stdexcept>
#include <string>
#include <vector>

namespace duckdb {

namespace {

std::string ToLowerCopy(const std::string &s) {
	std::string out;
	out.reserve(s.size());
	for (char c : s) {
		out.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
	}
	return out;
}

} // namespace

const char *FastqLayoutModeName(FastqLayoutMode mode) {
	switch (mode) {
	case FastqLayoutMode::AUTO:
		return "auto";
	case FastqLayoutMode::SINGLE:
		return "single";
	case FastqLayoutMode::PAIRED:
		return "paired";
	case FastqLayoutMode::PAIRED_INTERLEAVED:
		return "paired_interleaved";
	}
	// Unreachable when callers stay in-enum; throw rather than silently emitting
	// a default so an enum addition fails loudly.
	throw std::logic_error("FastqLayoutModeName: unhandled FastqLayoutMode");
}

FastqLayoutMode ParseFastqLayoutMode(const std::string &name) {
	const auto lower = ToLowerCopy(name);
	if (lower == "auto") {
		return FastqLayoutMode::AUTO;
	}
	if (lower == "single") {
		return FastqLayoutMode::SINGLE;
	}
	if (lower == "paired") {
		return FastqLayoutMode::PAIRED;
	}
	if (lower == "paired_interleaved") {
		return FastqLayoutMode::PAIRED_INTERLEAVED;
	}
	throw std::runtime_error("Unknown layout '" + name +
	                         "' — expected one of: auto, single, paired, paired_interleaved");
}

FastqLayoutMode ResolveLayoutFromCounts(const std::string &sample_ref, FastqLayoutMode requested, bool all_paired,
                                        bool any_paired) {
	// A sample reaching here always has >= 1 row (it came from a GROUP BY), so
	// the empty-set degenerate case (all_paired=true, any_paired=false) cannot
	// occur. The decision is driven by (bool_and, bool_or) of "sequence2 IS NOT
	// NULL" — so the errors name the sample but cannot point at an offending row.
	const bool all_single = !any_paired;

	switch (requested) {
	case FastqLayoutMode::AUTO:
		if (all_single) {
			return FastqLayoutMode::SINGLE;
		}
		if (all_paired) {
			return FastqLayoutMode::PAIRED;
		}
		throw std::runtime_error("Sample '" + sample_ref +
		                         "' mixes single-end and paired-end rows. Specify layout explicitly or "
		                         "split the input into per-sample groups.");

	case FastqLayoutMode::SINGLE:
		if (all_single) {
			return FastqLayoutMode::SINGLE;
		}
		throw std::runtime_error("Sample '" + sample_ref + "' has rows with non-null R2 but layout=single requested");

	case FastqLayoutMode::PAIRED:
		if (all_paired) {
			return FastqLayoutMode::PAIRED;
		}
		throw std::runtime_error("Sample '" + sample_ref + "' has rows missing R2 but layout=paired requested");

	case FastqLayoutMode::PAIRED_INTERLEAVED:
		if (all_paired) {
			return FastqLayoutMode::PAIRED_INTERLEAVED;
		}
		throw std::runtime_error("Sample '" + sample_ref +
		                         "' has rows missing R2 but layout=paired_interleaved requested");
	}
	throw std::logic_error("ResolveLayoutFromCounts: unhandled FastqLayoutMode");
}

std::vector<std::string> OutputFilenames(const std::string &sample_ref, FastqLayoutMode layout) {
	switch (layout) {
	case FastqLayoutMode::SINGLE:
	case FastqLayoutMode::PAIRED_INTERLEAVED:
		return {sample_ref + ".fastq.gz"};
	case FastqLayoutMode::PAIRED:
		return {sample_ref + "_1.fastq.gz", sample_ref + "_2.fastq.gz"};
	case FastqLayoutMode::AUTO:
		throw std::runtime_error("OutputFilenames: layout must be resolved before requesting filenames");
	}
	throw std::logic_error("OutputFilenames: unhandled FastqLayoutMode");
}

namespace {

// Parse a `scheme://host[/path]` URL into its host and (trailing-slash-
// normalised) path components. Throws when the host is empty.
void ParseHostAndPath(const std::string &url, const std::string &scheme, const std::string &remainder,
                      std::string &host_out, std::string &path_out) {
	const auto slash = remainder.find('/');
	std::string host;
	std::string path;
	if (slash == std::string::npos) {
		host = remainder;
		path = "/";
	} else {
		host = remainder.substr(0, slash);
		path = remainder.substr(slash);
	}
	if (host.empty()) {
		throw std::runtime_error("Upload target URL '" + url + "' is missing " + scheme + " host");
	}
	if (path.empty() || path.front() != '/') {
		path = "/" + path;
	}
	if (path.back() != '/') {
		path.push_back('/');
	}
	host_out = std::move(host);
	path_out = std::move(path);
}

bool IsCurlScheme(const std::string &scheme) {
	return scheme == "http" || scheme == "https" || scheme == "ftp" || scheme == "ftps";
}

} // namespace

UploadTargetURL ParseUploadTargetURL(const std::string &url) {
	if (url.empty()) {
		throw std::runtime_error("Upload target URL is empty");
	}
	const auto sep = url.find("://");
	if (sep == std::string::npos) {
		throw std::runtime_error("Upload target URL '" + url + "' is missing scheme separator '://'");
	}
	const std::string scheme = ToLowerCopy(url.substr(0, sep));
	const std::string remainder = url.substr(sep + 3);

	UploadTargetURL out;
	out.scheme = scheme;

	if (scheme == "aspera") {
		ParseHostAndPath(url, scheme, remainder, out.host, out.remote_dir);
		out.transport = UploadTransport::ASPERA;
		return out;
	}
	if (scheme == "file") {
		// file:/// — accept the strict three-slash absolute form
		// (`file:///tmp/x/`), and also a forgiving relative form
		// (`file://relative/dir/`) so test fixtures using DuckDB's
		// per-test relative `__TEST_DIR__` work without contortions.
		// The remainder after `file://` is taken verbatim as the path.
		std::string path = remainder;
		if (path.empty()) {
			throw std::runtime_error("Upload target URL '" + url + "' has no path component");
		}
		if (path.back() != '/') {
			path.push_back('/');
		}
		out.transport = UploadTransport::LOCAL_FILE;
		out.host.clear();
		out.remote_dir = path;
		return out;
	}
	if (IsCurlScheme(scheme)) {
		ParseHostAndPath(url, scheme, remainder, out.host, out.remote_dir);
		out.transport = UploadTransport::CURL;
		// libcurl wants the full URL; reconstruct it after normalisation
		// so the trailing slash is consistent with our other transports.
		out.url_for_curl = scheme + "://" + out.host + out.remote_dir;
		return out;
	}
	throw std::runtime_error(
	    "Unsupported URL scheme '" + scheme +
	    "' for upload target — expected aspera://, file://, ftp://, ftps://, http://, or https://");
}

std::vector<std::string> BuildAscpSendArgv(const AsperaSendOptions &opts) {
	if (opts.ascp_path.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: ascp_path is empty");
	}
	if (opts.key_path.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: key_path is empty");
	}
	if (opts.user.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: user is empty");
	}
	if (opts.host.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: host is empty");
	}
	if (opts.local_path.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: local_path is empty");
	}
	if (opts.remote_dir.empty()) {
		throw std::runtime_error("BuildAscpSendArgv: remote_dir is empty");
	}

	std::vector<std::string> argv;
	argv.reserve(16);
	argv.push_back(opts.ascp_path);
	argv.push_back("--mode=send");
	argv.push_back("--user=" + opts.user);
	argv.push_back("--host=" + opts.host);
	argv.push_back("-P");
	argv.push_back(std::to_string(opts.port));
	argv.push_back("-i");
	argv.push_back(opts.key_path);
	argv.push_back("-Q"); // adaptive rate
	argv.push_back("-d"); // create destination directory if missing
	if (!opts.max_rate.empty()) {
		argv.push_back("-l");
		argv.push_back(opts.max_rate);
	}
	argv.push_back(opts.local_path);
	argv.push_back(opts.remote_dir);
	return argv;
}

} // namespace duckdb
