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

FastqLayoutMode ResolveLayout(const std::string &sample_ref, FastqLayoutMode requested,
                              const std::vector<bool> &has_r2) {
	if (has_r2.empty()) {
		throw std::runtime_error("ResolveLayout: sample '" + sample_ref + "' has zero rows");
	}

	// Compute three per-mode "first offending row" indices. AUTO wants the
	// first row inconsistent with row 0 (either flip direction). SINGLE wants
	// the first row that *has* R2. PAIRED / PAIRED_INTERLEAVED want the first
	// row that *lacks* R2. The earlier code reused one counter and reported
	// the wrong row when the input started with R2 and later went single.
	std::size_t r2_count = 0;
	const std::size_t sentinel = has_r2.size();
	std::size_t first_auto_mismatch = sentinel;
	std::size_t first_with_r2 = sentinel;
	std::size_t first_without_r2 = sentinel;
	const bool first_has = has_r2[0];
	for (std::size_t i = 0; i < has_r2.size(); i++) {
		if (has_r2[i]) {
			r2_count++;
			if (first_with_r2 == sentinel) {
				first_with_r2 = i;
			}
		} else if (first_without_r2 == sentinel) {
			first_without_r2 = i;
		}
		if (has_r2[i] != first_has && first_auto_mismatch == sentinel) {
			first_auto_mismatch = i;
		}
	}
	const bool all_paired = r2_count == has_r2.size();
	const bool all_single = r2_count == 0;

	switch (requested) {
	case FastqLayoutMode::AUTO:
		if (all_single) {
			return FastqLayoutMode::SINGLE;
		}
		if (all_paired) {
			return FastqLayoutMode::PAIRED;
		}
		throw std::runtime_error("Sample '" + sample_ref +
		                         "' mixes single-end and paired-end rows (first mismatch at row " +
		                         std::to_string(first_auto_mismatch) +
		                         "). Specify layout explicitly or split the input into per-sample groups.");

	case FastqLayoutMode::SINGLE:
		if (all_single) {
			return FastqLayoutMode::SINGLE;
		}
		throw std::runtime_error("Sample '" + sample_ref + "' has rows with non-null R2 but layout=single requested" +
		                         " (first offending row " + std::to_string(first_with_r2) + ")");

	case FastqLayoutMode::PAIRED:
		if (all_paired) {
			return FastqLayoutMode::PAIRED;
		}
		throw std::runtime_error("Sample '" + sample_ref + "' has rows missing R2 but layout=paired requested" +
		                         " (first offending row " + std::to_string(first_without_r2) + ")");

	case FastqLayoutMode::PAIRED_INTERLEAVED:
		if (all_paired) {
			return FastqLayoutMode::PAIRED_INTERLEAVED;
		}
		throw std::runtime_error("Sample '" + sample_ref +
		                         "' has rows missing R2 but layout=paired_interleaved requested" +
		                         " (first offending row " + std::to_string(first_without_r2) + ")");
	}
	throw std::logic_error("ResolveLayout: unhandled FastqLayoutMode");
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
	if (scheme == "aspera") {
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
			throw std::runtime_error("Upload target URL '" + url + "' is missing aspera host");
		}
		if (path.empty() || path.front() != '/') {
			path = "/" + path;
		}
		if (path.back() != '/') {
			path.push_back('/');
		}
		out.transport = UploadTransport::ASPERA;
		out.host = host;
		out.remote_dir = path;
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
	throw std::runtime_error("Unsupported URL scheme '" + scheme +
	                         "' for upload target — expected aspera:// or file://");
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
