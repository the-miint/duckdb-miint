#pragma once

#include "duckdb/common/crypto/md5.hpp"
#include "duckdb/common/exception.hpp"
#include <string>

namespace miint {

// Incrementally hashes raw bytes as they stream past (see duckdb_seq_read's
// gzip/direct-read taps) and, once the underlying source has reached its true
// EOF, compares the accumulated digest against an ENA-reported fastq_md5.
//
// Mirrors the duckdb::MD5Context idiom used by GzipMd5Stream in
// ena_upload_reads.cpp (src/ena_upload_reads.cpp:364-368), but in the
// opposite direction: that struct hashes bytes on their way INTO a gzip
// encoder for upload; this one hashes bytes coming OUT of a download and
// compares against an externally supplied digest rather than just reporting
// the one it computed.
//
// An empty `expected_md5_hex` means ENA did not report an md5 for this file
// (or the caller chose not to verify it) -- both Add and VerifyOrThrow become
// no-ops so callers can construct a StreamMd5 unconditionally and let the
// empty case fall through cheaply rather than branching at every call site.
class StreamMd5 {
public:
	explicit StreamMd5(std::string expected_md5_hex) : expected_(std::move(expected_md5_hex)) {
	}

	void Add(duckdb::const_data_ptr_t data, duckdb::idx_t len) {
		if (expected_.empty() || len == 0) {
			return;
		}
		ctx_.Add(data, len);
	}

	// Finalizes the digest and throws duckdb::IOException on mismatch. `label`
	// identifies the run/file in the error message. No-op when the expected
	// md5 was empty, or when already finalized (safe to call from multiple
	// cleanup paths without double-throwing).
	//
	// The message intentionally avoids any of the transient-network wording
	// the ingest_ena_reads compute-orchestrator job's _TRANSIENT_ERROR_MARKERS
	// classifier scans for ("connection", "timed out", "timeout", "network",
	// "reset", "refused", "unreachable", "temporarily", "curl") -- an md5
	// mismatch is a data-integrity failure, not a retryable network blip, and
	// must not be misclassified as one.
	void VerifyOrThrow(const std::string &label) {
		if (expected_.empty() || finished_) {
			return;
		}
		finished_ = true;
		std::string actual = ctx_.FinishHex();
		if (actual != expected_) {
			throw duckdb::IOException(
			    "read_ena_sequences: md5 mismatch for '%s': ENA reported %s but downloaded bytes hash to %s", label,
			    expected_, actual);
		}
	}

private:
	std::string expected_;
	duckdb::MD5Context ctx_;
	bool finished_ = false;
};

} // namespace miint
