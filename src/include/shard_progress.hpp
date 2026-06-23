#pragma once

// Opt-in, human-readable, timestamped per-shard progress for the sharded
// aligners (align_minimap2_sharded / align_bowtie2_sharded). Gated by each
// function's `progress` named parameter (default false): when off, nothing is
// emitted and behavior is unchanged, so programmatic / library callers that
// never pass `progress := true` are unaffected. Lines go to stderr (like the
// `debug` facility in shard_debug.hpp) with a wall-clock timestamp.
//
// The `Format*` functions are PURE (no timing, no I/O) so they can be unit-
// tested in the standalone Catch2 binary, which links no libduckdb. This header
// is therefore intentionally duckdb-free: plain `uint64_t`/`int64_t`/`double`,
// no `idx_t` / `Value` / duckdb headers. `Emit` (timestamp + stderr write) is
// the only impure piece; it is exercised end-to-end by shaln's --verbose
// stderr-capture test, not the C++ unit test.

#include <chrono>
#include <cstdint>
#include <cstdio>
#include <ctime>
#include <mutex>
#include <string>

namespace duckdb {
namespace shard_progress {

// Group a non-negative integer with comma thousands separators: 48210 -> "48,210".
// Locale-independent (std::locale is avoided deliberately — it is not reliably
// available and would make the output environment-dependent).
inline std::string GroupThousands(uint64_t n) {
	const std::string digits = std::to_string(n);
	std::string out;
	const size_t len = digits.size();
	for (size_t i = 0; i < len; ++i) {
		if (i > 0 && (len - i) % 3 == 0) {
			out.push_back(',');
		}
		out.push_back(digits[i]);
	}
	return out;
}

// Body of a per-shard "started" line:
//   "shard 2/5 'shard_b': index loaded"           (n_reads < 0 -> count unknown)
//   "shard 2/5 'shard_b': index loaded, 50,000 reads"   (n_reads >= 0)
// n_reads is < 0 when the shard's read count isn't known at load time (the
// bowtie2 path streams reads; minimap2 pre-fetches them so it knows the count).
inline std::string FormatShardStart(uint64_t ordinal, uint64_t shard_count, const std::string &shard_name,
                                    int64_t n_reads) {
	std::string s =
	    "shard " + std::to_string(ordinal) + "/" + std::to_string(shard_count) + " '" + shard_name + "': index loaded";
	if (n_reads >= 0) {
		s += ", " + GroupThousands(static_cast<uint64_t>(n_reads)) + " reads";
	}
	return s;
}

// Body of a per-shard "done" line:
//   "shard 2/5 'shard_b': done - 50,000 reads, 48,210 alignments (2.34s)"
inline std::string FormatShardDone(uint64_t ordinal, uint64_t shard_count, const std::string &shard_name,
                                   uint64_t n_reads, uint64_t n_alignments, double elapsed_s) {
	char secs[32];
	std::snprintf(secs, sizeof(secs), "%.2fs", elapsed_s);
	return "shard " + std::to_string(ordinal) + "/" + std::to_string(shard_count) + " '" + shard_name + "': done - " +
	       GroupThousands(n_reads) + " reads, " + GroupThousands(n_alignments) + " alignments (" + secs + ")";
}

// Emit one formatted body line to stderr, prefixed with a local wall-clock
// timestamp (HH:MM:SS.t, tenths of a second) and an aligner `tag`
// (e.g. "minimap2" / "bowtie2") for greppability. Thread-safe: a process-wide
// mutex serializes the write so concurrent workers' lines never interleave.
// Only call when the caller's `progress` flag is set.
inline void Emit(const char *tag, const std::string &body) {
	using clock = std::chrono::system_clock;
	const auto now = clock::now();
	const std::time_t t = clock::to_time_t(now);
	const long long tenths = static_cast<long long>(
	    std::chrono::duration_cast<std::chrono::milliseconds>(now.time_since_epoch()).count() % 1000 / 100);
	std::tm tm_buf {};
#if defined(_WIN32)
	localtime_s(&tm_buf, &t);
#else
	localtime_r(&t, &tm_buf);
#endif
	char ts[16];
	std::snprintf(ts, sizeof(ts), "%02d:%02d:%02d.%01lld", tm_buf.tm_hour, tm_buf.tm_min, tm_buf.tm_sec, tenths);

	static std::mutex mtx;
	std::lock_guard<std::mutex> guard(mtx);
	std::fprintf(stderr, "%s [%s] %s\n", ts, tag, body.c_str());
	std::fflush(stderr);
}

} // namespace shard_progress
} // namespace duckdb
