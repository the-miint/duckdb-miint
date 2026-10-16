#pragma once

// Client for skb-arrow, a process host that runs curated scikit-bio
// operations over Arrow (https://github.com/the-miint/skb-arrow). The contract
// is skb-arrow's protocol spec at v0.1.0: docs/protocol.md (messages),
// docs/transport.md (segment files), docs/errors.md (error kinds, warnings).
//
// DuckDB-free, so the Catch2 `tests` binary (which has no libduckdb) can
// drive it against fake and real hosts.

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include "gpl_boundary/process.hpp"

struct ArrowSchema;
struct ArrowArray;

namespace duckdb {
namespace miint {
namespace skb_arrow {

/// The wire protocol this client speaks.
constexpr int kProtocolVersion = 1;

/// Target bytes of batches per input segment: skb-arrow's own default.
constexpr std::size_t kSegmentBytes = std::size_t(256) << 20;

/// A failure in skb-arrow's vocabulary: `host_incompatible`, `invalid_input`,
/// `invalid_param`, `unsupported`, `resource` or `internal`. The host reports
/// most of them; the client uses the same kinds for its own faults: a session
/// directory or segment it can't create or write is `resource`, a malformed
/// response is `host_incompatible`, and an unknown kind or a host that dies is
/// `internal`.
class Error : public std::runtime_error {
public:
	Error(std::string kind, const std::string &message, std::string traceback = {});
	const std::string &kind() const {
		return kind_;
	}
	/// The host's Python traceback. Only an `internal` error from the host has one.
	const std::string &traceback() const {
		return traceback_;
	}

private:
	std::string kind_;
	std::string traceback_;
};

/// The capabilities miint's SQL functions call, each at the `schema_version`
/// they were built against.
std::vector<std::pair<std::string, int64_t>> MiintCapabilities();

/// Called while waiting on the host. Throw to abandon the call: the host and
/// everything it started are killed, and the exception propagates.
using InterruptCheck = std::function<void()>;

/// Where `install_skb_arrow()` puts the binary: beside gpl-boundary in miint's
/// cache (`MiintGplBoundaryCacheDir()`). Empty when there is no cache dir.
std::string MiintSkbArrowCacheBinary();

/// Locate skb-arrow: `$MIINT_SKB_ARROW_PATH` (not validated), then miint's
/// cache, then PATH. Empty when none is found.
std::string FindSkbArrow();

/// Where session directories go: `$MIINT_SKB_ARROW_SEGMENT_DIR`; else, on
/// Linux, `/dev/shm` when it is a writable directory; else `$TMPDIR`; else
/// `/tmp`. Segments are memory: tmpfs keeps them off disk.
std::string SegmentDirParent();

/// A private (0700), empty directory for one host's segments, made with
/// mkdtemp under `parent`. Removed, with whatever is left in it, on
/// destruction.
class SessionDir {
public:
	explicit SessionDir(const std::string &parent);
	~SessionDir();
	SessionDir(const SessionDir &) = delete;
	SessionDir &operator=(const SessionDir &) = delete;

	const std::string &path() const {
		return path_;
	}
	/// Unlink every entry, then the directory. Idempotent: what is already
	/// gone (ENOENT) is done.
	void Remove() noexcept;

private:
	std::string path_;
};

/// Writes one input table as segment files `<prefix>-0`, `<prefix>-1`, … in
/// `dir`. Each is an uncompressed Arrow IPC stream of the schema and whole
/// batches, created exclusively (0600) and written with write(2), so a full
/// filesystem fails with ENOSPC rather than SIGBUS. A new segment starts once
/// the current one holds `cap` bytes of batches. With no batches, `Finish`
/// writes one schema-only segment: skb-arrow's empty table.
class SegmentWriter {
public:
	/// Copies `schema`; never releases it.
	SegmentWriter(std::string dir, std::string prefix, const ArrowSchema &schema, std::size_t cap = kSegmentBytes);
	/// Unlinks every segment written unless `Finish` succeeded.
	~SegmentWriter();
	SegmentWriter(const SegmentWriter &) = delete;
	SegmentWriter &operator=(const SegmentWriter &) = delete;

	/// Reads `batch`; never releases it.
	void Append(const ArrowArray &batch);
	/// Ends the last segment. Returns the segment names, in order.
	std::vector<std::string> Finish();

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

/// A capability's params: JSON scalars and arrays of strings. The host checks
/// names and values; an unset param takes the host's default.
class Params {
public:
	using Value = std::variant<int64_t, double, bool, std::string, std::vector<std::string>>;

	void SetInt(const std::string &name, int64_t value);
	/// Throws `Error("invalid_param")` unless `value` is finite: NaN and
	/// Infinity are not JSON.
	void SetNumber(const std::string &name, double value);
	void SetBool(const std::string &name, bool value);
	void SetString(const std::string &name, std::string value);
	void SetStrings(const std::string &name, std::vector<std::string> value);

	const std::vector<std::pair<std::string, Value>> &entries() const {
		return entries_;
	}

private:
	std::vector<std::pair<std::string, Value>> entries_;
};

/// One entry of a response's `warnings`.
struct Warning {
	std::string category;
	std::string message;
	int64_t count = 0;
};

/// A `result` response: the output table's segment names, in order.
struct CallResult {
	std::vector<std::string> output;
	std::vector<Warning> warnings;
};

/// One output segment, mapped read-only. Its file is already unlinked; the
/// bytes stay valid until this is destroyed.
class MappedSegment {
public:
	/// Takes ownership of an mmap(2) mapping.
	MappedSegment(void *data, std::size_t size);
	~MappedSegment();
	MappedSegment(MappedSegment &&other) noexcept;
	MappedSegment &operator=(MappedSegment &&other) noexcept;
	MappedSegment(const MappedSegment &) = delete;
	MappedSegment &operator=(const MappedSegment &) = delete;

	const uint8_t *data() const {
		return static_cast<const uint8_t *>(data_);
	}
	std::size_t size() const {
		return size_;
	}

private:
	void *data_ = nullptr;
	std::size_t size_ = 0;
};

struct HostOptions {
	/// skb-arrow's path, e.g. from `FindSkbArrow()`.
	std::string binary;
	/// Where the session directory goes, e.g. `SegmentDirParent()`.
	std::string segment_parent;
	/// Above 0, caps the host's threads: NUMBA_NUM_THREADS, OMP_NUM_THREADS,
	/// OPENBLAS_NUM_THREADS and VECLIB_MAXIMUM_THREADS.
	int threads = 0;
	/// Capabilities the caller uses, each with the `schema_version` it was
	/// built against. The handshake refuses a host without an exact match.
	std::vector<std::pair<std::string, int64_t>> required;
	/// How long `Shutdown` waits for the host to exit before killing it.
	int shutdown_grace_ms = 5000;
};

/// One skb-arrow host, in its own process group, with its own session
/// directory. Hosts are stateless, so any host serves any call. Not
/// thread-safe.
///
/// Order, per the contract: the handshake (in the constructor), then input
/// segments (a `SegmentWriter` in `dir()`), then `Call`, then `MapOutput` for
/// each output, then `Shutdown`. An error *response* leaves the host serving.
/// Any other failure (the host died or broke the protocol, or `check`
/// threw) leaves it gone: killed, reaped, and its directory removed.
class Host {
public:
	/// Creates the session directory, spawns `skb-arrow --segment-dir DIR`,
	/// and completes the `init` handshake. Throws `Error`.
	Host(const HostOptions &options, const InterruptCheck &check);
	/// Kills the host unless it was shut down.
	~Host();
	Host(const Host &) = delete;
	Host &operator=(const Host &) = delete;

	const std::string &host_version() const {
		return host_version_;
	}
	const std::string &dir() const {
		return dir_.path();
	}

	/// Sends one `call` and waits for its response, calling `check` about
	/// every 100 ms. `input` maps each input table to its segments, which
	/// the host unlinks. Throws `Error` for an error response, or for a host
	/// that dies or breaks the protocol.
	CallResult Call(const std::string &capability, const Params &params,
	                const std::vector<std::pair<std::string, std::vector<std::string>>> &input,
	                const InterruptCheck &check);

	/// Map output segment `name`, then unlink it.
	MappedSegment MapOutput(const std::string &name);

	/// Close the host's stdin, so it cleans up and exits; kill it if it
	/// hasn't within `shutdown_grace_ms`. Then remove the session directory.
	void Shutdown();

	/// SIGKILL the host and everything it started, reap it, and remove the
	/// session directory. Idempotent.
	void Kill() noexcept;

private:
	// Send one request line and wait for its response line.
	std::string Exchange(const std::string &request, const InterruptCheck &check, const char *awaiting);
	[[noreturn]] void ThrowHostDeath(const char *awaiting);

	HostOptions options_;
	// Declared before `child_`, so it is destroyed after: the host is gone
	// before its directory is removed.
	SessionDir dir_;
	gpl_boundary::ChildProcess child_;
	gpl_boundary::LineReader reader_;
	std::string host_version_;
	int64_t next_id_ = 1;
	bool done_ = false;
};

} // namespace skb_arrow
} // namespace miint
} // namespace duckdb
