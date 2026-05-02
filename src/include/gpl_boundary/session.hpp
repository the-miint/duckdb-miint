#pragma once

#include "gpl_boundary/process.hpp"
#include "gpl_boundary/shm.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace duckdb {
namespace miint {
namespace gpl_boundary {

/// One output segment produced by a batch — typically a single Arrow IPC
/// stream byte buffer (e.g., the FastTree SOA tree). The `OutputShmRegion`
/// owns the mmap; `bytes()` is a non-owning view over it. Keep the
/// `SubmitResult` alive until you've finished reading the bytes.
struct SubmitOutput {
	std::string label;   // matches gpl-boundary's `ShmOutput.label`
	OutputShmRegion shm; // owns mmap + auto-unlink on destroy

	const void *bytes() const {
		return shm.data();
	}
	std::size_t size_bytes() const {
		return shm.size_bytes();
	}
};

/// Response from a single `Session::Submit` call.
struct SubmitResult {
	int64_t batch_id = 0;
	uint32_t schema_version = 0;       // gpl-boundary per-tool schema version
	std::vector<SubmitOutput> outputs; // one entry per shm_outputs[]
	std::string result_json;           // raw JSON for the optional `result` field; empty if absent
};

/// Connection-scoped handle to a `gpl-boundary` daemon process.
///
/// Lifecycle (Phase 1.4):
///   1. Construct with a `ChildProcess` that's already spawned (caller's
///      responsibility — `ChildProcess(...)` does the fork+exec).
///   2. `Initialize()` sends the `{init: {...}}` line, blocks until the
///      daemon replies `{success:true, protocol_version:2}`, throws on
///      drift or error.
///   3. (Phase 4 will add `Submit()` for batches.)
///   4. `Shutdown()` sends `{shutdown:true}` and waits for the child to
///      exit. Idempotent. Safe to call after destruction would have
///      anyway reaped the child, but encouraged for the clean path.
///
/// Hard-codes `protocol_version: 2` per gpl-boundary `19306f6`. A version
/// drift fails fast at `Initialize()` rather than producing garbled batches
/// down the line.
class Session {
public:
	explicit Session(ChildProcess child);
	~Session();

	Session(const Session &) = delete;
	Session &operator=(const Session &) = delete;
	Session(Session &&) noexcept;
	Session &operator=(Session &&) noexcept;

	/// Send the init line and validate the daemon's response.
	/// Throws `std::runtime_error` on protocol mismatch, daemon error,
	/// premature EOF, or read/write failure.
	void Initialize();

	/// Send `{shutdown:true}` (best-effort) and wait for the child to exit.
	/// Idempotent: subsequent calls are no-ops.
	void Shutdown();

	/// Submit one batch to the daemon and block until the response arrives.
	///
	/// - `tool` — gpl-boundary tool name (e.g. `"fasttree"`).
	/// - `config_json` — already-serialized JSON object string for the
	///   tool's config block (e.g. `R"({"seed":42,"threads":1})"`).
	///   Caller is responsible for validity at the schema level
	///   (Phase 5's `phylogeny_fasttree` Bind builds + validates this).
	/// - `input_bytes` / `input_size` — Arrow IPC stream bytes of the input
	///   batch (typically produced by `EncodeIpcStream`).
	///
	/// Sends the wire-protocol-2 batch line including `shm_input_size`
	/// (mandatory at gpl-boundary `19306f6`+), reads the response, opens each
	/// output shm segment with the explicit size from `ShmOutput.size`, and
	/// returns the bundle. Input shm is unlinked before this function returns.
	///
	/// Throws `std::runtime_error` if `Initialize()` hasn't run yet, the
	/// daemon returns `{success:false}`, the response is malformed, or any
	/// I/O fails. The output segments are unlinked automatically when the
	/// returned `SubmitResult` is destroyed.
	SubmitResult Submit(const std::string &tool, const std::string &config_json, const void *input_bytes,
	                    std::size_t input_size);

	/// True iff Initialize() has completed successfully.
	bool initialized() const {
		return initialized_;
	}

	/// PID of the underlying daemon process. Useful for tests that verify
	/// connection-scoped daemon caching (same process across multiple
	/// `Submit` calls).
	pid_t daemon_pid() const {
		return child_.pid();
	}

private:
	ChildProcess child_;
	std::unique_ptr<LineReader> reader_;
	bool initialized_ = false;
	bool shut_down_ = false;
	int64_t next_batch_id_ = 1;
};

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
