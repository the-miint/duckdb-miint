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

/// Minimal RFC 8259 escape of a string for embedding inside a JSON string
/// literal. Escapes `"`, `\`, control characters (0x00–0x1F) using the
/// short forms (`\n` / `\r` / `\t`) where available and `\u00xx` otherwise.
/// The output does NOT include surrounding double quotes — caller emits those.
///
/// Use this any time a free-form string crosses into a hand-built JSON
/// message (e.g. `Session::Submit`'s batch envelope, the per-tool
/// `config_json` builders). Validated callers (e.g. `kKnownPresets`) don't
/// need it, but routing everything through this is cheaper than auditing
/// each builder for "is the source trusted today?".
std::string JsonEscape(const std::string &in);

/// One entry in the daemon's tool registry, as advertised on the init reply.
/// Populated by `Session::Initialize()` from the protocol-v3 `tools` field
/// (gpl-boundary commit 6b11337). Bind-time validators in table functions
/// can consult this to fail fast if the daemon doesn't speak the tool they
/// need or speaks the wrong schema version.
struct ToolEntry {
	std::string name;
	uint32_t schema_version = 0;
};

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
///      daemon replies `{success:true, protocol_version:3, tools:[...]}`,
///      throws on drift or error.
///   3. (Phase 4 will add `Submit()` for batches.)
///   4. `Shutdown()` sends `{shutdown:true}` and waits for the child to
///      exit. Idempotent. Safe to call after destruction would have
///      anyway reaped the child, but encouraged for the clean path.
///
/// Hard-codes `protocol_version: 3` per gpl-boundary `6b11337` (v0.2.0). A
/// version drift fails fast at `Initialize()` rather than producing garbled
/// batches down the line. The tool registry populated from the init reply
/// is reachable via `tools()` / `has_tool()` / `tool_schema_version()` so
/// bind-time validators can fail fast on unsupported tools.
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

	/// Tool registry as advertised by the daemon on the init reply. Empty
	/// before `Initialize()` returns or if the daemon omitted the field
	/// (older protocol versions would have already been rejected by the
	/// version check, so this is effectively non-empty post-Initialize on
	/// a v0.2.0+ daemon).
	const std::vector<ToolEntry> &tools() const {
		return tools_;
	}

	/// Linear search over `tools()` for `name`. O(n) but n ≤ ~10 in
	/// practice; not worth a hash map.
	bool has_tool(const std::string &name) const;

	/// Schema version of `name` in the daemon registry, or 0 if absent.
	/// Callers should compare against the schema version their code was
	/// written for and refuse to bind on mismatch.
	uint32_t tool_schema_version(const std::string &name) const;

	/// PID of the underlying daemon process. Useful for tests that verify
	/// connection-scoped daemon caching (same process across multiple
	/// `Submit` calls).
	pid_t daemon_pid() const {
		return child_.pid();
	}

private:
	/// Non-blocking drain of the daemon's stderr fd into `stderr_tail_`,
	/// truncating the front to the byte cap. Cheap no-op when nothing is
	/// pending. Called once per `Submit` (so the OS pipe buffer can't fill and
	/// wedge a verbose daemon on a blocked write) and again on the failure paths.
	void PumpStderr();

	ChildProcess child_;
	std::unique_ptr<LineReader> reader_;
	std::vector<ToolEntry> tools_;
	bool initialized_ = false;
	bool shut_down_ = false;
	int64_t next_batch_id_ = 1;

	/// Rolling tail of the daemon's stderr, accumulated across the session.
	/// Drained once per batch by `Submit` so a chatty tool (e.g. bowtie2 with
	/// quiet=false) can't fill the OS pipe buffer and deadlock the daemon on a
	/// blocked stderr write while we block waiting for its response. Capped to
	/// the last `kStderrTailCap` bytes; spliced into the exception message on
	/// the daemon-death / batch-failure paths.
	std::string stderr_tail_;
};

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
