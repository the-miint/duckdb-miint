#pragma once

#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

namespace duckdb {
namespace miint {
namespace gpl_boundary {

/// Locate the `gpl-boundary` binary. Lookup order:
///   1. `$MIINT_GPL_BOUNDARY_PATH` if set (explicit override; not validated to
///      exist — the caller can probe and surface a clearer error than "not
///      found in any location").
///   2. miint's install cache (`MiintGplBoundaryCacheBinary()`), where
///      `install_gpl_boundary()` deposits binaries.
///   3. `which gpl-boundary` (whatever's on PATH).
///
/// Returns the absolute path on success, or empty string if not found. Never
/// throws. Used by `phylogeny_fasttree_available()` and lazy session creation.
std::string FindGplBoundary();

/// Returns the directory miint installs gpl-boundary into via
/// `install_gpl_boundary()`. Resolves to (in order): `$XDG_CACHE_HOME/miint/bin`,
/// `$HOME/.cache/miint/bin`. Returns empty string if neither HOME nor
/// XDG_CACHE_HOME is set (extremely unlikely; signals "no install location"
/// to callers).
std::string MiintGplBoundaryCacheDir();

/// Returns `MiintGplBoundaryCacheDir() + "/gpl-boundary"`. Empty if the dir
/// itself is unresolvable.
std::string MiintGplBoundaryCacheBinary();

/// Generic helper: locate any executable on PATH using `which`. Returns the
/// absolute path on success, empty string otherwise. Uses a fork/exec/pipe
/// pattern with no-shell-injection guarantees. Mostly here so unit tests can
/// exercise the path-discovery plumbing against a known-present binary like
/// `bash`.
std::string FindExecutableInPath(const std::string &name);

/// How `ChildProcess` starts its child. The defaults are what gpl-boundary
/// relies on; a host that cleans up after itself when its caller goes away
/// (skb-arrow) wants `parent_death_signal = false` and `new_session = true`.
struct ChildProcessOptions {
	/// Linux only: the child sets `PR_SET_PDEATHSIG = SIGTERM`. Per prctl(2)
	/// the signal fires when the *thread* that forked exits, not the process,
	/// and DuckDB often runs table-function init on a worker thread.
	bool parent_death_signal = true;
	/// The child calls `setsid()`, so it leads its own process group and
	/// `KillGroup()` reaches everything it starts. It also leaves the
	/// terminal's process group, so a terminal Ctrl-C reaches only DuckDB.
	bool new_session = false;
	/// Set in the child's environment, replacing any inherited value.
	std::vector<std::pair<std::string, std::string>> set_env;
	/// Removed from the child's environment.
	std::vector<std::string> unset_env;
};

/// Long-lived child process with bidirectional pipes on stdin/stdout/stderr.
/// Used to host the `gpl-boundary` daemon for the entire lifetime of a DuckDB
/// connection.
///
/// On construction:
///   - close-on-exec pipes are created for stdin/stdout/stderr, so no other
///     child (spawned before or after) holds this child's pipe ends,
///   - the child is forked,
///   - argv[0] is execv'd (no shell, no PATH search — pass an absolute path,
///     or a path resolved by `FindExecutableInPath`),
///   - `options` applies (see `ChildProcessOptions`).
///
/// On destruction:
///   - any remaining pipe fds are closed,
///   - if the child is still running, SIGTERM is sent (with `new_session`,
///     SIGKILL to its whole process group) and the process is reaped (with a
///     short SIGKILL fallback after a grace period). The gpl-boundary daemon
///     has its own graceful-shutdown protocol that should already have run
///     before destruction.
///
/// Move-only.
class ChildProcess {
public:
	/// argv[0] must be an absolute path (or a path that resolves under the
	/// current working directory). Throws `std::runtime_error` if pipes,
	/// fork, or exec fail.
	explicit ChildProcess(const std::vector<std::string> &argv, const ChildProcessOptions &options = {});
	~ChildProcess();

	ChildProcess(const ChildProcess &) = delete;
	ChildProcess &operator=(const ChildProcess &) = delete;
	ChildProcess(ChildProcess &&) noexcept;
	ChildProcess &operator=(ChildProcess &&) noexcept;

	/// File descriptors. The caller may read from / write to these directly,
	/// but should NOT close them — `ChildProcess` owns the lifetime.
	int stdin_fd() const {
		return stdin_fd_;
	}
	int stdout_fd() const {
		return stdout_fd_;
	}
	int stderr_fd() const {
		return stderr_fd_;
	}

	/// Child PID. Useful for diagnostics; not for direct kill().
	pid_t pid() const {
		return pid_;
	}

	/// Close stdin (signals EOF to the child). Idempotent.
	void CloseStdin();

	/// Wait for the child to exit, returning the raw waitpid status (use
	/// WIFEXITED / WEXITSTATUS / etc.). Idempotent: subsequent calls return
	/// the cached status.
	int Wait();

	/// Reap the child with a bounded grace period and return a human-readable,
	/// signal-decoded description of how it terminated:
	///   - "killed by signal 15 (Terminated)"  — e.g. PR_SET_PDEATHSIG SIGTERM
	///     firing, or the OOM killer's SIGKILL (9), or a crash (SIGSEGV/SIGABRT)
	///   - "exited with code 7"
	///   - "still running after 200ms (no exit yet)" — closed a pipe yet alive
	///   - "already reaped or unwaitable (...)" — ECHILD: something else waited
	/// Non-throwing. Idempotent with Wait(): when the child is reaped here the
	/// status is cached, so the destructor and Wait() won't block again. Used on
	/// the daemon-death diagnostic paths in Session::Submit so a bare EPIPE
	/// ("Broken pipe") or stdout-EOF gains the actual cause.
	std::string ReapAndDescribe(int grace_ms = 200);

	/// SIGKILL the child's whole process group: the child and everything it
	/// started. Requires `new_session` (throws `std::logic_error` otherwise).
	/// Does nothing once the child is reaped, because its id may then name an
	/// unrelated group. A group already gone (ESRCH, or EPERM on macOS when
	/// only zombies remain) counts as killed. Reap with `Wait()` afterwards.
	void KillGroup();

private:
	pid_t pid_ = -1;
	int stdin_fd_ = -1;
	int stdout_fd_ = -1;
	int stderr_fd_ = -1;
	bool waited_ = false;
	int wait_status_ = 0;
	bool new_session_ = false;
};

/// Result of `RunCapture`.
struct RunResult {
	bool spawned = false;         // true iff fork+exec succeeded
	bool exited_normally = false; // true iff child returned via exit()
	int exit_code = -1;           // valid iff exited_normally
	int term_signal = 0;          // non-zero iff WIFSIGNALED — child was killed
	std::string stdout_text;
	std::string stderr_text;
	std::string spawn_error; // populated iff !spawned
};

/// Run argv[0] to completion with stdin closed, capturing stdout and stderr
/// separately. No shell, so no quoting. Never throws: a spawn failure comes
/// back in `spawn_error`.
RunResult RunCapture(const std::vector<std::string> &argv, const ChildProcessOptions &options = {});

/// Buffered line reader over a raw fd. One `ReadLine` call returns one line
/// (without the trailing `\n`).
///
/// EOF semantics: when the fd reaches EOF, any final unterminated line is
/// returned via one last `ReadLine(...) -> true` call. The next call returns
/// `false` (with `out` cleared). Callers that need to distinguish between
/// "got a real terminated line" and "got the trailing partial line" should
/// inspect what they're doing on the protocol side; for NDJSON this case
/// can only happen if the daemon crashed mid-line.
///
/// Not thread-safe; one reader per fd.
class LineReader {
public:
	explicit LineReader(int fd);

	LineReader(const LineReader &) = delete;
	LineReader &operator=(const LineReader &) = delete;

	/// Returns true and writes the next line into `out` (excluding the `\n`).
	/// Returns false on EOF (with `out` cleared). Throws on read errors.
	bool ReadLine(std::string &out);

	/// True when `ReadLine` would return without reading: a whole line is
	/// buffered, or the fd has reached EOF.
	bool HasLine() const;

	/// True when a whole, newline-terminated line is buffered: unlike
	/// `HasLine`, false at EOF after a partial last line.
	bool HasCompleteLine() const;

	/// One read() into the buffer: for a caller that polls the fd (to stay
	/// interruptible) and calls `ReadLine` only once `HasLine()`. Blocks only
	/// if the fd has nothing to read. Throws on read errors.
	void FillOnce();

private:
	int fd_;
	std::string buffer_; // unread bytes (stuff already consumed from the fd
	                     // but not yet returned via ReadLine)
	bool eof_ = false;
};

/// Trim `s` from the front to at most `cap` bytes, prefixing a marker so the
/// reader knows output was dropped. Keeps the tail — that's where the
/// failing call's output lives.
void TruncateHead(std::string &s, std::size_t cap);

/// Write `line` to `fd` followed by a single `\n`, flushing immediately.
/// Throws `std::runtime_error` on partial-write or write error. Used to send
/// NDJSON lines to the gpl-boundary daemon's stdin.
void WriteLine(int fd, const std::string &line);

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
