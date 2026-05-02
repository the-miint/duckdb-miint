#pragma once

#include <string>
#include <sys/types.h>
#include <vector>

namespace duckdb {
namespace miint {
namespace gpl_boundary {

/// Locate the `gpl-boundary` binary on PATH using `which`.
///
/// Returns the absolute path on success, or empty string if not found. Never
/// throws. Used by `phylogeny_fasttree_available()` and lazy session creation.
std::string FindGplBoundary();

/// Generic helper: locate any executable on PATH using `which`. Returns the
/// absolute path on success, empty string otherwise. Same fork/exec/pipe pattern
/// as `Bowtie2Aligner::find_executable` (`src/Bowtie2Aligner.cpp:26-81`); shares
/// the no-shell-injection guarantee. Mostly here so unit tests can exercise the
/// path-discovery plumbing against a known-present binary like `bash`.
std::string FindExecutableInPath(const std::string &name);

/// Long-lived child process with bidirectional pipes on stdin/stdout/stderr.
/// Used to host the `gpl-boundary` daemon for the entire lifetime of a DuckDB
/// connection.
///
/// On construction:
///   - pipes are created for stdin/stdout/stderr,
///   - the child is forked,
///   - argv[0] is execv'd (no shell, no PATH search — pass an absolute path,
///     or a path resolved by `FindExecutableInPath`),
///   - on Linux, the child sets `PR_SET_PDEATHSIG = SIGTERM` so it is reaped
///     if the parent dies abruptly.
///
/// On destruction:
///   - any remaining pipe fds are closed,
///   - if the child is still running, SIGTERM is sent and the process is
///     reaped (with a short SIGKILL fallback after a grace period — mirrors
///     `Bowtie2Aligner` (`src/Bowtie2Aligner.cpp:701-775`) but tighter
///     since the gpl-boundary daemon has its own graceful-shutdown protocol
///     that should already have run before destruction).
///
/// Move-only.
class ChildProcess {
public:
	/// argv[0] must be an absolute path (or a path that resolves under the
	/// current working directory). Throws `std::runtime_error` if pipes,
	/// fork, or exec fail.
	explicit ChildProcess(const std::vector<std::string> &argv);
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

private:
	pid_t pid_ = -1;
	int stdin_fd_ = -1;
	int stdout_fd_ = -1;
	int stderr_fd_ = -1;
	bool waited_ = false;
	int wait_status_ = 0;
};

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

private:
	int fd_;
	std::string buffer_; // unread bytes (stuff already consumed from the fd
	                     // but not yet returned via ReadLine)
	bool eof_ = false;
};

/// Write `line` to `fd` followed by a single `\n`, flushing immediately.
/// Throws `std::runtime_error` on partial-write or write error. Used to send
/// NDJSON lines to the gpl-boundary daemon's stdin.
void WriteLine(int fd, const std::string &line);

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
