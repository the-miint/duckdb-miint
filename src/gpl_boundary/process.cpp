#include "gpl_boundary/process.hpp"

#include <array>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <signal.h>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

#if defined(__linux__)
#include <sys/prctl.h>
#endif

namespace duckdb {
namespace miint {
namespace gpl_boundary {

namespace {
// RAII guard for an fd; only relevant in error paths during construction.
struct FdCloser {
	int fd;
	~FdCloser() {
		if (fd >= 0) {
			::close(fd);
		}
	}
	void release() {
		fd = -1;
	}
};

// Make a non-cloexec pipe. We deliberately leave fds inheritable across
// fork+exec because we wire them to the child's stdio via dup2.
void make_pipe(int fds[2]) {
	if (::pipe(fds) == -1) {
		throw std::runtime_error("gpl_boundary: pipe() failed: " + std::string(::strerror(errno)));
	}
}
} // namespace

std::string FindExecutableInPath(const std::string &name) {
	int pipefd[2];
	if (::pipe(pipefd) == -1) {
		return "";
	}

	pid_t pid = ::fork();
	if (pid == -1) {
		::close(pipefd[0]);
		::close(pipefd[1]);
		return "";
	}

	if (pid == 0) {
		// Child
		::close(pipefd[0]);
		::dup2(pipefd[1], STDOUT_FILENO);
		::close(pipefd[1]);

		int devnull = ::open("/dev/null", O_WRONLY);
		if (devnull != -1) {
			::dup2(devnull, STDERR_FILENO);
			::close(devnull);
		}

		::execlp("which", "which", name.c_str(), nullptr);
		::_exit(127);
	}

	::close(pipefd[1]);

	std::array<char, 256> buffer;
	std::string result;
	for (;;) {
		ssize_t n = ::read(pipefd[0], buffer.data(), buffer.size());
		if (n > 0) {
			result.append(buffer.data(), static_cast<size_t>(n));
			continue;
		}
		if (n == 0) {
			break; // EOF
		}
		if (errno == EINTR) {
			continue;
		}
		break; // any other error: bail with whatever we have
	}
	::close(pipefd[0]);

	int status = 0;
	::waitpid(pid, &status, 0);

	if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		return "";
	}

	while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) {
		result.pop_back();
	}
	return result;
}

std::string MiintGplBoundaryCacheDir() {
	if (const char *xdg = ::getenv("XDG_CACHE_HOME"); xdg && xdg[0] != '\0') {
		return std::string(xdg) + "/miint/bin";
	}
	if (const char *home = ::getenv("HOME"); home && home[0] != '\0') {
		return std::string(home) + "/.cache/miint/bin";
	}
	return {};
}

std::string MiintGplBoundaryCacheBinary() {
	const std::string dir = MiintGplBoundaryCacheDir();
	if (dir.empty()) {
		return {};
	}
	return dir + "/gpl-boundary";
}

std::string FindGplBoundary() {
	// 1. Explicit override.
	if (const char *override_path = ::getenv("MIINT_GPL_BOUNDARY_PATH"); override_path && override_path[0] != '\0') {
		return std::string(override_path);
	}
	// 2. miint's install cache (where install_gpl_boundary() deposits).
	std::string cached = MiintGplBoundaryCacheBinary();
	if (!cached.empty() && ::access(cached.c_str(), X_OK) == 0) {
		return cached;
	}
	// 3. PATH.
	return FindExecutableInPath("gpl-boundary");
}

// =============================================================================
// ChildProcess
// =============================================================================

ChildProcess::ChildProcess(const std::vector<std::string> &argv) {
	if (argv.empty()) {
		throw std::runtime_error("gpl_boundary: ChildProcess requires non-empty argv");
	}

	int in_pipe[2] = {-1, -1};
	int out_pipe[2] = {-1, -1};
	int err_pipe[2] = {-1, -1};
	int exec_status_pipe[2] = {-1, -1}; // self-pipe; child writes errno on exec failure

	// Track all fds for cleanup on error.
	auto close_all = [&]() {
		for (int *p : {in_pipe, out_pipe, err_pipe, exec_status_pipe}) {
			if (p[0] >= 0) {
				::close(p[0]);
				p[0] = -1;
			}
			if (p[1] >= 0) {
				::close(p[1]);
				p[1] = -1;
			}
		}
	};

	try {
		make_pipe(in_pipe);
		make_pipe(out_pipe);
		make_pipe(err_pipe);
		make_pipe(exec_status_pipe);
	} catch (...) {
		close_all();
		throw;
	}

	// CLOEXEC on the exec-status pipe: parent reads "" (EOF) when execv
	// succeeds (the kernel closes both fds across the exec), or reads the
	// raw errno value when the child wrote one before _exit. This avoids
	// the race where waitpid(WNOHANG) returns 0 because the child hasn't
	// yet had time to fail exec.
	if (::fcntl(exec_status_pipe[1], F_SETFD, FD_CLOEXEC) == -1) {
		const int err = errno;
		close_all();
		throw std::runtime_error("gpl_boundary: fcntl(FD_CLOEXEC) failed: " + std::string(::strerror(err)));
	}

	// Build argv as char* array with terminating nullptr.
	std::vector<std::string> owned = argv;
	std::vector<char *> c_argv;
	c_argv.reserve(owned.size() + 1);
	for (auto &s : owned) {
		c_argv.push_back(s.data());
	}
	c_argv.push_back(nullptr);

	pid_t pid = ::fork();
	if (pid == -1) {
		const int err = errno;
		close_all();
		throw std::runtime_error("gpl_boundary: fork() failed: " + std::string(::strerror(err)));
	}

	if (pid == 0) {
		// Child. Async-signal-safe path only — execv soon clears the slate.
		::dup2(in_pipe[0], STDIN_FILENO);
		::dup2(out_pipe[1], STDOUT_FILENO);
		::dup2(err_pipe[1], STDERR_FILENO);

		// Close every fd we no longer need; the parent must see EOF on its
		// read ends after the child's natural exit. exec_status_pipe[1] keeps
		// CLOEXEC on it and dies when execv succeeds.
		::close(in_pipe[0]);
		::close(in_pipe[1]);
		::close(out_pipe[0]);
		::close(out_pipe[1]);
		::close(err_pipe[0]);
		::close(err_pipe[1]);
		::close(exec_status_pipe[0]);

#if defined(__linux__)
		::prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif

		::execv(c_argv[0], c_argv.data());

		// exec failed — write errno into the status pipe so the parent knows
		// why. write/_exit are async-signal-safe.
		const int err = errno;
		(void)::write(exec_status_pipe[1], &err, sizeof(err));
		::_exit(127);
	}

	// Parent. Close the child-side ends so EOF works correctly.
	::close(in_pipe[0]);
	::close(out_pipe[1]);
	::close(err_pipe[1]);
	::close(exec_status_pipe[1]);

	pid_ = pid;
	stdin_fd_ = in_pipe[1];
	stdout_fd_ = out_pipe[0];
	stderr_fd_ = err_pipe[0];

	// Block until the child either execs (CLOEXEC closes the pipe → EOF) or
	// fails exec (writes errno). Worst case we sleep until execv runs in the
	// child, which is microseconds.
	//
	// Short-read loop: write(2) of 4 bytes is atomic on Linux pipes so this
	// almost certainly returns sizeof(int) in one shot, but POSIX permits a
	// short read after a signal; loop defensively rather than misinterpret a
	// partial errno.
	int exec_errno = 0;
	ssize_t total = 0;
	bool got_eof = false;
	while (total < static_cast<ssize_t>(sizeof(exec_errno))) {
		ssize_t n = ::read(exec_status_pipe[0], reinterpret_cast<char *>(&exec_errno) + total,
		                   sizeof(exec_errno) - static_cast<size_t>(total));
		if (n == 0) {
			got_eof = true;
			break;
		}
		if (n < 0) {
			if (errno == EINTR) {
				continue;
			}
			break;
		}
		total += n;
	}
	::close(exec_status_pipe[0]);
	if (!got_eof && total > 0) {
		// exec failed
		const std::string program = c_argv[0];
		// Reap so we don't leak a zombie.
		int status = 0;
		::waitpid(pid_, &status, 0);
		waited_ = true;
		wait_status_ = status;
		if (stdin_fd_ >= 0) {
			::close(stdin_fd_);
			stdin_fd_ = -1;
		}
		if (stdout_fd_ >= 0) {
			::close(stdout_fd_);
			stdout_fd_ = -1;
		}
		if (stderr_fd_ >= 0) {
			::close(stderr_fd_);
			stderr_fd_ = -1;
		}
		throw std::runtime_error("gpl_boundary: failed to exec '" + program +
		                         "': " + std::string(::strerror(exec_errno)));
	}
	// n == 0 (EOF): execv succeeded.
}

ChildProcess::~ChildProcess() {
	// Order matters: SIGTERM the daemon FIRST so its signal handler can run
	// while its stdin pipe is still open. Closing stdin first would cause it
	// to wake from `read()` with EOF and possibly take a different (uglier)
	// shutdown path — for tools that hold POSIX shm handles (Phase 3+) we
	// want them to receive SIGTERM and run their unlink-on-shutdown handler,
	// not exit through an "unexpected EOF" branch.
	if (pid_ > 0 && !waited_) {
		::kill(pid_, SIGTERM);
	}
	if (stdin_fd_ >= 0) {
		::close(stdin_fd_);
	}
	if (stdout_fd_ >= 0) {
		::close(stdout_fd_);
	}
	if (stderr_fd_ >= 0) {
		::close(stderr_fd_);
	}
	if (pid_ > 0 && !waited_) {
		// 30 × 10 ms grace, then SIGKILL.
		for (int i = 0; i < 30; ++i) {
			int status = 0;
			const pid_t r = ::waitpid(pid_, &status, WNOHANG);
			if (r == pid_) {
				return;
			}
			if (r == -1) {
				// Already gone (ECHILD) or error — give up.
				return;
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(10));
		}
		::kill(pid_, SIGKILL);
		int status = 0;
		::waitpid(pid_, &status, 0);
	}
}

ChildProcess::ChildProcess(ChildProcess &&other) noexcept
    : pid_(other.pid_), stdin_fd_(other.stdin_fd_), stdout_fd_(other.stdout_fd_), stderr_fd_(other.stderr_fd_),
      waited_(other.waited_), wait_status_(other.wait_status_) {
	other.pid_ = -1;
	other.stdin_fd_ = -1;
	other.stdout_fd_ = -1;
	other.stderr_fd_ = -1;
	other.waited_ = true;
}

ChildProcess &ChildProcess::operator=(ChildProcess &&other) noexcept {
	if (this != &other) {
		this->~ChildProcess();
		new (this) ChildProcess(std::move(other));
	}
	return *this;
}

void ChildProcess::CloseStdin() {
	if (stdin_fd_ >= 0) {
		::close(stdin_fd_);
		stdin_fd_ = -1;
	}
}

int ChildProcess::Wait() {
	if (waited_) {
		return wait_status_;
	}
	int status = 0;
	::waitpid(pid_, &status, 0);
	wait_status_ = status;
	waited_ = true;
	return status;
}

// =============================================================================
// LineReader
// =============================================================================

LineReader::LineReader(int fd) : fd_(fd) {
}

bool LineReader::ReadLine(std::string &out) {
	out.clear();
	for (;;) {
		// Pull the next full line out of the buffer if possible.
		const auto nl = buffer_.find('\n');
		if (nl != std::string::npos) {
			out.assign(buffer_, 0, nl);
			buffer_.erase(0, nl + 1);
			return true;
		}
		if (eof_) {
			// No newline left, no more bytes coming. If there's anything in
			// the buffer, surface it as the final unterminated line.
			if (!buffer_.empty()) {
				out = std::move(buffer_);
				buffer_.clear();
				return true;
			}
			return false;
		}
		// Need more bytes.
		char chunk[4096];
		ssize_t n;
		do {
			n = ::read(fd_, chunk, sizeof(chunk));
		} while (n == -1 && errno == EINTR);
		if (n < 0) {
			throw std::runtime_error("gpl_boundary: LineReader::ReadLine read() failed: " +
			                         std::string(::strerror(errno)));
		}
		if (n == 0) {
			eof_ = true;
			continue;
		}
		buffer_.append(chunk, static_cast<size_t>(n));
	}
}

// =============================================================================
// WriteLine
// =============================================================================

void WriteLine(int fd, const std::string &line) {
	// One contiguous payload (line + \n) so we have at most one short-write
	// recovery loop. write() can return short on regular files; on pipes it
	// usually returns full but EINTR is possible.
	std::string payload = line;
	payload.push_back('\n');
	const char *data = payload.data();
	size_t remaining = payload.size();
	while (remaining > 0) {
		ssize_t n;
		do {
			n = ::write(fd, data, remaining);
		} while (n == -1 && errno == EINTR);
		if (n < 0) {
			throw std::runtime_error("gpl_boundary: WriteLine write() failed: " + std::string(::strerror(errno)));
		}
		data += n;
		remaining -= static_cast<size_t>(n);
	}
}

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
