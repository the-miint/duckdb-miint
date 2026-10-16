#include "gpl_boundary/process.hpp"

#include <array>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <mutex>
#include <poll.h>
#include <signal.h>
#include <stdexcept>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

#if defined(__linux__)
#include <sys/prctl.h>
#endif
#if defined(__APPLE__)
#include <crt_externs.h>
#else
extern char **environ; // NOLINT(cppcoreguidelines-avoid-non-const-global-variables): POSIX's own
#endif

namespace duckdb {
namespace miint {
namespace gpl_boundary {

namespace {
#if defined(__linux__)
// pipe2 creates both ends close-on-exec atomically, so no lock is needed.
struct SpawnLock {};
#else
// Without pipe2, a pipe is inheritable between pipe() and fcntl(). Holding
// this across pipe creation and fork keeps our own spawns out of that window;
// a fork elsewhere in the process is not covered.
std::mutex &SpawnMutex() {
	static std::mutex mutex;
	return mutex;
}
struct SpawnLock {
	std::lock_guard<std::mutex> lock {SpawnMutex()};
};
#endif

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

// Make a pipe whose ends are both close-on-exec, so no child inherits them:
// not this one's (it gets its ends on fds 0-2 via WireStdFd) and not any
// other child, spawned before or after. An inherited end keeps a pipe open,
// so a child would never see EOF on its stdin while another child lives.
void make_pipe(int fds[2]) {
#if defined(__linux__)
	if (::pipe2(fds, O_CLOEXEC) == -1) {
		throw std::runtime_error("gpl_boundary: pipe2() failed: " + std::string(::strerror(errno)));
	}
#else
	if (::pipe(fds) == -1) {
		throw std::runtime_error("gpl_boundary: pipe() failed: " + std::string(::strerror(errno)));
	}
	if (::fcntl(fds[0], F_SETFD, FD_CLOEXEC) == -1 || ::fcntl(fds[1], F_SETFD, FD_CLOEXEC) == -1) {
		const int err = errno;
		::close(fds[0]);
		::close(fds[1]);
		throw std::runtime_error("gpl_boundary: fcntl(FD_CLOEXEC) failed: " + std::string(::strerror(err)));
	}
#endif
}

// In a forked child (async-signal-safe): put `fd` on `target` (0, 1 or 2) so
// it survives exec. Copying above fd 2 first matters when the parent runs
// with a standard fd closed: pipes then get numbers 0-2, and a direct dup2
// could overwrite one pipe end with another, or be a no-op that leaves
// FD_CLOEXEC set. The copy closes at exec; dup2 clears the flag on `target`.
bool WireStdFd(int fd, int target) {
	const int copy = ::fcntl(fd, F_DUPFD_CLOEXEC, 3);
	return copy >= 0 && ::dup2(copy, target) >= 0;
}

char **CurrentEnviron() {
#if defined(__APPLE__)
	return *::_NSGetEnviron(); // a dylib can't reference `environ` directly
#else
	return environ;
#endif
}

// The parent's environment with `options` applied, as "NAME=value" entries.
// Built before fork: allocating in the child is not async-signal-safe.
std::vector<std::string> BuildEnvironment(const ChildProcessOptions &options) {
	auto replaced = [&options](const std::string &name) {
		for (const auto &unset : options.unset_env) {
			if (unset == name) {
				return true;
			}
		}
		for (const auto &kv : options.set_env) {
			if (kv.first == name) {
				return true;
			}
		}
		return false;
	};
	std::vector<std::string> env;
	for (char **entry = CurrentEnviron(); entry && *entry; ++entry) {
		const std::string text(*entry);
		if (!replaced(text.substr(0, text.find('=')))) {
			env.push_back(text);
		}
	}
	for (const auto &kv : options.set_env) {
		env.push_back(kv.first + "=" + kv.second);
	}
	return env;
}

// Decode a raw waitpid status into a human-readable, signal-named string.
// Kept next to make_pipe so the process-control vocabulary lives in one place.
std::string DescribeWaitStatus(int status) {
	if (WIFSIGNALED(status)) {
		const int sig = WTERMSIG(status);
		const char *name = ::strsignal(sig);
		return "killed by signal " + std::to_string(sig) + " (" + (name ? name : "unknown") + ")";
	}
	if (WIFEXITED(status)) {
		return "exited with code " + std::to_string(WEXITSTATUS(status));
	}
	return "terminated with raw wait status " + std::to_string(status);
}
} // namespace

std::string FindExecutableInPath(const std::string &name) {
	int pipefd[2];
	pid_t pid;
	{
		SpawnLock spawn_lock;
		try {
			make_pipe(pipefd);
		} catch (const std::runtime_error &) {
			return "";
		}
		pid = ::fork();
	}
	if (pid == -1) {
		::close(pipefd[0]);
		::close(pipefd[1]);
		return "";
	}

	if (pid == 0) {
		// Child. Every pipe end is close-on-exec; only the copies wired onto
		// stdout and stderr survive exec.
		if (!WireStdFd(pipefd[1], STDOUT_FILENO)) {
			::_exit(127);
		}
		int devnull = ::open("/dev/null", O_WRONLY | O_CLOEXEC);
		if (devnull != -1) {
			(void)WireStdFd(devnull, STDERR_FILENO);
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

ChildProcess::ChildProcess(const std::vector<std::string> &argv, const ChildProcessOptions &options) {
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

	SpawnLock spawn_lock;
	try {
		make_pipe(in_pipe);
		make_pipe(out_pipe);
		make_pipe(err_pipe);
		// Close-on-exec like every pipe here: the parent reads "" (EOF) when
		// execv succeeds (the kernel closes both fds across the exec), or the
		// raw errno value when the child wrote one before _exit. This avoids
		// the race where waitpid(WNOHANG) returns 0 because the child hasn't
		// yet had time to fail exec.
		make_pipe(exec_status_pipe);
	} catch (...) {
		close_all();
		throw;
	}

	// Build argv (and, when options change it, the environment) as char*
	// arrays with terminating nullptrs, before fork.
	std::vector<std::string> owned = argv;
	std::vector<char *> c_argv;
	c_argv.reserve(owned.size() + 1);
	for (auto &s : owned) {
		c_argv.push_back(s.data());
	}
	c_argv.push_back(nullptr);
	const bool custom_env = !options.set_env.empty() || !options.unset_env.empty();
	std::vector<std::string> env_owned;
	std::vector<char *> c_envp;
	if (custom_env) {
		env_owned = BuildEnvironment(options);
		c_envp.reserve(env_owned.size() + 1);
		for (auto &s : env_owned) {
			c_envp.push_back(s.data());
		}
		c_envp.push_back(nullptr);
	}

	pid_t pid = ::fork();
	if (pid == -1) {
		const int err = errno;
		close_all();
		throw std::runtime_error("gpl_boundary: fork() failed: " + std::string(::strerror(err)));
	}

	if (pid == 0) {
		// Child. Async-signal-safe path only — execv soon clears the slate.
		// Every pipe end is close-on-exec, so the only ones that survive exec
		// are the copies wired onto 0-2, and the parent sees EOF on its read
		// ends once the child exits.
		bool ok = WireStdFd(in_pipe[0], STDIN_FILENO) && WireStdFd(out_pipe[1], STDOUT_FILENO) &&
		          WireStdFd(err_pipe[1], STDERR_FILENO);
		if (ok && options.new_session) {
			ok = ::setsid() != -1;
		}
#if defined(__linux__)
		if (ok && options.parent_death_signal) {
			::prctl(PR_SET_PDEATHSIG, SIGTERM);
		}
#endif
		if (ok) {
			if (custom_env) {
				::execve(c_argv[0], c_argv.data(), c_envp.data());
			} else {
				::execv(c_argv[0], c_argv.data());
			}
		}

		// Setup or exec failed — write errno into the status pipe so the
		// parent knows why. write/_exit are async-signal-safe.
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
	new_session_ = options.new_session;
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
		// Setup (fd wiring, setsid) or exec itself failed.
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
		throw std::runtime_error("gpl_boundary: failed to start '" + program +
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
	//
	// A session leader's whole group gets SIGKILL instead: SIGTERM to the
	// leader alone would leave whatever it started running. Such a host's
	// graceful exit is stdin EOF, which its owner sends before destruction.
	if (pid_ > 0 && !waited_) {
		if (new_session_) {
			::killpg(pid_, SIGKILL);
		} else {
			::kill(pid_, SIGTERM);
		}
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
      waited_(other.waited_), wait_status_(other.wait_status_), new_session_(other.new_session_) {
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

std::string ChildProcess::ReapAndDescribe(int grace_ms) {
	if (waited_) {
		return DescribeWaitStatus(wait_status_);
	}
	if (pid_ <= 0) {
		return "no child process";
	}
	// Poll WNOHANG up to grace_ms. On the daemon-death paths the child is a
	// zombie already (the pipe EOF/EPIPE that brought us here implies it exited),
	// so the first poll usually succeeds; the grace window only covers the rare
	// race where the kernel hasn't transitioned it yet.
	for (int elapsed = 0;; elapsed += 10) {
		int status = 0;
		const pid_t r = ::waitpid(pid_, &status, WNOHANG);
		if (r == pid_) {
			wait_status_ = status;
			waited_ = true;
			return DescribeWaitStatus(status);
		}
		if (r == -1) {
			// ECHILD: someone else already reaped it (e.g. a process-wide
			// SIGCHLD handler) — itself a useful tell. Mark waited_ so the
			// destructor and Wait() don't block on a pid that's gone.
			const std::string err = ::strerror(errno);
			waited_ = true;
			return "already reaped or unwaitable (" + err + ")";
		}
		// r == 0: still running.
		if (elapsed >= grace_ms) {
			return "still running after " + std::to_string(grace_ms) + "ms (no exit yet)";
		}
		std::this_thread::sleep_for(std::chrono::milliseconds(10));
	}
}

void ChildProcess::KillGroup() {
	if (!new_session_) {
		throw std::logic_error("gpl_boundary: KillGroup requires ChildProcessOptions::new_session");
	}
	if (pid_ <= 0 || waited_) {
		return;
	}
	// ESRCH (or EPERM on macOS, when only zombies remain) means already gone.
	(void)::killpg(pid_, SIGKILL);
}

namespace {
// Read both fds until both reach EOF, appending to the matching sink. Reading
// one to EOF before the other deadlocks once the child fills the other pipe.
void DrainBoth(int out_fd, int err_fd, std::string &out, std::string &err) {
	pollfd fds[2] = {{out_fd, POLLIN, 0}, {err_fd, POLLIN, 0}};
	std::string *sinks[2] = {&out, &err};
	std::array<char, 4096> buf {};
	int open = 2;
	while (open > 0) {
		if (::poll(fds, 2, -1) < 0) {
			if (errno == EINTR) {
				continue;
			}
			break; // any other error: bail with whatever we got
		}
		for (int i = 0; i < 2; ++i) {
			if (fds[i].fd < 0 || fds[i].revents == 0) {
				continue;
			}
			const ssize_t n = ::read(fds[i].fd, buf.data(), buf.size());
			if (n > 0) {
				sinks[i]->append(buf.data(), static_cast<size_t>(n));
			} else if (n == 0 || errno != EINTR) {
				fds[i].fd = -1; // EOF or error: poll ignores a negative fd
				--open;
			}
		}
	}
}
} // namespace

RunResult RunCapture(const std::vector<std::string> &argv, const ChildProcessOptions &options) {
	RunResult r;
	try {
		ChildProcess child(argv, options);
		// Close stdin so the child immediately sees EOF and won't block
		// waiting for input on a one-shot probe.
		child.CloseStdin();
		DrainBoth(child.stdout_fd(), child.stderr_fd(), r.stdout_text, r.stderr_text);
		const int status = child.Wait();
		r.spawned = true;
		if (WIFEXITED(status)) {
			r.exited_normally = true;
			r.exit_code = WEXITSTATUS(status);
		} else if (WIFSIGNALED(status)) {
			r.term_signal = WTERMSIG(status);
		}
	} catch (const std::exception &e) {
		r.spawn_error = e.what();
	}
	return r;
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
		FillOnce();
	}
}

bool LineReader::HasLine() const {
	return eof_ || HasCompleteLine();
}

bool LineReader::HasCompleteLine() const {
	return buffer_.find('\n') != std::string::npos;
}

void LineReader::FillOnce() {
	char chunk[4096];
	ssize_t n;
	do {
		n = ::read(fd_, chunk, sizeof(chunk));
	} while (n == -1 && errno == EINTR);
	if (n < 0) {
		throw std::runtime_error("gpl_boundary: LineReader read() failed: " + std::string(::strerror(errno)));
	}
	if (n == 0) {
		eof_ = true;
		return;
	}
	buffer_.append(chunk, static_cast<size_t>(n));
}

// =============================================================================
// stderr capture
// =============================================================================

void TruncateHead(std::string &s, std::size_t cap) {
	if (s.size() > cap) {
		s = "...(truncated head)...\n" + s.substr(s.size() - cap);
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
