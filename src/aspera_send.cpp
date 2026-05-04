// SPDX-License-Identifier: MIT
//
// See aspera_send.hpp.

#include "aspera_send.hpp"

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)

#include "duckdb/common/exception.hpp"

#include <array>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __linux__
#include <sys/prctl.h>
#endif

extern char **environ; // NOLINT(readability-redundant-declaration)

namespace miint {

AsperaSendResult RunAsperaSend(const std::vector<std::string> &argv, const std::string &password) {
	if (argv.empty()) {
		throw duckdb::IOException("RunAsperaSend: argv is empty");
	}

	// Build the C-style argv. Pointers reference storage owned by `argv` so
	// they stay valid until execvp.
	std::vector<char *> c_argv;
	c_argv.reserve(argv.size() + 1);
	for (const auto &a : argv) {
		c_argv.push_back(const_cast<char *>(a.c_str()));
	}
	c_argv.push_back(nullptr);

	int stderr_pipe[2];
	if (pipe(stderr_pipe) == -1) {
		throw duckdb::IOException("RunAsperaSend: failed to create stderr pipe: %s", std::strerror(errno));
	}

	pid_t pid = fork();
	if (pid == -1) {
		close(stderr_pipe[0]);
		close(stderr_pipe[1]);
		throw duckdb::IOException("RunAsperaSend: fork failed: %s", std::strerror(errno));
	}

	if (pid == 0) {
		// Child.
#ifdef __linux__
		prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif
		// Redirect stderr to the parent-readable pipe; route stdout to
		// /dev/null (ascp progress lines are noise we don't capture); close
		// stdin. If /dev/null can't be opened, refuse to inherit the
		// parent's stdout into the child rather than leak it.
		dup2(stderr_pipe[1], STDERR_FILENO);
		close(stderr_pipe[0]);
		close(stderr_pipe[1]);
		int devnull = open("/dev/null", O_RDWR);
		if (devnull == -1) {
			_exit(127);
		}
		dup2(devnull, STDOUT_FILENO);
		dup2(devnull, STDIN_FILENO);
		close(devnull);

		// Inject the password via ASPERA_SCP_PASS so it never appears in
		// argv (visible via /proc/<pid>/cmdline). `setenv` copies the
		// string into the child's environ block, so it's fine for the
		// caller's `password` to be destroyed after fork.
		setenv("ASPERA_SCP_PASS", password.c_str(), /*overwrite=*/1);

		execvp(c_argv[0], c_argv.data());
		// execvp only returns on failure.
		_exit(127);
	}

	// Parent: drain stderr until child closes its end, then waitpid.
	close(stderr_pipe[1]);
	std::string stderr_buf;
	std::array<char, 4096> buf;
	while (true) {
		ssize_t n = read(stderr_pipe[0], buf.data(), buf.size());
		if (n > 0) {
			stderr_buf.append(buf.data(), static_cast<std::size_t>(n));
		} else if (n == 0) {
			break;
		} else if (errno == EINTR) {
			continue; // retry, signal interrupted the read
		} else {
			break; // genuine error; we'll surface it via waitpid status
		}
	}
	close(stderr_pipe[0]);

	int status = 0;
	while (true) {
		const pid_t r = waitpid(pid, &status, 0);
		if (r == pid) {
			break;
		}
		if (r == -1 && errno == EINTR) {
			continue;
		}
		throw duckdb::IOException("RunAsperaSend: waitpid failed: %s", std::strerror(errno));
	}
	int exit_code = -1;
	if (WIFEXITED(status)) {
		exit_code = WEXITSTATUS(status);
	} else if (WIFSIGNALED(status)) {
		// Surface signal-killed children as a distinct non-zero code so the
		// caller's diagnostic mentions it without us needing a separate enum.
		exit_code = 128 + WTERMSIG(status);
	}
	return AsperaSendResult {exit_code, std::move(stderr_buf)};
}

} // namespace miint

#endif // MIINT_ASPERA_SUPPORTED && MIINT_STATIC_BUILD
