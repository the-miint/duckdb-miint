// SPDX-License-Identifier: MIT
//
// Tests that a child process, and everything it started, has exited.

#pragma once

#include <catch2/catch_test_macros.hpp>

#include <cerrno>
#include <chrono>
#include <poll.h>
#include <unistd.h>

namespace miint_test {

// True once `fd` reaches EOF within `ms`; false on timeout. Discards data.
inline bool EofWithin(int fd, int ms) {
	const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(ms);
	char buf[256];
	for (;;) {
		const auto left =
		    std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::steady_clock::now()).count();
		if (left <= 0) {
			return false;
		}
		pollfd p {fd, POLLIN, 0};
		const int ready = ::poll(&p, 1, static_cast<int>(left));
		if (ready < 0 && errno == EINTR) {
			continue;
		}
		if (ready <= 0) {
			return false;
		}
		const ssize_t n = ::read(fd, buf, sizeof(buf));
		if (n == 0) {
			return true;
		}
		if (n < 0 && errno != EINTR) {
			return false;
		}
	}
}

// A pipe whose write end every child (and grandchild) inherits: plain pipe(),
// not close-on-exec. Once the test closes its own write end, EOF on the read
// end means every process holding it has exited. Detecting death this way,
// not with kill(pid, 0), holds under a non-reaping PID 1 and on macOS.
struct Witness {
	int read_fd = -1;
	int write_fd = -1;
	Witness() {
		int fds[2];
		REQUIRE(::pipe(fds) == 0);
		read_fd = fds[0];
		write_fd = fds[1];
	}
	void CloseWriteEnd() {
		::close(write_fd);
		write_fd = -1;
	}
	~Witness() {
		::close(read_fd);
		if (write_fd >= 0) {
			::close(write_fd);
		}
	}
	Witness(const Witness &) = delete;
	Witness &operator=(const Witness &) = delete;
};

} // namespace miint_test
