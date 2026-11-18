// An OpenMP critical section keeps working after the thread that first used OpenMP has exited
// (the-miint/duckdb-miint#482).
//
// Built only where libomp is linked statically (cmake/openmp.cmake, MIINT_OPENMP_STATIC_LIBOMP).
// That runtime shuts down when the last thread that used it exits and initializes again on the next
// call; the shutdown frees the lock an `omp critical` left in a variable in the caller's binary, the
// variable keeps the stale pointer, and the next critical section segfaults. DuckDB runs table
// functions on threads that exit when their database closes, so krepp's criticals hit exactly this.
// Loading miint has to keep the runtime up (src/openmp_runtime_keepalive.cpp).
//
// The critical section here is the test's own, so the test does not depend on which libraries in
// miint happen to use one.
//
// It needs a process of its own. Once any thread that stays alive has used OpenMP, the runtime never
// shuts down and the sequence below cannot fail, and other tests in this binary may already have
// done that on the main thread. So the test re-runs this binary on a hidden case and reads the
// child's exit status; Catch2's macros report from the child as usual.

#include <catch2/catch_test_macros.hpp>

#include "duckdb/main/database.hpp"

#include <mach-o/dyld.h>
#include <sys/wait.h>
#include <unistd.h>

#include <atomic>
#include <cstdint>
#include <string>
#include <thread>
#include <vector>

namespace {

std::atomic<int> critical_entries {0};

//! Enters an unnamed critical section on a new thread, which then exits.
void EnterCriticalOnThreadThatExits() {
	std::thread worker([] {
#pragma omp critical
		{ critical_entries++; }
	});
	worker.join();
}

std::string ExecutablePath() {
	uint32_t size = 0;
	_NSGetExecutablePath(nullptr, &size);
	std::vector<char> buffer(size);
	REQUIRE(_NSGetExecutablePath(buffer.data(), &size) == 0);
	return std::string(buffer.data());
}

} // namespace

TEST_CASE("an OpenMP critical section works again after the first OpenMP thread exits", "[openmp]") {
	const std::string self = ExecutablePath();
	std::string test_spec = "[openmp_runtime_child]";
	std::vector<char *> argv {const_cast<char *>(self.c_str()), test_spec.data(), nullptr};
	// fork + execv, as src/gpl_boundary/process.cpp does: the exec gives the child a fresh runtime.
	const pid_t child = fork();
	REQUIRE(child >= 0);
	if (child == 0) {
		execv(self.c_str(), argv.data());
		_exit(127);
	}
	int status = 0;
	REQUIRE(waitpid(child, &status, 0) == child);
	INFO("child " << (WIFSIGNALED(status) ? "killed by signal " + std::to_string(WTERMSIG(status))
	                                      : "exited with status " + std::to_string(WEXITSTATUS(status))));
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
}

// Hidden: run only by the case above, in its own process.
TEST_CASE("child: two critical sections, each on a thread that exits", "[.openmp_runtime_child]") {
	// Loading miint (statically linked into this binary) is what has to keep the runtime alive.
	duckdb::DuckDB db(nullptr);
	EnterCriticalOnThreadThatExits();
	EnterCriticalOnThreadThatExits();
	REQUIRE(critical_entries == 2);
}
