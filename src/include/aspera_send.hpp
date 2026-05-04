// SPDX-License-Identifier: MIT
//
// Aspera write-side: fork+execvp `ascp --mode=send ...`. Used by the
// `ena_upload_reads` table function. Argv is built via
// `duckdb::BuildAscpSendArgv` (DuckDB-free, unit-tested in
// test/cpp/test_ena_upload_reads.cpp); this file owns the process side —
// fork, env injection of ASPERA_SCP_PASS so the password never lands in
// /proc/self/cmdline, stderr capture, and waitpid.

#pragma once

#include "aspera_utils.hpp" // for MIINT_ASPERA_SUPPORTED

#include <string>
#include <vector>

namespace miint {

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)

struct AsperaSendResult {
	int exit_code; // 0 = success
	std::string stderr_output;
};

// Run `ascp` with the given argv and password. The password is injected via
// the ASPERA_SCP_PASS environment variable in the child process so it never
// appears in argv (visible to anyone via /proc/<pid>/cmdline).
//
// Returns the child's exit code and the captured stderr. The caller decides
// whether to throw — typically a non-zero exit code → IOException with the
// stderr output. Errors that prevent the child from starting (fork/pipe
// failure) raise duckdb::IOException directly from this function.
AsperaSendResult RunAsperaSend(const std::vector<std::string> &argv, const std::string &password);

#endif // MIINT_ASPERA_SUPPORTED && MIINT_STATIC_BUILD

} // namespace miint
