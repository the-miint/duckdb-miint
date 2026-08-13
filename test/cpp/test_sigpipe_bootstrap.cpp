// SIGPIPE backstop for the Catch2 `tests` binary.
//
// The suite exercises gpl-boundary Session / ChildProcess (and Aspera) code that
// writes to subprocess pipes which can close early — e.g. the
// "Session::Initialize throws when child closes stdin without responding" test
// deliberately makes the child exit mid-handshake. In the loaded extension those
// writes are backstopped by SetupSignalHandling()'s process-wide
// `signal(SIGPIPE, SIG_IGN)` (src/miint_extension.cpp:141), installed at load, so
// a broken-pipe write surfaces as EPIPE instead of delivering SIGPIPE. The test
// binary never calls LoadInternal, so without this backstop a broken-pipe write
// hits SIGPIPE's default disposition and kills the whole runner (exit 141),
// aborting the run before Catch2 can report a summary.
//
// Install the same process-wide backstop before any test runs (namespace-scope
// object with static storage duration -> constructed before main()). This
// mirrors production exactly; Session's per-thread pthread_sigmask guards remain
// the finer-grained defense and are unaffected.

#include <csignal>

namespace {

struct SigpipeBackstop {
	SigpipeBackstop() {
		std::signal(SIGPIPE, SIG_IGN);
	}
};

[[maybe_unused]] const SigpipeBackstop g_sigpipe_backstop;

} // namespace
