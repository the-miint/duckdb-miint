#include "gpl_boundary/session.hpp"

#include <array>
#include <cstring>
#include <fcntl.h>
#include <memory>
#include <pthread.h>
#include <signal.h>
#include <sstream>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

#include "yyjson.hpp"

namespace duckdb {
namespace miint {
namespace gpl_boundary {

namespace yj = duckdb_yyjson;

namespace {
// gpl-boundary's protocol version at v0.2.0 (commit 6b11337).
// Bumped 2→3 in that commit when the init reply began advertising the
// registered `tools` array. If the daemon reports a different value, fail
// fast at session boot.
constexpr int kRequiredProtocolVersion = 3;

// Pin the daemon's idle (stdin-silence) auto-shutdown to 60 minutes. Without
// an explicit value the daemon defaults to 60s (GPL-boundary main.rs:
// DEFAULT_IDLE_TIMEOUT_MS), which kills a live session whenever miint takes
// >60s to hand off the next batch — e.g. the per-shard reads⋈read_to_shard
// join blocking on its first chunk, or a slow per-shard index load. That
// surfaces misleadingly as "WriteLine ... Broken pipe [daemon exited with
// code 0]" on the *next* Submit. The daemon is already reaped deterministically
// via PR_SET_PDEATHSIG (parent death) and Session::Shutdown / ~ChildProcess
// SIGTERM, so the idle timer is only a backstop against a wedged-but-alive
// miint; 60 min keeps that backstop while tolerating any realistic batch gap.
constexpr const char *kInitLine = R"({"init":{"idle_timeout_ms":3600000}})";
constexpr const char *kShutdownLine = R"({"shutdown":true})";

// Get a string field from a yyjson object; returns empty string if absent
// or not a string.
std::string get_str(yj::yyjson_val *obj, const char *key) {
	yj::yyjson_val *v = yj::yyjson_obj_get(obj, key);
	if (!v || !yj::yyjson_is_str(v)) {
		return {};
	}
	const char *s = yj::yyjson_get_str(v);
	return s ? std::string(s) : std::string();
}

// RAII guard around a yyjson_doc. yyjson_doc_free is the only resource we
// need to release on every exit path; do it via unique_ptr so refactors can't
// accidentally drop a free.
using YyjsonDocPtr = std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)>;
YyjsonDocPtr make_doc(yj::yyjson_doc *p) {
	return YyjsonDocPtr(p, &yj::yyjson_doc_free);
}

// Cap on how much daemon stderr we retain (and splice into a thrown
// exception). Big enough to carry a typical bowtie2 assertion or stack trace,
// small enough not to pollute SQL error displays. Excess is truncated from the
// front (we want the tail — that's where the failing-batch output lives).
constexpr std::size_t kStderrTailCap = 4096;

// Non-blocking drain of an fd, appending everything currently available to
// `sink`. Used per-batch (Session::PumpStderr — keeps the daemon's stderr pipe
// from filling, which would wedge the daemon on a blocked write while we wait
// for its batch response) and on the failure paths (to grab the daemon's last
// words). MUST be non-blocking — the daemon is typically still alive, so the
// pipe has no EOF and a blocking read would hang.
//
// Restores the original O_NONBLOCK state on exit so other readers (none today,
// but defensive) aren't surprised. Truncation is the caller's job (see
// TruncateHead) so the rolling buffer is only trimmed once after appending.
void AppendAvailableStderr(int fd, std::string &sink) {
	if (fd < 0) {
		return;
	}
	const int orig_flags = ::fcntl(fd, F_GETFL, 0);
	if (orig_flags < 0) {
		return;
	}
	if (!(orig_flags & O_NONBLOCK)) {
		if (::fcntl(fd, F_SETFL, orig_flags | O_NONBLOCK) < 0) {
			return;
		}
	}
	std::array<char, 4096> buf {};
	for (;;) {
		ssize_t n = ::read(fd, buf.data(), buf.size());
		if (n > 0) {
			sink.append(buf.data(), static_cast<size_t>(n));
			continue;
		}
		// n == 0 (EOF) or -1 (EAGAIN/EWOULDBLOCK/other): stop draining.
		break;
	}
	if (!(orig_flags & O_NONBLOCK)) {
		(void)::fcntl(fd, F_SETFL, orig_flags); // best effort restore
	}
}

// Trim `s` from the front to at most `cap` bytes, prefixing a marker so the
// reader knows output was dropped. Keeps the tail — that's where the
// failing-batch output lives.
void TruncateHead(std::string &s, std::size_t cap) {
	if (s.size() > cap) {
		s = "...(truncated head)...\n" + s.substr(s.size() - cap);
	}
}

// Compose a "why did the daemon vanish" suffix for Submit's failure paths.
// Reaps the child to decode its termination — a SIGTERM (e.g. PR_SET_PDEATHSIG
// firing when the DuckDB worker thread that forked the daemon unwinds) reads
// very differently from a SIGKILL (OOM killer) or a SIGSEGV (crash) — and
// splices the daemon's stderr tail. Without this the EPIPE path throws a bare
// "Broken pipe" and the EOF path a bare "closed stdout", neither of which tells
// the operator what actually happened. Reap first so the stderr pipe has hit
// EOF and AppendAvailableStderr captures everything the daemon last wrote;
// `stderr_tail` is the session's rolling buffer (already holding the prior
// batches' output) so the dropped-head marker is applied once after appending.
std::string DaemonDeathSuffix(ChildProcess &child, std::string &stderr_tail) {
	std::string suffix = " [daemon " + child.ReapAndDescribe() + "]";
	AppendAvailableStderr(child.stderr_fd(), stderr_tail);
	TruncateHead(stderr_tail, kStderrTailCap);
	if (!stderr_tail.empty()) {
		suffix += "\n--- daemon stderr ---\n" + stderr_tail;
	}
	return suffix;
}

// RAII guard that blocks SIGPIPE on the calling thread only. `Session::Shutdown`
// uses this to write to a possibly-broken pipe without disturbing other threads
// (DuckDB's scheduler, other sessions) that might rely on the process-wide
// SIGPIPE disposition. `pthread_sigmask` is thread-local; `sigaction` is not.
//
// The destructor drains any pending SIGPIPE on this thread before restoring
// the mask: a blocked-but-pending signal is otherwise delivered the moment
// we unblock, defeating the whole point of the guard.
class ScopedBlockSigpipe {
public:
	ScopedBlockSigpipe() {
		// sigemptyset / sigaddset / sigismember are function-like macros on
		// macOS (`#define sigemptyset(set) (*(set) = 0, 0)`), so we can't
		// qualify them with `::` — the qualifier prevents macro expansion.
		// glibc declares them as both functions and macros; both spellings
		// resolve there. Drop the qualifier for portability.
		sigemptyset(&pipeonly_);
		sigaddset(&pipeonly_, SIGPIPE);
		// Check whether SIGPIPE is already blocked on this thread. If yes,
		// our restore won't change anything; if no, we need to block + drain.
		sigset_t current;
		sigemptyset(&current);
		if (::pthread_sigmask(SIG_BLOCK, &pipeonly_, &current) != 0) {
			active_ = false;
			return;
		}
		active_ = true;
		previously_blocked_ = sigismember(&current, SIGPIPE) == 1;
	}
	~ScopedBlockSigpipe() {
		if (!active_) {
			return;
		}
#if !defined(__APPLE__)
		// Drain any SIGPIPE that landed on this thread while it was blocked.
		// `sigtimedwait` with a zero timeout is the standard way to consume
		// a single pending signal.
		//
		// macOS doesn't implement sigtimedwait (it's a POSIX-RT extension).
		// Skipping the drain is safe there: SetupSignalHandling() registers a
		// process-wide SIG_IGN for SIGPIPE on builds where MIINT_HAS_GPL_BOUNDARY
		// is set, so a pending signal that lands at unblock time is harmlessly
		// ignored rather than delivered to a default-disposition kill handler.
		// The drain is a defense-in-depth optimization on Linux; on macOS it
		// isn't load-bearing.
		struct timespec zero {
			0, 0
		};
		while (::sigtimedwait(&pipeonly_, nullptr, &zero) >= 0) {
			// keep draining
		}
#endif
		if (!previously_blocked_) {
			::pthread_sigmask(SIG_UNBLOCK, &pipeonly_, nullptr);
		}
	}
	ScopedBlockSigpipe(const ScopedBlockSigpipe &) = delete;
	ScopedBlockSigpipe &operator=(const ScopedBlockSigpipe &) = delete;

private:
	sigset_t pipeonly_;
	bool active_ = false;
	bool previously_blocked_ = false;
};
} // namespace

Session::Session(ChildProcess child)
    : child_(std::move(child)), reader_(std::make_unique<LineReader>(child_.stdout_fd())) {
}

Session::~Session() {
	// Best-effort graceful close. Constructor of ChildProcess already arranges
	// SIGTERM-on-destroy as a safety net.
	if (!shut_down_) {
		try {
			Shutdown();
		} catch (...) {
			// Destructor must not throw. The child will be reaped by ~ChildProcess.
		}
	}
}

Session::Session(Session &&other) noexcept
    : child_(std::move(other.child_)), reader_(std::move(other.reader_)), tools_(std::move(other.tools_)),
      initialized_(other.initialized_), shut_down_(other.shut_down_), stderr_tail_(std::move(other.stderr_tail_)) {
	other.initialized_ = false;
	other.shut_down_ = true;
}

Session &Session::operator=(Session &&other) noexcept {
	if (this != &other) {
		this->~Session();
		new (this) Session(std::move(other));
	}
	return *this;
}

void Session::PumpStderr() {
	AppendAvailableStderr(child_.stderr_fd(), stderr_tail_);
	TruncateHead(stderr_tail_, kStderrTailCap);
}

void Session::Initialize() {
	if (initialized_) {
		return;
	}
	WriteLine(child_.stdin_fd(), kInitLine);

	std::string reply;
	if (!reader_->ReadLine(reply)) {
		throw std::runtime_error("gpl_boundary: daemon closed stdout before answering init handshake");
	}

	auto doc = make_doc(yj::yyjson_read(reply.data(), reply.size(), 0));
	if (!doc) {
		throw std::runtime_error("gpl_boundary: init reply was not valid JSON: " + reply);
	}
	yj::yyjson_val *root = yj::yyjson_doc_get_root(doc.get());
	if (!yj::yyjson_is_obj(root)) {
		throw std::runtime_error("gpl_boundary: init reply was not a JSON object: " + reply);
	}
	yj::yyjson_val *success = yj::yyjson_obj_get(root, "success");
	if (!success || !yj::yyjson_is_true(success)) {
		const std::string err = get_str(root, "error");
		throw std::runtime_error("gpl_boundary: init handshake failed: " +
		                         (err.empty() ? std::string("(no error message): ") + reply : err));
	}
	yj::yyjson_val *pv = yj::yyjson_obj_get(root, "protocol_version");
	if (!pv || !yj::yyjson_is_int(pv)) {
		throw std::runtime_error("gpl_boundary: init reply missing or non-integer protocol_version: " + reply);
	}
	const int got = static_cast<int>(yj::yyjson_get_int(pv));
	if (got != kRequiredProtocolVersion) {
		throw std::runtime_error("gpl_boundary: protocol_version mismatch — daemon reports " + std::to_string(got) +
		                         ", miint requires " + std::to_string(kRequiredProtocolVersion) +
		                         ". Update the gpl-boundary binary to a compatible release.");
	}

	// Parse the `tools` array (protocol v3+). Each entry is `{"name": str,
	// "schema_version": int}`. Missing/malformed entries are skipped rather
	// than throwing — the version check above already guarantees we're
	// talking to a v3 daemon, so this is robustness against future
	// additive changes to per-entry fields, not a backward-compat hatch.
	yj::yyjson_val *tools_arr = yj::yyjson_obj_get(root, "tools");
	if (tools_arr && yj::yyjson_is_arr(tools_arr)) {
		const size_t n = yj::yyjson_arr_size(tools_arr);
		tools_.reserve(n);
		for (size_t i = 0; i < n; ++i) {
			yj::yyjson_val *item = yj::yyjson_arr_get(tools_arr, i);
			if (!yj::yyjson_is_obj(item)) {
				continue;
			}
			const std::string name = get_str(item, "name");
			if (name.empty()) {
				continue;
			}
			yj::yyjson_val *sv = yj::yyjson_obj_get(item, "schema_version");
			uint32_t schema_version = 0;
			if (sv && yj::yyjson_is_int(sv)) {
				const int64_t v = yj::yyjson_get_int(sv);
				if (v > 0) {
					schema_version = static_cast<uint32_t>(v);
				}
			}
			tools_.push_back(ToolEntry {name, schema_version});
		}
	}

	initialized_ = true;
}

bool Session::has_tool(const std::string &name) const {
	for (const auto &t : tools_) {
		if (t.name == name) {
			return true;
		}
	}
	return false;
}

uint32_t Session::tool_schema_version(const std::string &name) const {
	for (const auto &t : tools_) {
		if (t.name == name) {
			return t.schema_version;
		}
	}
	return 0;
}

// Minimal RFC 8259 escape for embedding strings inside hand-built JSON. Exposed
// in the public header (session.hpp) so per-tool `config_json` builders share
// one canonical escaper rather than each defining their own. Output does NOT
// include surrounding quotes.
std::string JsonEscape(const std::string &in) {
	std::string out;
	out.reserve(in.size() + 2);
	for (char c : in) {
		switch (c) {
		case '"':
			out += "\\\"";
			break;
		case '\\':
			out += "\\\\";
			break;
		case '\n':
			out += "\\n";
			break;
		case '\r':
			out += "\\r";
			break;
		case '\t':
			out += "\\t";
			break;
		default:
			if (static_cast<unsigned char>(c) < 0x20) {
				char buf[8];
				std::snprintf(buf, sizeof(buf), "\\u%04x", static_cast<unsigned char>(c));
				out += buf;
			} else {
				out += c;
			}
		}
	}
	return out;
}

namespace {
// Build the BatchRequest JSON line. Fields ordered to match gpl-boundary's
// `README.md` "Batch request fields" section. We hand-format rather than
// going through yyjson's writer because the schema is fixed and the cost of
// pulling in a full JSON builder for one writer-side message isn't justified.
std::string build_batch_line(const std::string &tool, const std::string &config_json, const std::string &shm_input_name,
                             std::size_t shm_input_size, int64_t batch_id, bool request_metrics) {
	std::ostringstream os;
	os << R"({"tool":")" << JsonEscape(tool) << R"(",)"
	   << R"("config":)" << config_json << ","
	   << R"("shm_input":")" << JsonEscape(shm_input_name) << R"(",)"
	   << R"("shm_input_size":)" << shm_input_size << ","
	   << R"("batch_id":)" << batch_id;
	// Opt-in only: emit the flag solely when requested, so a non-telemetry batch
	// is byte-for-byte the pre-0.4.2 request shape (older daemons ignore unknown
	// fields, but emitting nothing keeps the contract explicit and free).
	if (request_metrics) {
		os << R"(,"metrics":true)";
	}
	os << "}";
	return os.str();
}
} // namespace

SubmitResult Session::Submit(const std::string &tool, const std::string &config_json, const void *input_bytes,
                             std::size_t input_size, bool request_metrics) {
	if (!initialized_) {
		throw std::runtime_error("gpl_boundary: Session::Submit called before Initialize()");
	}
	if (input_size == 0) {
		throw std::runtime_error("gpl_boundary: Session::Submit requires input_size > 0");
	}

	// Validate `config_json` parses as a single JSON object value before we
	// splice it into the outer batch envelope. Without this, a buggy/malicious
	// caller could pass `{"seed":42},"injected":"value"` and inject extra
	// top-level fields into the request — gpl-boundary would silently honor
	// them. Phase 5's Bind ought to be doing this too; defense-in-depth here.
	{
		auto cfg_doc = make_doc(yj::yyjson_read(config_json.data(), config_json.size(), 0));
		if (!cfg_doc) {
			throw std::runtime_error("gpl_boundary: Submit's config_json is not valid JSON: " + config_json);
		}
		yj::yyjson_val *cfg_root = yj::yyjson_doc_get_root(cfg_doc.get());
		if (!cfg_root || !yj::yyjson_is_obj(cfg_root)) {
			throw std::runtime_error("gpl_boundary: Submit's config_json must be a JSON object: " + config_json);
		}
	}

	// 1. Allocate input shm and copy the caller's IPC bytes in. The region's
	//    destructor unlinks the segment automatically when this function
	//    returns (whether normal exit or exception).
	InputShmRegion in_shm = InputShmRegion::Create(input_size);
	std::memcpy(in_shm.data(), input_bytes, input_size);

	// 2. Reserve a batch_id but do NOT increment until we've successfully sent
	//    the line. If WriteLine throws (EPIPE on a dead daemon), the daemon
	//    never saw this batch_id and we shouldn't burn it for the next call.
	const int64_t batch_id = next_batch_id_;
	const std::string batch_line =
	    build_batch_line(tool, config_json, in_shm.name(), in_shm.size_bytes(), batch_id, request_metrics);
	try {
		WriteLine(child_.stdin_fd(), batch_line);
	} catch (const std::runtime_error &e) {
		// EPIPE here means the daemon's stdin read-end is gone — it died before
		// we finished handing off this batch. Enrich the bare "Broken pipe" with
		// how it died + its stderr so the cause is actionable.
		throw std::runtime_error(std::string(e.what()) + DaemonDeathSuffix(child_, stderr_tail_));
	}
	++next_batch_id_;

	// 3. Read response.
	std::string reply;
	if (!reader_->ReadLine(reply)) {
		// The daemon consumed our batch line but closed stdout without replying:
		// it died between read and write. Same enrichment as the EPIPE path.
		throw std::runtime_error("gpl_boundary: daemon closed stdout while waiting "
		                         "for batch response" +
		                         DaemonDeathSuffix(child_, stderr_tail_));
	}

	// Drain the daemon's stderr now that this batch is answered. Verbose tools
	// (e.g. bowtie2 with quiet=false) print a per-invocation summary to stderr;
	// undrained across many batches that fills the OS pipe buffer and wedges the
	// daemon on a blocked write while we block in the ReadLine above. Draining
	// every batch keeps the pipe near-empty so it can never back up. Cheap
	// no-op when the tool is quiet.
	PumpStderr();

	// 4. Parse response.
	auto doc = make_doc(yj::yyjson_read(reply.data(), reply.size(), 0));
	if (!doc) {
		throw std::runtime_error("gpl_boundary: batch reply was not valid JSON: " + reply);
	}
	yj::yyjson_val *root = yj::yyjson_doc_get_root(doc.get());
	if (!yj::yyjson_is_obj(root)) {
		throw std::runtime_error("gpl_boundary: batch reply was not a JSON object: " + reply);
	}
	yj::yyjson_val *success = yj::yyjson_obj_get(root, "success");
	if (!success || !yj::yyjson_is_true(success)) {
		const std::string err = get_str(root, "error");
		// Fold any final stderr into the rolling tail. The per-batch PumpStderr
		// above already captured this batch's output; this catches bytes the
		// daemon wrote between that drain and the failure. The daemon's own
		// error JSON often says only "subprocess exited unexpectedly"; the
		// actual segfault / assertion / OOM message is in the stderr stream we'd
		// otherwise discard.
		PumpStderr();
		std::string msg = "gpl_boundary: batch failed (batch_id=" + std::to_string(batch_id) + "): ";
		msg += (err.empty() ? std::string("(no error message): ") + reply : err);
		if (!stderr_tail_.empty()) {
			msg += "\n--- daemon stderr ---\n";
			msg += stderr_tail_;
		}
		throw std::runtime_error(msg);
	}

	// MVP correlation: at most one Submit is outstanding at a time (caller
	// blocks for the response before issuing the next), so we simply assert
	// the daemon echoed our batch_id rather than maintaining a map of
	// in-flight requests. If/when miint ever pipelines multiple Submits on
	// the same connection (gpl-boundary supports parallel batches across
	// distinct fingerprints), the assert below has to be replaced with a
	// pid_t→Promise dispatch. Until then: assert.
	SubmitResult result;
	yj::yyjson_val *bid = yj::yyjson_obj_get(root, "batch_id");
	if (bid && yj::yyjson_is_int(bid)) {
		result.batch_id = yj::yyjson_get_int(bid);
		if (result.batch_id != batch_id) {
			throw std::runtime_error("gpl_boundary: batch_id mismatch — sent " + std::to_string(batch_id) +
			                         " but daemon replied with " + std::to_string(result.batch_id));
		}
	} else {
		// daemon should always echo batch_id back; missing it is a protocol violation
		throw std::runtime_error("gpl_boundary: batch reply missing batch_id field: " + reply);
	}

	yj::yyjson_val *sv = yj::yyjson_obj_get(root, "schema_version");
	if (sv && yj::yyjson_is_int(sv)) {
		result.schema_version = static_cast<uint32_t>(yj::yyjson_get_int(sv));
	}

	// shm_outputs is the meat. Open each segment with the explicit size.
	// Manual loop instead of yyjson_arr_foreach macro — the macro expands
	// without namespace qualification and breaks our `yj::` alias.
	yj::yyjson_val *outs = yj::yyjson_obj_get(root, "shm_outputs");
	if (outs && yj::yyjson_is_arr(outs)) {
		const size_t arr_size = yj::yyjson_arr_size(outs);
		for (size_t i = 0; i < arr_size; ++i) {
			yj::yyjson_val *item = yj::yyjson_arr_get(outs, i);
			if (!yj::yyjson_is_obj(item)) {
				throw std::runtime_error("gpl_boundary: shm_outputs entry is not an object");
			}
			const std::string out_name = get_str(item, "name");
			const std::string out_label = get_str(item, "label");
			yj::yyjson_val *size_v = yj::yyjson_obj_get(item, "size");
			if (out_name.empty() || !size_v || !yj::yyjson_is_int(size_v)) {
				throw std::runtime_error("gpl_boundary: shm_outputs entry missing name/size");
			}
			const int64_t size_signed = yj::yyjson_get_int(size_v);
			if (size_signed <= 0) {
				throw std::runtime_error("gpl_boundary: shm_outputs entry has non-positive size: " +
				                         std::to_string(size_signed));
			}
			SubmitOutput out {out_label, OutputShmRegion::Open(out_name, static_cast<std::size_t>(size_signed))};
			result.outputs.push_back(std::move(out));
		}
	}

	// Optional `result` field — capture the raw JSON for callers that want it
	// (e.g., `phylogeny_fasttree` will extract `n_nodes`/`n_leaves` from here).
	yj::yyjson_val *res = yj::yyjson_obj_get(root, "result");
	if (res) {
		size_t len = 0;
		char *json_str = yj::yyjson_val_write(res, 0, &len);
		if (json_str) {
			result.result_json.assign(json_str, len);
			std::free(json_str);
		}
	}

	// Optional `metrics` field — additive, opt-in batch metrics (e.g. the
	// bowtie2 worker's getrusage). Captured raw exactly like `result`; old
	// daemons omit it and callers see an empty string. Folded into
	// align_bowtie2_sharded's per-batch telemetry when present.
	yj::yyjson_val *met = yj::yyjson_obj_get(root, "metrics");
	if (met) {
		size_t len = 0;
		char *json_str = yj::yyjson_val_write(met, 0, &len);
		if (json_str) {
			result.metrics_json.assign(json_str, len);
			std::free(json_str);
		}
	}

	return result;
}

void Session::Shutdown() {
	if (shut_down_) {
		return;
	}
	shut_down_ = true;

	// Best-effort send of the shutdown line. The daemon may already be gone
	// (e.g., crashed mid-batch). EPIPE on write is benign here — the child
	// destructor will SIGTERM anything left running.
	//
	// Block SIGPIPE on this thread only so a closed pipe surfaces as
	// -1/EPIPE rather than killing the whole process. Thread-scoped (vs
	// process-wide `sigaction`) is critical: DuckDB scheduler threads and
	// other miint subsystems may be writing to other pipes concurrently and
	// must keep the host's SIGPIPE disposition.
	if (child_.stdin_fd() >= 0) {
		ScopedBlockSigpipe block_sigpipe;
		try {
			WriteLine(child_.stdin_fd(), kShutdownLine);
		} catch (const std::runtime_error &) {
			// Pipe closed / write failed. Daemon is likely already exiting.
		}
	}
	child_.CloseStdin();
	child_.Wait();
}

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
