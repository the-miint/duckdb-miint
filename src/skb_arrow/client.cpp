#include "skb_arrow/client.hpp"

#include <cerrno>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <dirent.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <type_traits>
#include <unistd.h>

#include "nanoarrow/nanoarrow.h"
#include "nanoarrow/nanoarrow_ipc.h"
#include "yyjson.hpp"

namespace duckdb {
namespace miint {
namespace skb_arrow {

namespace gb = ::duckdb::miint::gpl_boundary;
namespace yj = duckdb_yyjson;

namespace {

// How often a call waiting on the host checks for an interrupt.
constexpr int kPollMs = 100;
// How long to collect a dead host's stderr. skb-arrow's stderr drainer holds
// an unterminated last line for up to 50 ms, and reaches EOF only once the
// host and the drainer have both exited.
constexpr int kStderrEofMs = 2000;
// How much of a dead host's stderr goes into the error: the tail.
constexpr std::size_t kStderrTailCap = 4096;

bool IsResourceErrno(int err) {
	return err == ENOSPC || err == EDQUOT || err == ENOMEM;
}

std::string Errno(int err) {
	return ::strerror(err);
}

// docs/transport.md#segments: [A-Za-z0-9_-][A-Za-z0-9._-]{0,127}.
bool IsSegmentName(const std::string &name) {
	if (name.empty() || name.size() > 128) {
		return false;
	}
	for (std::size_t i = 0; i < name.size(); ++i) {
		const char c = name[i];
		const bool alnum = (c >= 'A' && c <= 'Z') || (c >= 'a' && c <= 'z') || (c >= '0' && c <= '9');
		if (!alnum && c != '_' && c != '-' && !(c == '.' && i > 0)) {
			return false;
		}
	}
	return true;
}

// Read `fd` until EOF or `ms` pass, whichever comes first.
std::string ReadToEof(int fd, int ms) {
	std::string out;
	const auto deadline = std::chrono::steady_clock::now() + std::chrono::milliseconds(ms);
	char buf[4096];
	for (;;) {
		const auto left =
		    std::chrono::duration_cast<std::chrono::milliseconds>(deadline - std::chrono::steady_clock::now()).count();
		if (left <= 0) {
			return out;
		}
		pollfd p {fd, POLLIN, 0};
		const int ready = ::poll(&p, 1, static_cast<int>(left));
		if (ready < 0 && errno == EINTR) {
			continue;
		}
		if (ready <= 0) {
			return out;
		}
		const ssize_t n = ::read(fd, buf, sizeof(buf));
		if (n > 0) {
			out.append(buf, static_cast<std::size_t>(n));
		} else if (n == 0 || errno != EINTR) {
			return out;
		}
	}
}

// --- JSON ------------------------------------------------------------------

using DocPtr = std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)>;
using MutDocPtr = std::unique_ptr<yj::yyjson_mut_doc, decltype(&yj::yyjson_mut_doc_free)>;

// A response line parsed as a JSON object, or a host_incompatible Error.
DocPtr ParseResponse(const std::string &line) {
	DocPtr doc(yj::yyjson_read(line.data(), line.size(), 0), &yj::yyjson_doc_free);
	if (!doc || !yj::yyjson_is_obj(yj::yyjson_doc_get_root(doc.get()))) {
		std::string shown = line.substr(0, 200);
		if (line.size() > 200) {
			shown += "...";
		}
		throw Error("host_incompatible", "skb-arrow sent a response that is not a JSON object: " + shown);
	}
	return doc;
}

// The string at `key`; `present` says whether it was a string at all.
std::string GetString(yj::yyjson_val *obj, const char *key, bool *present = nullptr) {
	yj::yyjson_val *v = yj::yyjson_obj_get(obj, key);
	const bool is_str = v && yj::yyjson_is_str(v);
	if (present) {
		*present = is_str;
	}
	return is_str ? std::string(yj::yyjson_get_str(v), yj::yyjson_get_len(v)) : std::string();
}

[[noreturn]] void ThrowMalformed(const std::string &what) {
	throw Error("host_incompatible", "skb-arrow sent a malformed response: " + what);
}

// docs/errors.md: a kind this client doesn't know is `internal`, and only an
// `internal` error carries a traceback.
Error ErrorFromResponse(yj::yyjson_val *root) {
	static const char *const kKinds[] = {"host_incompatible", "invalid_input", "invalid_param",
	                                     "unsupported",       "resource",      "internal"};
	std::string kind = GetString(root, "kind");
	std::string message = GetString(root, "message");
	bool known = false;
	for (const char *k : kKinds) {
		known = known || kind == k;
	}
	if (!known) {
		message = "skb-arrow reported an error of unknown kind '" + kind + "': " + message;
		return Error("internal", message);
	}
	return Error(kind, message, kind == "internal" ? GetString(root, "traceback") : std::string());
}

std::vector<Warning> ParseWarnings(yj::yyjson_val *root) {
	std::vector<Warning> warnings;
	yj::yyjson_val *arr = yj::yyjson_obj_get(root, "warnings");
	if (!arr) {
		return warnings;
	}
	if (!yj::yyjson_is_arr(arr)) {
		ThrowMalformed("'warnings' is not an array");
	}
	for (std::size_t i = 0; i < yj::yyjson_arr_size(arr); ++i) {
		yj::yyjson_val *w = yj::yyjson_arr_get(arr, i);
		yj::yyjson_val *count = yj::yyjson_obj_get(w, "count");
		bool has_category = false;
		bool has_message = false;
		Warning warning;
		warning.category = GetString(w, "category", &has_category);
		warning.message = GetString(w, "message", &has_message);
		if (!has_category || !has_message || !count || !yj::yyjson_is_int(count)) {
			ThrowMalformed("a warning lacks a string 'category', a string 'message' or an integer 'count'");
		}
		warning.count = yj::yyjson_get_sint(count);
		warnings.push_back(std::move(warning));
	}
	return warnings;
}

// A `result` response: its output segment names and warnings.
CallResult ParseResult(yj::yyjson_val *root) {
	CallResult result;
	yj::yyjson_val *output = yj::yyjson_obj_get(root, "output");
	if (!output || !yj::yyjson_is_arr(output) || yj::yyjson_arr_size(output) == 0) {
		ThrowMalformed("'output' is not a non-empty array");
	}
	for (std::size_t i = 0; i < yj::yyjson_arr_size(output); ++i) {
		yj::yyjson_val *name = yj::yyjson_arr_get(output, i);
		const std::string segment =
		    yj::yyjson_is_str(name) ? std::string(yj::yyjson_get_str(name), yj::yyjson_get_len(name)) : "";
		if (!IsSegmentName(segment)) {
			ThrowMalformed("'output' names an invalid segment");
		}
		result.output.push_back(segment);
	}
	result.warnings = ParseWarnings(root);
	return result;
}

// --- spawning --------------------------------------------------------------

gb::ChildProcess SpawnHost(const HostOptions &options, const std::string &dir) {
	gb::ChildProcessOptions spawn;
	// The host cleans up after itself when its caller goes away (stdin EOF,
	// or stdout losing its reader), and PDEATHSIG would fire when the
	// spawning DuckDB *thread* exits, mid-query.
	spawn.parent_death_signal = false;
	// Its own process group, so cancelling stops everything it started.
	spawn.new_session = true;
	// A Python-hosted DuckDB's interpreter settings must not reach the host.
	spawn.unset_env = {"PYTHONHOME", "PYTHONPATH"};
	spawn.set_env = {{"PYTHONNOUSERSITE", "1"}};
	if (options.threads > 0) {
		const std::string n = std::to_string(options.threads);
		for (const char *var :
		     {"NUMBA_NUM_THREADS", "OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"}) {
			spawn.set_env.emplace_back(var, n);
		}
	}
	try {
		return gb::ChildProcess({options.binary, "--segment-dir", dir}, spawn);
	} catch (const std::runtime_error &e) {
		// Most often a launcher whose interpreter is gone: the installation's
		// problem, not miint's.
		throw Error("host_incompatible", "cannot start skb-arrow at '" + options.binary + "': " + e.what() +
		                                     ". If its installation is broken, reinstall it: "
		                                     "SELECT install_skb_arrow(true);");
	}
}

} // namespace

// =============================================================================
// Error, discovery, placement
// =============================================================================

Error::Error(std::string kind, const std::string &message, std::string traceback)
    : std::runtime_error(message), kind_(std::move(kind)), traceback_(std::move(traceback)) {
}

std::vector<std::pair<std::string, int64_t>> MiintCapabilities() {
	return {{"ancombc", 1}, {"mantel", 1}};
}

std::string MiintSkbArrowCacheBinary() {
	const std::string dir = gb::MiintGplBoundaryCacheDir();
	return dir.empty() ? std::string() : dir + "/skb-arrow";
}

std::string FindSkbArrow() {
	if (const char *override_path = ::getenv("MIINT_SKB_ARROW_PATH"); override_path && override_path[0] != '\0') {
		return override_path;
	}
	std::string cached = MiintSkbArrowCacheBinary();
	if (!cached.empty() && ::access(cached.c_str(), X_OK) == 0) {
		return cached;
	}
	return gb::FindExecutableInPath("skb-arrow");
}

std::string SegmentDirParent() {
	if (const char *dir = ::getenv("MIINT_SKB_ARROW_SEGMENT_DIR"); dir && dir[0] != '\0') {
		return dir;
	}
#if defined(__linux__)
	struct stat st;
	if (::stat("/dev/shm", &st) == 0 && S_ISDIR(st.st_mode) && ::access("/dev/shm", W_OK | X_OK) == 0) {
		return "/dev/shm";
	}
#endif
	if (const char *tmp = ::getenv("TMPDIR"); tmp && tmp[0] != '\0') {
		return tmp;
	}
	return "/tmp";
}

// =============================================================================
// SessionDir
// =============================================================================

SessionDir::SessionDir(const std::string &parent) {
	std::string tmpl = parent + "/miint-skb-XXXXXX";
	if (!::mkdtemp(tmpl.data())) {
		const int err = errno;
		throw Error("resource", "cannot create an skb-arrow session directory in '" + parent + "': " + Errno(err) +
		                            " (the parent comes from MIINT_SKB_ARROW_SEGMENT_DIR, else /dev/shm or TMPDIR)");
	}
	path_ = tmpl;
}

SessionDir::~SessionDir() {
	Remove();
}

void SessionDir::Remove() noexcept {
	// List first, then unlink: removing entries while readdir walks the
	// directory can skip some (it does on macOS).
	std::vector<std::string> names;
	if (DIR *d = ::opendir(path_.c_str())) {
		while (const dirent *e = ::readdir(d)) {
			if (std::strcmp(e->d_name, ".") != 0 && std::strcmp(e->d_name, "..") != 0) {
				names.emplace_back(e->d_name);
			}
		}
		::closedir(d);
	}
	for (const auto &name : names) {
		(void)::unlink((path_ + "/" + name).c_str());
	}
	(void)::rmdir(path_.c_str());
}

// =============================================================================
// SegmentWriter
// =============================================================================

struct SegmentWriter::Impl {
	std::string dir;
	std::string prefix;
	std::size_t cap;
	ArrowSchema schema {};
	std::vector<std::string> names; // segments this writer created
	std::string path;               // the open segment
	int fd = -1;
	ArrowIpcWriter writer {};
	bool writer_live = false;
	std::size_t written = 0; // bytes the open segment's writer has written
	int write_errno = 0;
	bool finished = false;

	~Impl() {
		CloseQuietly();
		if (!finished) {
			for (const auto &name : names) {
				(void)::unlink((dir + "/" + name).c_str());
			}
		}
		if (schema.release) {
			schema.release(&schema);
		}
	}

	// An ArrowIpcOutputStream over the open segment's fd, via write(2).
	static ArrowErrorCode Write(ArrowIpcOutputStream *stream, const void *buf, int64_t size, int64_t *written_out,
	                            ArrowError *error) {
		auto &self = *static_cast<Impl *>(stream->private_data);
		const char *p = static_cast<const char *>(buf);
		int64_t left = size;
		while (left > 0) {
			const ssize_t n = ::write(self.fd, p, static_cast<std::size_t>(left));
			if (n < 0) {
				const int err = errno;
				if (err == EINTR) {
					continue;
				}
				self.write_errno = err;
				ArrowErrorSet(error, "%s", ::strerror(err));
				return err;
			}
			p += n;
			left -= n;
		}
		*written_out = size;
		self.written += static_cast<std::size_t>(size);
		return NANOARROW_OK;
	}
	static void Release(ArrowIpcOutputStream *stream) {
		stream->release = nullptr; // the fd belongs to Impl
	}

	[[noreturn]] void Fail(const std::string &what, const ArrowError &err) {
		if (write_errno != 0) {
			throw Error(IsResourceErrno(write_errno) ? "resource" : "internal",
			            "cannot write skb-arrow segment " + path + ": " + Errno(write_errno));
		}
		throw Error("internal", "cannot encode skb-arrow segment " + path + " (" + what + "): " + err.message);
	}

	void Open() {
		const std::string name = prefix + "-" + std::to_string(names.size());
		if (!IsSegmentName(name)) {
			throw Error("internal", "'" + name + "' is not a valid skb-arrow segment name");
		}
		path = dir + "/" + name;
		fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC, 0600);
		if (fd < 0) {
			const int err = errno;
			throw Error(IsResourceErrno(err) ? "resource" : "internal",
			            "cannot create skb-arrow segment " + path + ": " + Errno(err));
		}
		names.push_back(name);
		ArrowIpcOutputStream stream {};
		stream.write = &Impl::Write;
		stream.release = &Impl::Release;
		stream.private_data = this;
		if (ArrowIpcWriterInit(&writer, &stream) != NANOARROW_OK) {
			if (stream.release) {
				stream.release(&stream);
			}
			throw Error("internal", "cannot start an Arrow IPC writer for " + path);
		}
		writer_live = true;
		ArrowError err {};
		if (ArrowIpcWriterWriteSchema(&writer, &schema, &err) != NANOARROW_OK) {
			Fail("schema", err);
		}
		written = 0; // the cap counts batches, not the schema
	}

	void Close() {
		ArrowError err {};
		if (ArrowIpcWriterWriteArrayView(&writer, nullptr, &err) != NANOARROW_OK) {
			Fail("end of stream", err);
		}
		ArrowIpcWriterReset(&writer);
		writer_live = false;
		const int rc = ::close(fd);
		fd = -1;
		if (rc != 0) {
			const int close_err = errno;
			throw Error(IsResourceErrno(close_err) ? "resource" : "internal",
			            "cannot write skb-arrow segment " + path + ": " + Errno(close_err));
		}
	}

	void CloseQuietly() noexcept {
		if (writer_live) {
			ArrowIpcWriterReset(&writer);
			writer_live = false;
		}
		if (fd >= 0) {
			::close(fd);
			fd = -1;
		}
	}
};

SegmentWriter::SegmentWriter(std::string dir, std::string prefix, const ArrowSchema &schema, std::size_t cap)
    : impl_(std::make_unique<Impl>()) {
	impl_->dir = std::move(dir);
	impl_->prefix = std::move(prefix);
	impl_->cap = cap;
	if (ArrowSchemaDeepCopy(&schema, &impl_->schema) != NANOARROW_OK) {
		throw Error("resource", "cannot copy an Arrow schema for an skb-arrow segment");
	}
}

SegmentWriter::~SegmentWriter() = default;

void SegmentWriter::Append(const ArrowArray &batch) {
	Impl &s = *impl_;
	if (s.writer_live && s.written >= s.cap) {
		s.Close();
	}
	if (!s.writer_live) {
		s.Open();
	}
	ArrowArrayView view;
	std::memset(&view, 0, sizeof(view));
	ArrowError err {};
	ArrowErrorCode rc = ArrowArrayViewInitFromSchema(&view, &s.schema, &err);
	if (rc == NANOARROW_OK) {
		rc = ArrowArrayViewSetArray(&view, &batch, &err);
	}
	if (rc == NANOARROW_OK) {
		rc = ArrowIpcWriterWriteArrayView(&s.writer, &view, &err);
	}
	ArrowArrayViewReset(&view);
	if (rc != NANOARROW_OK) {
		s.Fail("batch", err);
	}
}

std::vector<std::string> SegmentWriter::Finish() {
	Impl &s = *impl_;
	if (!s.writer_live) {
		s.Open(); // no batches: one schema-only segment, an empty table
	}
	s.Close();
	s.finished = true;
	return s.names;
}

// =============================================================================
// Params
// =============================================================================

void Params::SetInt(const std::string &name, int64_t value) {
	entries_.emplace_back(name, value);
}

void Params::SetNumber(const std::string &name, double value) {
	if (!std::isfinite(value)) {
		throw Error("invalid_param",
		            "skb-arrow param '" + name + "' must be a finite number, not " + std::to_string(value));
	}
	entries_.emplace_back(name, value);
}

void Params::SetBool(const std::string &name, bool value) {
	entries_.emplace_back(name, value);
}

void Params::SetString(const std::string &name, std::string value) {
	entries_.emplace_back(name, std::move(value));
}

void Params::SetStrings(const std::string &name, std::vector<std::string> value) {
	entries_.emplace_back(name, std::move(value));
}

// =============================================================================
// MappedSegment
// =============================================================================

MappedSegment::MappedSegment(void *data, std::size_t size) : data_(data), size_(size) {
}

MappedSegment::~MappedSegment() {
	if (data_) {
		::munmap(data_, size_);
	}
}

MappedSegment::MappedSegment(MappedSegment &&other) noexcept : data_(other.data_), size_(other.size_) {
	other.data_ = nullptr;
	other.size_ = 0;
}

MappedSegment &MappedSegment::operator=(MappedSegment &&other) noexcept {
	if (this != &other) {
		this->~MappedSegment();
		new (this) MappedSegment(std::move(other));
	}
	return *this;
}

// =============================================================================
// Host
// =============================================================================

Host::Host(const HostOptions &options, const InterruptCheck &check)
    : options_(options), dir_(options.segment_parent), child_(SpawnHost(options, dir_.path())),
      reader_(child_.stdout_fd()) {
	try {
		const std::string line = Exchange(R"({"type":"init","protocol_version":1})", check, "the handshake");
		const DocPtr doc = ParseResponse(line);
		yj::yyjson_val *root = yj::yyjson_doc_get_root(doc.get());
		const std::string type = GetString(root, "type");
		if (type == "error") {
			throw ErrorFromResponse(root);
		}
		if (type != "ready") {
			ThrowMalformed("expected 'ready' to init, got '" + type + "'");
		}
		yj::yyjson_val *version = yj::yyjson_obj_get(root, "protocol_version");
		if (!version || !yj::yyjson_is_int(version) || yj::yyjson_get_sint(version) != kProtocolVersion) {
			ThrowMalformed("'ready' does not carry protocol_version " + std::to_string(kProtocolVersion));
		}
		bool has_version = false;
		host_version_ = GetString(root, "host_version", &has_version);
		yj::yyjson_val *capabilities = yj::yyjson_obj_get(root, "capabilities");
		if (!has_version || !capabilities || !yj::yyjson_is_obj(capabilities)) {
			ThrowMalformed("'ready' lacks a string 'host_version' or an object 'capabilities'");
		}
		const std::string fix = ". Install the version miint is tested with: SELECT install_skb_arrow();";
		for (const auto &required : options_.required) {
			const std::string &name = required.first;
			yj::yyjson_val *v = yj::yyjson_obj_get(capabilities, name.c_str());
			const std::string which = "skb-arrow " + host_version_ + " at " + options_.binary;
			if (!v || !yj::yyjson_is_int(v)) {
				throw Error("host_incompatible", which + " has no capability '" + name +
				                                     "', which miint needs at schema_version " +
				                                     std::to_string(required.second) + fix);
			}
			if (yj::yyjson_get_sint(v) != required.second) {
				throw Error("host_incompatible", which + " offers " + name + " at schema_version " +
				                                     std::to_string(yj::yyjson_get_sint(v)) +
				                                     ", but miint was built against schema_version " +
				                                     std::to_string(required.second) + fix);
			}
		}
	} catch (...) {
		Kill();
		throw;
	}
}

Host::~Host() {
	Kill();
}

CallResult Host::Call(const std::string &capability, const Params &params,
                      const std::vector<std::pair<std::string, std::vector<std::string>>> &input,
                      const InterruptCheck &check) {
	if (done_) {
		throw Error("internal", "skb-arrow host is gone; start a new one");
	}
	const int64_t id = next_id_++;

	// {"type":"call","id":…,"capability":…,"params":{…},"input":{table:[segments]}}
	MutDocPtr doc(yj::yyjson_mut_doc_new(nullptr), &yj::yyjson_mut_doc_free);
	yj::yyjson_mut_doc *d = doc.get();
	yj::yyjson_mut_val *root = yj::yyjson_mut_obj(d);
	yj::yyjson_mut_doc_set_root(d, root);
	yj::yyjson_mut_obj_add_str(d, root, "type", "call");
	yj::yyjson_mut_obj_add_sint(d, root, "id", id);
	yj::yyjson_mut_obj_add_strncpy(d, root, "capability", capability.data(), capability.size());
	if (!params.entries().empty()) {
		yj::yyjson_mut_val *p = yj::yyjson_mut_obj_add_obj(d, root, "params");
		for (const auto &entry : params.entries()) {
			yj::yyjson_mut_val *value = std::visit(
			    [d](const auto &v) -> yj::yyjson_mut_val * {
				    using T = std::decay_t<decltype(v)>;
				    if constexpr (std::is_same_v<T, int64_t>) {
					    return yj::yyjson_mut_sint(d, v);
				    } else if constexpr (std::is_same_v<T, double>) {
					    return yj::yyjson_mut_real(d, v);
				    } else if constexpr (std::is_same_v<T, bool>) {
					    return yj::yyjson_mut_bool(d, v);
				    } else if constexpr (std::is_same_v<T, std::string>) {
					    return yj::yyjson_mut_strncpy(d, v.data(), v.size());
				    } else {
					    yj::yyjson_mut_val *arr = yj::yyjson_mut_arr(d);
					    for (const auto &s : v) {
						    yj::yyjson_mut_arr_add_strncpy(d, arr, s.data(), s.size());
					    }
					    return arr;
				    }
			    },
			    entry.second);
			yj::yyjson_mut_obj_add(p, yj::yyjson_mut_strncpy(d, entry.first.data(), entry.first.size()), value);
		}
	}
	yj::yyjson_mut_val *in = yj::yyjson_mut_obj_add_obj(d, root, "input");
	for (const auto &table : input) {
		yj::yyjson_mut_val *segments = yj::yyjson_mut_arr(d);
		for (const auto &segment : table.second) {
			yj::yyjson_mut_arr_add_strncpy(d, segments, segment.data(), segment.size());
		}
		yj::yyjson_mut_obj_add(in, yj::yyjson_mut_strncpy(d, table.first.data(), table.first.size()), segments);
	}
	std::size_t len = 0;
	char *json = yj::yyjson_mut_write(d, 0, &len);
	if (!json) {
		throw Error("internal", "cannot encode the skb-arrow call to " + capability);
	}
	const std::string request(json, len);
	std::free(json);

	std::string line;
	try {
		line = Exchange(request, check, "a call");
	} catch (...) {
		Kill();
		throw;
	}

	// A response that breaks the protocol means the host can't be trusted
	// with another call; an error response (thrown below) leaves it serving.
	DocPtr response(nullptr, &yj::yyjson_doc_free);
	yj::yyjson_val *r = nullptr;
	try {
		response = ParseResponse(line);
		r = yj::yyjson_doc_get_root(response.get());
		yj::yyjson_val *rid = yj::yyjson_obj_get(r, "id");
		if (!rid || !yj::yyjson_is_int(rid) || yj::yyjson_get_sint(rid) != id) {
			ThrowMalformed("expected the response to call id " + std::to_string(id));
		}
		const std::string type = GetString(r, "type");
		if (type == "result") {
			return ParseResult(r);
		}
		if (type != "error") {
			ThrowMalformed("expected 'result' or 'error', got '" + type + "'");
		}
	} catch (...) {
		Kill();
		throw;
	}
	throw ErrorFromResponse(r);
}

MappedSegment Host::MapOutput(const std::string &name) {
	if (!IsSegmentName(name)) {
		throw Error("host_incompatible", "skb-arrow named an invalid output segment '" + name + "'");
	}
	const std::string path = dir_.path() + "/" + name;
	const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
	if (fd < 0) {
		const int err = errno;
		throw Error("host_incompatible", "skb-arrow output segment " + path + " cannot be opened: " + Errno(err));
	}
	struct stat st;
	if (::fstat(fd, &st) != 0 || st.st_size <= 0) {
		::close(fd);
		::unlink(path.c_str());
		throw Error("host_incompatible", "skb-arrow output segment " + path + " is empty");
	}
	const auto size = static_cast<std::size_t>(st.st_size);
	void *p = ::mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0);
	const int map_err = errno;
	::close(fd);
	// Receiver disposes: the mapping keeps the bytes; the name is free again.
	::unlink(path.c_str());
	if (p == MAP_FAILED) {
		throw Error(IsResourceErrno(map_err) ? "resource" : "internal",
		            "cannot map skb-arrow output segment " + path + ": " + Errno(map_err));
	}
	return MappedSegment(p, size);
}

void Host::Shutdown() {
	if (done_) {
		return;
	}
	child_.CloseStdin();
	(void)child_.ReapAndDescribe(options_.shutdown_grace_ms);
	Kill(); // a host still running after the grace period
}

void Host::Kill() noexcept {
	if (done_) {
		return;
	}
	done_ = true;
	try {
		child_.KillGroup();
	} catch (...) { // only without new_session, which SpawnHost always sets
	}
	child_.Wait();
	dir_.Remove();
}

std::string Host::Exchange(const std::string &request, const InterruptCheck &check, const char *awaiting) {
	try {
		gb::WriteLine(child_.stdin_fd(), request);
	} catch (const std::runtime_error &) {
		ThrowHostDeath(awaiting);
	}
	while (!reader_.HasLine()) {
		pollfd p {child_.stdout_fd(), POLLIN, 0};
		const int ready = ::poll(&p, 1, kPollMs);
		if (ready < 0 && errno != EINTR) {
			const int err = errno;
			throw Error("internal", "cannot wait on skb-arrow: poll: " + Errno(err));
		}
		if (ready > 0) {
			reader_.FillOnce();
		}
		if (check) {
			check();
		}
	}
	// At EOF without a whole line, the host died, perhaps while writing: a
	// partial line is not a response.
	if (!reader_.HasCompleteLine()) {
		ThrowHostDeath(awaiting);
	}
	std::string line;
	(void)reader_.ReadLine(line);
	return line;
}

void Host::ThrowHostDeath(const char *awaiting) {
	// stderr first: its drainer relays the host's last words, then exits.
	// Killing the group before reading them would lose them.
	std::string tail = ReadToEof(child_.stderr_fd(), kStderrEofMs);
	gb::TruncateHead(tail, kStderrTailCap);
	child_.KillGroup(); // anything the host left running
	const std::string how = child_.ReapAndDescribe(1000);
	Kill();
	std::string message = std::string("skb-arrow exited during ") + awaiting + " (" + how + ")";
	if (!tail.empty()) {
		message += "\n--- skb-arrow stderr ---\n" + tail;
	}
	throw Error("internal", message);
}

} // namespace skb_arrow
} // namespace miint
} // namespace duckdb
