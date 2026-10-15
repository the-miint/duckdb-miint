#include "aspera_stream.hpp"

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)

#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include <algorithm>
#include <cctype>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __linux__
#include <sys/prctl.h>
#endif

namespace miint {

// ===========================================================================
// AsperaProcess
// ===========================================================================

AsperaProcess::AsperaProcess(const AsperaConfig &config, const std::vector<std::string> &remote_paths)
    : child_pid_(-1), pipe_fd_(-1), remaining_(0), tar_mode_(remote_paths.size() > 1), first_file_read_(false),
      pipe_eof_(false), stderr_fd_(-1), has_peek_(false), peek_byte_(0) {
	if (remote_paths.empty()) {
		throw duckdb::IOException("AsperaProcess: no remote paths provided");
	}

	// Build argv
	std::vector<std::string> args_storage;
	args_storage.reserve(20 + remote_paths.size());
	std::vector<char *> args;

	auto push_arg = [&](const std::string &arg) {
		args_storage.push_back(arg);
		args.push_back(const_cast<char *>(args_storage.back().c_str()));
	};

	push_arg(config.ascp_path);
	// No progress or session status on stdout, which carries the data under
	// stdio://. Without -q, a failed login still writes "Session Stop (Error:
	// ...)" there, which reads as the start of the file.
	push_arg("-q");
	push_arg("-QT"); // Adaptive rate + disable encryption
	push_arg("-l");  // Max transfer rate
	push_arg(config.max_rate);
	push_arg("-P"); // SSH port
	push_arg(std::to_string(config.port));
	push_arg("--mode=recv");
	push_arg("--user=" + config.user);
	push_arg("--host=" + config.host);
	push_arg("-i"); // SSH key
	push_arg(config.key_path);

	// For large file lists, write to a temp file to avoid exceeding argv limits
	size_t total_path_chars = 0;
	for (const auto &p : remote_paths) {
		total_path_chars += p.size();
	}

	if (total_path_chars > 100000) {
		// Write paths to temp file
		std::string tmpl_str = miint::GetTempDir() + "/miint_aspera_XXXXXX";
		std::vector<char> tmpl_buf(tmpl_str.begin(), tmpl_str.end());
		tmpl_buf.push_back('\0');
		int fd = mkstemp(tmpl_buf.data());
		if (fd == -1) {
			throw duckdb::IOException("AsperaProcess: failed to create temp file for --file-list");
		}
		file_list_path_ = std::string(tmpl_buf.data());
		for (const auto &p : remote_paths) {
			auto line = p + "\n";
			auto written = write(fd, line.data(), line.size());
			(void)written;
		}
		close(fd);
		push_arg("--file-list=" + file_list_path_);
	} else {
		for (const auto &p : remote_paths) {
			push_arg(p);
		}
	}

	// Target: stdio:// for single file, stdio-tar:// for multiple
	push_arg(tar_mode_ ? "stdio-tar://" : "stdio://");

	args.push_back(nullptr);

	// ascp's stderr goes to an unlinked temp file rather than a pipe: nothing
	// reads it while the transfer runs, so a pipe could fill and stall ascp, and
	// the file still holds the reason when ascp fails. Without the file, stderr
	// goes to /dev/null. Made before the pipe, so nothing between pipe() and
	// fork() can throw.
	{
		std::string tmpl_str = miint::GetTempDir() + "/miint_ascp_stderr_XXXXXX";
		std::vector<char> tmpl_buf(tmpl_str.begin(), tmpl_str.end());
		tmpl_buf.push_back('\0');
		stderr_fd_ = mkostemp(tmpl_buf.data(), O_CLOEXEC);
		if (stderr_fd_ != -1) {
			unlink(tmpl_buf.data());
		}
	}

	// Every descriptor made here is close-on-exec, so only this child gets them,
	// through the dup2s below, not a process another thread forks meanwhile. A
	// stray copy of the pipe's write end would keep it open after ascp exits,
	// and AwaitData would wait on that process instead. pipe2 sets it at
	// creation; without pipe2 (macOS) a fork between pipe() and fcntl() can
	// still copy the ends.
	int stdout_pipe[2];
#ifdef __linux__
	int pipe_rc = pipe2(stdout_pipe, O_CLOEXEC);
#else
	int pipe_rc = pipe(stdout_pipe);
	if (pipe_rc == 0) {
		fcntl(stdout_pipe[0], F_SETFD, FD_CLOEXEC);
		fcntl(stdout_pipe[1], F_SETFD, FD_CLOEXEC);
	}
#endif
	if (pipe_rc == -1) {
		int err = errno;
		if (stderr_fd_ != -1) {
			close(stderr_fd_);
			stderr_fd_ = -1;
		}
		if (!file_list_path_.empty()) {
			unlink(file_list_path_.c_str());
		}
		throw duckdb::IOException("AsperaProcess: failed to create pipe: %s", strerror(err));
	}

	child_pid_ = fork();
	if (child_pid_ == -1) {
		int err = errno;
		close(stdout_pipe[0]);
		close(stdout_pipe[1]);
		if (stderr_fd_ != -1) {
			close(stderr_fd_);
		}
		if (!file_list_path_.empty()) {
			unlink(file_list_path_.c_str());
		}
		throw duckdb::IOException("AsperaProcess: failed to fork: %s", strerror(err));
	}

	if (child_pid_ == 0) {
		// Child process
#ifdef __linux__
		prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif
		// When the server rejects the key, ascp asks for a password: on its
		// controlling terminal if it has one, else on stdin, and waits for an
		// answer. A new session has no controlling terminal, and stdin is
		// /dev/null, never closed (a closed fd 0 goes to whatever ascp opens
		// next), so the prompt sees end of input and ascp exits.
		setsid();

		// This process may itself run with fd 0, 1 or 2 closed, which puts a pipe
		// end, the stderr file or /dev/null on one of them. So each source is
		// first copied above fd 2, close-on-exec: no dup2 below then overwrites a
		// source a later one needs, and nothing closes a standard descriptor.
		int devnull = open("/dev/null", O_RDWR | O_CLOEXEC);
		int in_fd = devnull == -1 ? -1 : fcntl(devnull, F_DUPFD_CLOEXEC, STDERR_FILENO + 1);
		int out_fd = fcntl(stdout_pipe[1], F_DUPFD_CLOEXEC, STDERR_FILENO + 1);
		int err_src = stderr_fd_ != -1 ? stderr_fd_ : devnull;
		int err_fd = err_src == -1 ? -1 : fcntl(err_src, F_DUPFD_CLOEXEC, STDERR_FILENO + 1);
		if (in_fd == -1 || out_fd == -1 || err_fd == -1) {
			_exit(127);
		}
		if (dup2(in_fd, STDIN_FILENO) == -1 || dup2(out_fd, STDOUT_FILENO) == -1 || dup2(err_fd, STDERR_FILENO) == -1) {
			_exit(127);
		}

		execv(config.ascp_path.c_str(), args.data());
		_exit(127);
	}

	// Parent
	close(stdout_pipe[1]);
	pipe_fd_ = stdout_pipe[0];
}

AsperaProcess::~AsperaProcess() {
	if (pipe_fd_ >= 0) {
		close(pipe_fd_);
		pipe_fd_ = -1;
	}
	if (stderr_fd_ >= 0) {
		close(stderr_fd_);
		stderr_fd_ = -1;
	}
	if (child_pid_ > 0) {
		kill(child_pid_, SIGTERM);
		// Give ascp a moment to exit cleanly, then force-kill.
		// ascp catches SIGTERM and may block during network cleanup.
		int status;
		for (int i = 0; i < 10; i++) {
			pid_t ret = waitpid(child_pid_, &status, WNOHANG);
			if (ret > 0) {
				goto reaped;
			}
			usleep(100000); // 100ms
		}
		// Still alive after 1s — force kill
		kill(child_pid_, SIGKILL);
		waitpid(child_pid_, &status, 0);
	reaped:
		child_pid_ = -1;
	}
	if (!file_list_path_.empty()) {
		unlink(file_list_path_.c_str());
	}
}

ssize_t AsperaProcess::ReadRaw(void *dst, size_t n) {
	if (has_peek_ && n > 0) {
		static_cast<char *>(dst)[0] = peek_byte_;
		has_peek_ = false;
		return 1;
	}
	ssize_t r;
	do {
		r = read(pipe_fd_, dst, n);
	} while (r == -1 && errno == EINTR);
	return r;
}

bool AsperaProcess::ReadByte(char &c) {
	if (pipe_eof_) {
		return false;
	}
	ssize_t n = ReadRaw(&c, 1);
	if (n <= 0) {
		pipe_eof_ = true;
		return false;
	}
	return true;
}

bool AsperaProcess::ReadExact(void *dst, size_t n) {
	auto *buf = static_cast<char *>(dst);
	size_t total = 0;
	while (total < n) {
		ssize_t r = ReadRaw(buf + total, n - total);
		if (r <= 0) {
			pipe_eof_ = true;
			return false;
		}
		total += static_cast<size_t>(r);
	}
	return true;
}

bool AsperaProcess::ReadUntil(std::string &out, char delimiter) {
	out.clear();
	char c;
	while (ReadByte(c)) {
		if (c == delimiter) {
			return true;
		}
		out += c;
	}
	return false;
}

bool AsperaProcess::NextFile(std::string &out_filename, size_t &out_size) {
	if (pipe_eof_) {
		return false;
	}

	if (!tar_mode_) {
		// stdio:// mode: single file, no framing
		if (first_file_read_) {
			return false;
		}
		first_file_read_ = true;
		out_filename = "";
		out_size = 0; // Unknown size in stdio:// mode
		remaining_ = SIZE_MAX;
		return true;
	}

	// stdio-tar:// mode: parse framing header
	// Format: \nFile: <filename>\0\nSize: <decimal_bytes>\n\n

	// Read leading \n
	char c;
	if (!ReadByte(c)) {
		return false;
	}
	if (c != '\n') {
		// Unexpected byte — could be end of stream or corruption
		return false;
	}

	// Read "File: "
	std::string prefix;
	for (int i = 0; i < 6; i++) {
		if (!ReadByte(c)) {
			return false;
		}
		prefix += c;
	}
	if (prefix != "File: ") {
		return false;
	}

	// Read filename until NUL
	if (!ReadUntil(out_filename, '\0')) {
		return false;
	}

	// Read \nSize:
	std::string size_prefix;
	for (int i = 0; i < 7; i++) {
		if (!ReadByte(c)) {
			return false;
		}
		size_prefix += c;
	}
	if (size_prefix != "\nSize: ") {
		return false;
	}

	// Read size until \n
	std::string size_str;
	if (!ReadUntil(size_str, '\n')) {
		return false;
	}

	// Read the trailing \n (empty line before data)
	if (!ReadByte(c) || c != '\n') {
		return false;
	}

	try {
		out_size = std::stoull(size_str);
	} catch (const std::exception &) {
		throw duckdb::IOException("AsperaProcess: malformed stdio-tar size header: '%s'", size_str);
	}
	remaining_ = out_size;
	first_file_read_ = true;
	return true;
}

int AsperaProcess::ReadBounded(void *dst, size_t len) {
	if (remaining_ == 0) {
		return 0;
	}
	if (pipe_eof_) {
		return 0;
	}

	size_t to_read = len;
	if (tar_mode_ && remaining_ != SIZE_MAX) {
		to_read = std::min(len, remaining_);
	}

	ssize_t n = ReadRaw(dst, to_read);
	if (n <= 0) {
		pipe_eof_ = true;
		return (n == 0) ? 0 : -1;
	}

	if (tar_mode_ && remaining_ != SIZE_MAX) {
		remaining_ -= static_cast<size_t>(n);
	}
	return static_cast<int>(n);
}

void AsperaProcess::ReadCurrentFileToBuffer(std::vector<char> &out) {
	if (tar_mode_ && remaining_ != SIZE_MAX) {
		out.resize(remaining_);
		size_t total = 0;
		while (total < out.size()) {
			int n = ReadBounded(out.data() + total, out.size() - total);
			if (n <= 0) {
				out.resize(total);
				return;
			}
			total += static_cast<size_t>(n);
		}
	} else {
		// stdio:// mode: read until EOF
		out.clear();
		char buf[65536];
		while (true) {
			int n = ReadBounded(buf, sizeof(buf));
			if (n <= 0) {
				break;
			}
			out.insert(out.end(), buf, buf + n);
		}
	}
}

int AsperaProcess::WaitForExit() {
	if (child_pid_ <= 0) {
		return -1;
	}

	// Drain remaining pipe data to avoid blocking the child
	if (pipe_fd_ >= 0) {
		char discard[4096];
		ssize_t n;
		while ((n = read(pipe_fd_, discard, sizeof(discard))) > 0 || (n == -1 && errno == EINTR)) {
		}
		close(pipe_fd_);
		pipe_fd_ = -1;
	}

	int status = 0;
	pid_t reaped;
	do {
		reaped = waitpid(child_pid_, &status, 0);
	} while (reaped == -1 && errno == EINTR);
	child_pid_ = -1;

	if (reaped != -1 && WIFEXITED(status)) {
		return WEXITSTATUS(status);
	}
	return -1;
}

bool AsperaProcess::AwaitData() {
	if (has_peek_) {
		return true;
	}
	char c;
	if (!ReadByte(c)) {
		return false;
	}
	peek_byte_ = c;
	has_peek_ = true;
	return true;
}

std::string AsperaProcess::StderrText() const {
	if (stderr_fd_ < 0) {
		return "";
	}
	// The reason ascp failed is at the end, so keep the last 4 KiB.
	constexpr off_t KEEP = 4096;
	struct stat st;
	if (fstat(stderr_fd_, &st) != 0) {
		return "";
	}
	off_t offset = st.st_size > KEEP ? st.st_size - KEEP : 0;
	const bool truncated = offset > 0;
	std::string raw;
	char buf[4096];
	ssize_t n;
	while ((n = pread(stderr_fd_, buf, sizeof(buf), offset)) > 0) {
		raw.append(buf, static_cast<size_t>(n));
		offset += n;
	}
	std::string text = truncated ? "[...] " : "";
	for (char ch : raw) {
		if (std::isspace(static_cast<unsigned char>(ch))) {
			if (!text.empty() && text.back() != ' ') {
				text += ' ';
			}
		} else {
			text += ch;
		}
	}
	if (!text.empty() && text.back() == ' ') {
		text.pop_back();
	}
	return text;
}

// ===========================================================================
// AsperaSeqStream
// ===========================================================================

AsperaSeqStream::AsperaSeqStream()
    : process(nullptr), is_gzipped(false), zs({}), zs_initialized(false), compressed_avail(0), compressed_next(nullptr),
      input_eof(false), stream_end(false) {
}

AsperaSeqStream::~AsperaSeqStream() {
	if (zs_initialized) {
		inflateEnd(&zs);
	}
}

int aspera_seq_read(AsperaSeqStream *stream, void *dst, unsigned int len) {
	// On error: set the thread-local error channel and return -1; SequenceReader
	// picks it up after kseq++ returns and raises IOException. See
	// g_seq_read_error docs in duckdb_seq_stream.hpp.
	try {
		if (!stream->is_gzipped) {
			int n = stream->process->ReadBounded(dst, len);
			if (n < 0) {
				g_seq_read_error = "aspera_seq_read: ascp pipe read failed";
				return -1;
			}
			return n;
		}

		auto read_raw = [stream](void *buf, size_t sz) -> int {
			return stream->process->ReadBounded(buf, sz);
		};

		int result = InflateFromSource(stream->zs, stream->compressed_buf, AsperaSeqStream::COMPRESSED_BUF_SIZE,
		                               stream->compressed_avail, stream->compressed_next, stream->input_eof,
		                               stream->stream_end, read_raw, dst, len);
		if (result < 0) {
			g_seq_read_error = stream->stream_end ? "aspera_seq_read: gzip decompression error"
			                                      : "aspera_seq_read: gzip stream truncated before end marker "
			                                        "(ascp transfer likely incomplete)";
			return -1;
		}
		return result;
	} catch (const std::exception &e) {
		g_seq_read_error = std::string("aspera_seq_read: ") + e.what();
		return -1;
	} catch (...) {
		g_seq_read_error = "aspera_seq_read: unknown error";
		return -1;
	}
}

int aspera_seq_close(AsperaSeqStream *stream) {
	delete stream;
	return 0;
}

AsperaSeqStream *CreateAsperaSeqStream(AsperaProcess *process, bool is_gzipped) {
	auto stream = new AsperaSeqStream();
	stream->process = process;
	stream->is_gzipped = is_gzipped;
	if (is_gzipped) {
		if (inflateInit2(&stream->zs, 16 + MAX_WBITS) != Z_OK) {
			delete stream;
			throw duckdb::IOException("Failed to initialize zlib for Aspera stream");
		}
		stream->zs_initialized = true;
	}
	return stream;
}

} // namespace miint

#endif // MIINT_ASPERA_SUPPORTED && MIINT_STATIC_BUILD
