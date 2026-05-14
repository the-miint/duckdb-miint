#include "aspera_stream.hpp"

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)

#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include <algorithm>
#include <csignal>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <sstream>
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
      pipe_eof_(false) {
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

	// Create pipe for stdout
	int stdout_pipe[2];
	if (pipe(stdout_pipe) == -1) {
		if (!file_list_path_.empty()) {
			unlink(file_list_path_.c_str());
		}
		throw duckdb::IOException("AsperaProcess: failed to create pipe: %s", strerror(errno));
	}

	child_pid_ = fork();
	if (child_pid_ == -1) {
		close(stdout_pipe[0]);
		close(stdout_pipe[1]);
		if (!file_list_path_.empty()) {
			unlink(file_list_path_.c_str());
		}
		throw duckdb::IOException("AsperaProcess: failed to fork: %s", strerror(errno));
	}

	if (child_pid_ == 0) {
		// Child process
#ifdef __linux__
		prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif
		// Redirect stdout to pipe
		dup2(stdout_pipe[1], STDOUT_FILENO);
		close(stdout_pipe[0]);
		close(stdout_pipe[1]);

		// Redirect stderr to /dev/null
		int devnull = open("/dev/null", O_WRONLY);
		if (devnull != -1) {
			dup2(devnull, STDERR_FILENO);
			close(devnull);
		}

		// Close stdin (ascp doesn't need it for recv)
		close(STDIN_FILENO);

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

bool AsperaProcess::ReadByte(char &c) {
	if (pipe_eof_) {
		return false;
	}
	ssize_t n = read(pipe_fd_, &c, 1);
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
		ssize_t r = read(pipe_fd_, buf + total, n - total);
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
	if (pipe_eof_ || remaining_ == 0) {
		return 0;
	}

	size_t to_read = len;
	if (tar_mode_ && remaining_ != SIZE_MAX) {
		to_read = std::min(len, remaining_);
	}

	ssize_t n = read(pipe_fd_, dst, to_read);
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
		while (read(pipe_fd_, discard, sizeof(discard)) > 0) {
		}
		close(pipe_fd_);
		pipe_fd_ = -1;
	}

	int status;
	waitpid(child_pid_, &status, 0);
	child_pid_ = -1;

	if (WIFEXITED(status)) {
		return WEXITSTATUS(status);
	}
	return -1;
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
	if (!stream->is_gzipped) {
		int n = stream->process->ReadBounded(dst, len);
		if (n < 0) {
			// ascp pipe returned a hard error (not EOF). Throw so the
			// read_ena_sequences mid-stream catch records the run as truncated
			// rather than silently treating the empty batch as completion.
			throw duckdb::IOException("aspera_seq_read: ascp pipe read failed");
		}
		return n;
	}

	auto read_raw = [stream](void *buf, size_t sz) -> int {
		// ReadBounded already returns -1 on error; propagate so
		// InflateFromSource can surface it.
		return stream->process->ReadBounded(buf, sz);
	};

	int result = InflateFromSource(stream->zs, stream->compressed_buf, AsperaSeqStream::COMPRESSED_BUF_SIZE,
	                               stream->compressed_avail, stream->compressed_next, stream->input_eof,
	                               stream->stream_end, read_raw, dst, len);
	if (result < 0) {
		throw duckdb::IOException(stream->stream_end ? "aspera_seq_read: gzip decompression error"
		                                             : "aspera_seq_read: gzip stream truncated before end marker "
		                                               "(ascp transfer likely incomplete)");
	}
	return result;
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
