#pragma once

#include "aspera_utils.hpp"
#include "duckdb_seq_stream.hpp"
#include <string>
#include <vector>

#if MIINT_ASPERA_SUPPORTED && defined(MIINT_STATIC_BUILD)

namespace miint {

// Manages a single ascp subprocess reading from stdio:// or stdio-tar://.
// Owns the child process lifetime and the stdout pipe fd.
class AsperaProcess {
public:
	// Launch ascp with the given config and remote paths.
	// Uses stdio:// for a single path, stdio-tar:// for multiple.
	AsperaProcess(const AsperaConfig &config, const std::vector<std::string> &remote_paths);
	~AsperaProcess();

	AsperaProcess(const AsperaProcess &) = delete;
	AsperaProcess &operator=(const AsperaProcess &) = delete;

	// Advance to the next file in the tar stream.
	// Parses "\nFile: <name>\0\nSize: <nbytes>\n\n" header.
	// Returns false if no more files (EOF on pipe).
	// For stdio:// (single file), the first call sets out_filename="" and out_size=0
	// (unknown size), and subsequent calls return false.
	bool NextFile(std::string &out_filename, size_t &out_size);

	// Read up to len bytes from the current file's data segment.
	// Returns bytes read (>0), 0 = current file exhausted, -1 = error.
	// For stdio-tar://, respects the bounded size from NextFile().
	// For stdio://, reads until pipe EOF.
	int ReadBounded(void *dst, size_t len);

	// Read ALL remaining bytes for the current file into out.
	// Used for paired-end temp file buffering of file 1.
	void ReadCurrentFileToBuffer(std::vector<char> &out);

	// Wait for child process and return exit code. -1 on signal death.
	int WaitForExit();

	// Check if the process is using stdio-tar:// (multiple files) vs stdio:// (single).
	bool IsTarMode() const {
		return tar_mode_;
	}

private:
	pid_t child_pid_;
	int pipe_fd_;                // stdout of child
	size_t remaining_;           // bytes left in current file segment (tar mode)
	bool tar_mode_;              // true = stdio-tar://, false = stdio://
	bool first_file_read_;       // has NextFile() been called yet?
	bool pipe_eof_;              // pipe reached EOF
	std::string file_list_path_; // temp file for --file-list (cleaned up in destructor)

	// Read exactly n bytes from pipe. Returns false on EOF/error before n bytes.
	bool ReadExact(void *dst, size_t n);

	// Read bytes until delimiter is found. The delimiter is consumed but not included in out.
	// Returns false on EOF/error before delimiter.
	bool ReadUntil(std::string &out, char delimiter);

	// Read a single byte. Returns false on EOF.
	bool ReadByte(char &c);
};

// Stream adapter for kseq++ that reads from an AsperaProcess
// (bounded pipe segment) with optional gzip decompression.
struct AsperaSeqStream {
	AsperaProcess *process; // Non-owning; process outlives stream
	bool is_gzipped;

	z_stream zs;
	bool zs_initialized;
	static constexpr size_t COMPRESSED_BUF_SIZE = 65536;
	char compressed_buf[COMPRESSED_BUF_SIZE];
	int compressed_avail;
	char *compressed_next;
	bool input_eof;
	bool stream_end; // Z_STREAM_END observed → subsequent reads are legitimately EOF

	AsperaSeqStream();
	~AsperaSeqStream();

	AsperaSeqStream(const AsperaSeqStream &) = delete;
	AsperaSeqStream &operator=(const AsperaSeqStream &) = delete;
};

// kseq++ callbacks
int aspera_seq_read(AsperaSeqStream *stream, void *dst, unsigned int len);
int aspera_seq_close(AsperaSeqStream *stream);

// kseq++ KStreamIn instantiation for Aspera streams
using AsperaSeqStreamIn = klibpp::KStreamIn<AsperaSeqStream *, int (*)(AsperaSeqStream *, void *, unsigned int)>;

// Owning handle for a raw AsperaSeqStream until it is transferred into a kseq++
// KStreamIn (see DuckDBSeqStreamHandle for rationale).
using AsperaSeqStreamHandle = std::unique_ptr<AsperaSeqStream, int (*)(AsperaSeqStream *)>;

// Factory: create an AsperaSeqStream wrapping the current file segment of an AsperaProcess.
// is_gzipped should be determined from the filename (e.g., ends with .gz).
// Caller takes ownership (kseq++ close callback deletes it).
AsperaSeqStream *CreateAsperaSeqStream(AsperaProcess *process, bool is_gzipped);

} // namespace miint

#endif // MIINT_ASPERA_SUPPORTED && MIINT_STATIC_BUILD
