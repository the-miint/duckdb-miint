#include "per_run_reader.hpp"
#include "duckdb_seq_stream.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <unistd.h>

namespace miint {

PerRunReader::PerRunReader(duckdb::FileSystem &fs, ENARunInfo run, bool use_aspera, bool trim, std::mutex &open_mutex,
                           const AsperaConfig *aspera_config)
    : fs_(fs), run_(std::move(run)), use_aspera_(use_aspera), trim_(trim), open_mutex_(open_mutex),
      aspera_config_(aspera_config) {
}

PerRunReader::~PerRunReader() {
	// Readers / processes are auto-released by their unique_ptrs.
	// Temp files: best-effort cleanup on the way out.
	if (!sff_temp_path_.empty()) {
		std::remove(sff_temp_path_.c_str());
		sff_temp_path_.clear();
	}
#if MIINT_ASPERA_SUPPORTED
	if (!aspera_temp_path_.empty()) {
		std::remove(aspera_temp_path_.c_str());
		aspera_temp_path_.clear();
	}
	// AsperaProcess destructor handles SIGTERM → SIGKILL → waitpid
#endif
}

void PerRunReader::Open() {
	if (opened_) {
		return;
	}
	if (run_.format == ENASequenceFormat::SFF) {
		OpenSFF();
	} else {
#if MIINT_ASPERA_SUPPORTED
		if (use_aspera_) {
			OpenAspera();
		} else {
			OpenHTTP();
		}
#else
		OpenHTTP();
#endif
	}
	opened_ = true;
}

SequenceRecordBatch PerRunReader::ReadBatch(size_t max_size) {
	if (exhausted_) {
		return SequenceRecordBatch {};
	}
	SequenceRecordBatch batch;
	if (run_.format == ENASequenceFormat::SFF) {
		batch = sff_reader_->read(max_size);
	} else {
		batch = fastx_reader_->read(max_size);
	}
	if (batch.empty()) {
		exhausted_ = true;
	}
	return batch;
}

void PerRunReader::Finish() {
	if (finished_) {
		return;
	}
	finished_ = true;
#if MIINT_ASPERA_SUPPORTED
	if (use_aspera_ && aspera_process_) {
		int exit_code = aspera_process_->WaitForExit();
		aspera_process_.reset();
		if (exit_code != 0 && exit_code != -1) {
			throw duckdb::IOException("read_ena_sequences: ascp exited with code %d for run '%s'", exit_code,
			                          run_.run_accession);
		}
	}
#endif
}

void PerRunReader::OpenHTTP() {
	if (run_.fastq_urls.empty()) {
		throw duckdb::IOException("read_ena_sequences: no FASTQ URLs available for run '%s'", run_.run_accession);
	}

	// Serialize file opens to avoid overwhelming remote servers with simultaneous
	// HTTP HEAD requests. Reads still proceed in parallel.
	std::lock_guard<std::mutex> open_guard(open_mutex_);

	// CreateDuckDBSeqStream returns a raw `new`'d pointer. Ownership transfers
	// to SequenceReader only on successful construction. Until that call
	// succeeds, we own the streams and MUST clean them up on any throw path —
	// SequenceReader's constructor does not deallocate on failure, so without
	// these RAII guards any exception (not just std::runtime_error) leaks.
	auto stream_deleter = [](DuckDBSeqStream *p) {
		if (p) {
			duckdb_seq_close(p);
		}
	};
	std::unique_ptr<DuckDBSeqStream, decltype(stream_deleter)> s1(
	    duckdb::CreateDuckDBSeqStream(fs_, run_.fastq_urls[0]), stream_deleter);
	std::unique_ptr<DuckDBSeqStream, decltype(stream_deleter)> s2(nullptr, stream_deleter);
	if (run_.is_paired && run_.fastq_urls.size() >= 2) {
		s2.reset(duckdb::CreateDuckDBSeqStream(fs_, run_.fastq_urls[1]));
	}

	try {
		fastx_reader_ = std::make_unique<SequenceReader>(s1.get(), s2.get(), true);
		// SequenceReader took ownership — stop tracking here.
		s1.release();
		s2.release();
	} catch (const std::runtime_error &e) {
		if (s2 && std::string(e.what()) == "Empty stream (sequence2)") {
			// ENA metadata may report PAIRED layout but the second FASTQ file
			// is actually empty (single-end deposited with PAIRED metadata).
			// The failed constructor may have consumed the streams; regardless,
			// our RAII guards will release whichever ones aren't owned by it.
			// Reopen just the first file and retry as single-end.
			s1.reset();
			s2.reset();
			std::unique_ptr<DuckDBSeqStream, decltype(stream_deleter)> s1_retry(
			    duckdb::CreateDuckDBSeqStream(fs_, run_.fastq_urls[0]), stream_deleter);
			fastx_reader_ =
			    std::make_unique<SequenceReader>(s1_retry.get(), static_cast<DuckDBSeqStream *>(nullptr), true);
			s1_retry.release();
		} else {
			// s1 / s2 guards clean up automatically on the rethrow.
			throw;
		}
	}
	// Any other exception type (e.g., IOException) propagates and RAII cleans up.
}

void PerRunReader::OpenSFF() {
	if (run_.sff_url.empty()) {
		throw duckdb::IOException("read_ena_sequences: no SFF URL available for run '%s'", run_.run_accession);
	}

	auto temp_dir = GetTempDir();
	auto temp_path = temp_dir + "/miint_ena_sff_XXXXXX";
	std::vector<char> tmpl_buf(temp_path.begin(), temp_path.end());
	tmpl_buf.push_back('\0');
	int fd = mkstemp(tmpl_buf.data());
	if (fd == -1) {
		throw duckdb::IOException("read_ena_sequences: failed to create temp file for SFF download");
	}
	temp_path = std::string(tmpl_buf.data());
	sff_temp_path_ = temp_path;

	bool fd_closed = false;
	auto cleanup_on_error = [&]() {
		if (!fd_closed) {
			close(fd);
			fd_closed = true;
		}
		std::remove(temp_path.c_str());
		sff_temp_path_.clear();
	};

	try {
		duckdb::unique_ptr<duckdb::FileHandle> source;
		{
			// Serialize only the HTTP connection setup, not the full download.
			std::lock_guard<std::mutex> open_guard(open_mutex_);
			source = fs_.OpenFile(run_.sff_url, duckdb::FileOpenFlags(duckdb::FileOpenFlags::FILE_FLAGS_READ));
		}

		static constexpr size_t DOWNLOAD_BUF_SIZE = 1048576; // 1 MB
		std::vector<char> buf(DOWNLOAD_BUF_SIZE);
		while (true) {
			auto bytes_read = source->Read(buf.data(), DOWNLOAD_BUF_SIZE);
			if (bytes_read == 0) {
				break;
			}
			size_t written = 0;
			while (written < bytes_read) {
				ssize_t w = write(fd, buf.data() + written, bytes_read - written);
				if (w < 0) {
					if (errno == EINTR) {
						continue;
					}
					int saved_errno = errno;
					cleanup_on_error();
					throw duckdb::IOException(
					    "read_ena_sequences: failed writing SFF temp file for run '%s' (errno=%d: %s)",
					    run_.run_accession, saved_errno, strerror(saved_errno));
				}
				written += static_cast<size_t>(w);
			}
		}
		close(fd);
		fd_closed = true;
	} catch (...) {
		cleanup_on_error();
		throw;
	}

	sff_reader_ = std::make_unique<SFFReader>(temp_path, trim_);
}

#if MIINT_ASPERA_SUPPORTED

void PerRunReader::OpenAspera() {
	std::vector<std::string> remote_paths;
	for (const auto &ap : run_.aspera_paths) {
		remote_paths.push_back(ap.remote_path);
	}

	AsperaConfig config = aspera_config_ ? *aspera_config_ : AsperaConfig {};
	if (!run_.aspera_paths.empty()) {
		config.host = run_.aspera_paths[0].host;
	}

	// Serialize ascp launches (fork/exec); pipe reads use thread-owned descriptors.
	{
		std::lock_guard<std::mutex> open_guard(open_mutex_);
		aspera_process_ = std::make_unique<AsperaProcess>(config, remote_paths);
	}
	auto *proc = aspera_process_.get();

	std::string filename;
	size_t file_size;
	if (!proc->NextFile(filename, file_size)) {
		throw duckdb::IOException("read_ena_sequences: Aspera stream ended unexpectedly (expected files for run '%s')",
		                          run_.run_accession);
	}

	std::string gz_hint = filename;
	if (gz_hint.empty() && !run_.aspera_paths.empty()) {
		gz_hint = run_.aspera_paths[0].remote_path;
	}
	bool is_gz = duckdb::IsGzipped(gz_hint);

	if (!run_.is_paired) {
		auto *s1 = CreateAsperaSeqStream(proc, is_gz);
		try {
			fastx_reader_ = std::make_unique<SequenceReader>(s1, static_cast<AsperaSeqStream *>(nullptr), true);
		} catch (...) {
			delete s1;
			throw;
		}
		return;
	}

	// Paired-end: stream R1 from pipe to temp file, then stream R2 live.
	static constexpr size_t COPY_BUF_SIZE = 1024 * 1024; // 1 MB

	std::string tmpl = GetTempDir() + "/miint_aspera_r1_XXXXXX";
	std::vector<char> tmpl_buf(tmpl.begin(), tmpl.end());
	tmpl_buf.push_back('\0');
	int fd = mkstemp(tmpl_buf.data());
	if (fd == -1) {
		throw duckdb::IOException("read_ena_sequences: failed to create temp file for paired-end Aspera buffering");
	}
	aspera_temp_path_ = std::string(tmpl_buf.data());

	std::vector<char> copy_buf(COPY_BUF_SIZE);
	while (true) {
		int n = proc->ReadBounded(copy_buf.data(), COPY_BUF_SIZE);
		if (n == 0) {
			break;
		}
		if (n < 0) {
			close(fd);
			std::remove(aspera_temp_path_.c_str());
			aspera_temp_path_.clear();
			throw duckdb::IOException("read_ena_sequences: pipe read error while buffering R1 for run '%s'",
			                          run_.run_accession);
		}
		size_t written = 0;
		while (written < static_cast<size_t>(n)) {
			ssize_t w = write(fd, copy_buf.data() + written, static_cast<size_t>(n) - written);
			if (w < 0) {
				if (errno == EINTR) {
					continue;
				}
				int err = errno;
				close(fd);
				std::remove(aspera_temp_path_.c_str());
				aspera_temp_path_.clear();
				throw duckdb::IOException(
				    "read_ena_sequences: failed to write temp file for paired-end Aspera buffering "
				    "(run '%s', errno=%d: %s)",
				    run_.run_accession, err, strerror(err));
			}
			written += static_cast<size_t>(w);
		}
	}
	close(fd);

	std::string filename2;
	size_t file_size2;
	if (!proc->NextFile(filename2, file_size2)) {
		std::remove(aspera_temp_path_.c_str());
		aspera_temp_path_.clear();
		throw duckdb::IOException("read_ena_sequences: Aspera stream ended unexpectedly (expected R2 for run '%s')",
		                          run_.run_accession);
	}

	std::string gz_hint2 = filename2;
	if (gz_hint2.empty() && run_.aspera_paths.size() >= 2) {
		gz_hint2 = run_.aspera_paths[1].remote_path;
	}
	bool is_gz2 = duckdb::IsGzipped(gz_hint2);

	auto *s1 = duckdb::CreateDuckDBSeqStream(fs_, aspera_temp_path_, is_gz);
	AsperaSeqStream *s2 = nullptr;
	try {
		s2 = CreateAsperaSeqStream(proc, is_gz2);
		fastx_reader_ = std::make_unique<SequenceReader>(s1, s2, true);
	} catch (...) {
		delete s2;
		delete s1;
		// On throw here the SequenceReader never took ownership of s1, so the
		// R1 temp file we already wrote will be removed by the destructor via
		// aspera_temp_path_. Clear it only after Unix-level `std::remove` runs
		// from ~PerRunReader to keep the cleanup path in one place.
		throw;
	}
}

#endif // MIINT_ASPERA_SUPPORTED

} // namespace miint
