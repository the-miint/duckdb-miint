#include "per_run_reader.hpp"
#include "duckdb_seq_stream.hpp"
#include "miint_log.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/common/string_util.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <unistd.h>

namespace miint {

PerRunReader::PerRunReader(duckdb::FileSystem &fs, ENARunInfo run, bool use_aspera, bool trim, std::mutex &open_mutex,
                           const AsperaConfig *aspera_config, uint64_t max_sequences,
                           duckdb::ClientContext *log_context)
    : fs_(fs), run_(std::move(run)), use_aspera_(use_aspera), trim_(trim), open_mutex_(open_mutex),
      aspera_config_(aspera_config), max_sequences_(max_sequences), log_context_(log_context) {
}

PerRunReader::~PerRunReader() {
	// Readers / processes are auto-released by their unique_ptrs.
	// Temp files: best-effort cleanup on the way out.
	// (Aspera: AsperaProcess dtor handles SIGTERM → SIGKILL → waitpid for
	// both aspera_process_ and aspera_process_paired_, no explicit work here.)
	if (!sff_temp_path_.empty()) {
		std::remove(sff_temp_path_.c_str());
		sff_temp_path_.clear();
	}
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
	if (max_sequences_ > 0) {
		if (sequences_emitted_ >= max_sequences_) {
			exhausted_ = true;
			return SequenceRecordBatch {};
		}
		// Clamp this call so the cap fires on an exact boundary instead of
		// spilling into the next batch. batch.size() counts pairs as 1, which
		// matches the contract documented on the ctor.
		uint64_t remaining = max_sequences_ - sequences_emitted_;
		if (remaining < max_size) {
			max_size = static_cast<size_t>(remaining);
		}
	}
	SequenceRecordBatch batch;
	if (run_.format == ENASequenceFormat::SFF) {
		batch = sff_reader_->read(max_size);
	} else {
		batch = fastx_reader_->read(max_size);
	}
	sequences_emitted_ += batch.size();
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
	if (!use_aspera_) {
		return;
	}
	// Wait on both processes. If the first throws, still wait on the second to
	// avoid leaving an orphan ascp around; re-throw after both are reaped.
	auto wait_one = [this](std::unique_ptr<AsperaProcess> &proc, const char *tag) -> std::string {
		if (!proc) {
			return {};
		}
		int exit_code = proc->WaitForExit();
		proc.reset();
		if (exit_code != 0 && exit_code != -1) {
			return duckdb::StringUtil::Format("read_ena_sequences: %s ascp exited with code %d for run '%s'", tag,
			                                  exit_code, run_.run_accession);
		}
		return {};
	};
	std::string err_r1 = wait_one(aspera_process_, "R1");
	std::string err_r2 = wait_one(aspera_process_paired_, "R2");
	// When both processes errored we surface R1 via the thrown exception;
	// log R2's message first so the user doesn't lose half the story.
	if (!err_r1.empty() && !err_r2.empty()) {
		if (log_context_) {
			miint::EmitWarning(*log_context_, err_r2);
		} else {
			duckdb::Printer::Print(err_r2);
		}
	}
	if (!err_r1.empty()) {
		throw duckdb::IOException(err_r1);
	}
	if (!err_r2.empty()) {
		throw duckdb::IOException(err_r2);
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

	// SFF format requires the entire file on disk before any record can be
	// parsed — we download-then-read rather than stream. `max_sequences` still
	// caps emitted rows downstream, but the bandwidth has already been spent
	// by the time ReadBatch() starts enforcing the cap. Warn loudly once per
	// SFF run so users don't silently assume the cap saved them a download.
	if (max_sequences_ > 0) {
		auto msg =
		    duckdb::StringUtil::Format("read_ena_sequences: WARNING: max_sequences=%llu on SFF run '%s' — SFF requires "
		                               "downloading the full file before any record can be parsed; the row cap applies "
		                               "but bandwidth savings do not",
		                               static_cast<unsigned long long>(max_sequences_), run_.run_accession);
		if (log_context_) {
			miint::EmitWarning(*log_context_, msg);
		} else {
			duckdb::Printer::Print(msg);
		}
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
	if (run_.aspera_paths.empty()) {
		throw duckdb::IOException("read_ena_sequences: no Aspera paths available for run '%s'", run_.run_accession);
	}

	AsperaConfig config = aspera_config_ ? *aspera_config_ : AsperaConfig {};
	config.host = run_.aspera_paths[0].host;

	// Single-end: one ascp process in stdio:// (single-file) mode.
	if (!run_.is_paired) {
		{
			std::lock_guard<std::mutex> open_guard(open_mutex_);
			aspera_process_ =
			    std::make_unique<AsperaProcess>(config, std::vector<std::string> {run_.aspera_paths[0].remote_path});
		}
		std::string filename;
		size_t file_size;
		if (!aspera_process_->NextFile(filename, file_size)) {
			throw duckdb::IOException(
			    "read_ena_sequences: Aspera stream ended unexpectedly (expected files for run '%s')",
			    run_.run_accession);
		}
		std::string gz_hint = filename.empty() ? run_.aspera_paths[0].remote_path : filename;
		bool is_gz = duckdb::IsGzipped(gz_hint);

		auto *s1 = CreateAsperaSeqStream(aspera_process_.get(), is_gz);
		try {
			fastx_reader_ = std::make_unique<SequenceReader>(s1, static_cast<AsperaSeqStream *>(nullptr), true);
		} catch (...) {
			delete s1;
			throw;
		}
		return;
	}

	// Paired-end: two ascp processes, one per remote path, each in stdio://
	// mode. Both pipes feed the SequenceReader concurrently — no R1 temp file.
	// Early termination (e.g., max_sequences, LIMIT-inside-LATERAL) tears both
	// down via AsperaProcess destructors, so the cap saves real bandwidth.
	if (run_.aspera_paths.size() < 2) {
		throw duckdb::IOException("read_ena_sequences: paired-end run '%s' is missing its second Aspera path",
		                          run_.run_accession);
	}
	{
		std::lock_guard<std::mutex> open_guard(open_mutex_);
		aspera_process_ =
		    std::make_unique<AsperaProcess>(config, std::vector<std::string> {run_.aspera_paths[0].remote_path});
		aspera_process_paired_ =
		    std::make_unique<AsperaProcess>(config, std::vector<std::string> {run_.aspera_paths[1].remote_path});
	}

	std::string f1;
	size_t sz1;
	std::string f2;
	size_t sz2;
	if (!aspera_process_->NextFile(f1, sz1)) {
		throw duckdb::IOException(
		    "read_ena_sequences: Aspera R1 stream ended unexpectedly (expected files for run '%s')",
		    run_.run_accession);
	}
	if (!aspera_process_paired_->NextFile(f2, sz2)) {
		throw duckdb::IOException(
		    "read_ena_sequences: Aspera R2 stream ended unexpectedly (expected files for run '%s')",
		    run_.run_accession);
	}

	bool is_gz1 = duckdb::IsGzipped(f1.empty() ? run_.aspera_paths[0].remote_path : f1);
	bool is_gz2 = duckdb::IsGzipped(f2.empty() ? run_.aspera_paths[1].remote_path : f2);

	auto *s1 = CreateAsperaSeqStream(aspera_process_.get(), is_gz1);
	AsperaSeqStream *s2 = nullptr;
	try {
		s2 = CreateAsperaSeqStream(aspera_process_paired_.get(), is_gz2);
		fastx_reader_ = std::make_unique<SequenceReader>(s1, s2, true);
	} catch (...) {
		delete s2;
		delete s1;
		throw;
	}
}

#endif // MIINT_ASPERA_SUPPORTED

} // namespace miint
