#pragma once

#include "SequenceRecord.hpp"
#include "SequenceReader.hpp"
#include "SFFReader.hpp"
#include "aspera_stream.hpp"
#include "aspera_utils.hpp"
#include "ena_run_info_extractor.hpp"
#include "duckdb/common/file_system.hpp"

#include <memory>
#include <mutex>
#include <string>

namespace miint {

// Per-run streaming reader for ENA FASTX / SFF / Aspera runs.
// Owns the reader, ascp process, and any temp files for a SINGLE run.
// Destructor is responsible for all teardown (reader destroy, ascp kill, temp unlink).
//
// Usage (typical):
//   reader = std::make_unique<PerRunReader>(fs, run, use_aspera, trim, open_mutex, aspera_cfg);
//   reader->Open();                                       // throws on failure
//   while (!reader->Exhausted()) {
//       auto batch = reader->ReadBatch(STANDARD_VECTOR_SIZE);
//       if (batch.empty()) break;
//       // consume batch
//   }
//   reader->Finish();   // for Aspera: wait for ascp + throw on bad exit
//   reader.reset();     // destructor cleans up any remaining state
//
// Retry policy stays with the caller: destroy the reader on Open/Read failure,
// create a fresh one, call Open() again. Caller also resets any per-run
// sequence counter between attempts — the counter is not owned here.
class PerRunReader {
public:
	// `max_sequences == 0` means unlimited. For paired-end runs, `max_sequences`
	// counts pairs (one batch row per pair), not underlying FASTQ records.
	PerRunReader(duckdb::FileSystem &fs, ENARunInfo run, bool use_aspera, bool trim, std::mutex &open_mutex,
	             const AsperaConfig *aspera_config, uint64_t max_sequences = 0);
	~PerRunReader();

	PerRunReader(const PerRunReader &) = delete;
	PerRunReader &operator=(const PerRunReader &) = delete;

	// Open the reader for this run. Single attempt; throws on failure.
	// Safe to call only once; subsequent calls are a no-op.
	void Open();

	// Read up to max_size records. An empty batch signals EOF and marks the
	// reader Exhausted(). Calling ReadBatch again after EOF returns empty.
	SequenceRecordBatch ReadBatch(size_t max_size);

	// Call after a successful EOF (empty batch). For Aspera runs, waits for
	// ascp to exit and throws IOException on non-zero, non-killed exit code.
	// No-op for FASTX/SFF. Idempotent.
	void Finish();

	bool Exhausted() const {
		return exhausted_;
	}

	const ENARunInfo &Run() const {
		return run_;
	}

private:
	duckdb::FileSystem &fs_;
	ENARunInfo run_;
	bool use_aspera_;
	bool trim_;
	std::mutex &open_mutex_;
	const AsperaConfig *aspera_config_; // nullable; non-null only when use_aspera_ is true

	std::unique_ptr<SequenceReader> fastx_reader_;
	std::unique_ptr<SFFReader> sff_reader_;
	std::string sff_temp_path_;
#if MIINT_ASPERA_SUPPORTED
	// Single-end and paired-end Aspera both stream directly: one ascp process
	// per remote file, each in stdio:// (single-file) mode. Paired-end uses
	// both processes concurrently; R2 is nullptr for single-end.
	std::unique_ptr<AsperaProcess> aspera_process_;
	std::unique_ptr<AsperaProcess> aspera_process_paired_;
#endif
	bool exhausted_ = false;
	bool opened_ = false;
	bool finished_ = false;
	uint64_t max_sequences_ = 0; // 0 = unlimited
	uint64_t sequences_emitted_ = 0;

	void OpenHTTP();
	void OpenSFF();
#if MIINT_ASPERA_SUPPORTED
	void OpenAspera();
#endif
};

} // namespace miint
