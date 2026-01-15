#pragma once

/*
 * Bowtie2Aligner - C++ wrapper for the bowtie2 alignment tool
 *
 * ARCHITECTURE OVERVIEW:
 * =====================
 * Unlike Minimap2Aligner which uses minimap2's C library API directly,
 * Bowtie2Aligner wraps the external bowtie2 binary because:
 * 1. Bowtie2 is GPL licensed (cannot include as submodule)
 * 2. Bowtie2 has no public C/C++ API
 *
 * STREAMING DESIGN:
 * ================
 * Bowtie2 supports Unix pipeline streaming: `cat seqs.fa | bowtie2 -x idx -f - | ...`
 * As sequences arrive on stdin, bowtie2 processes them and outputs SAM to stdout.
 * This class maintains a persistent bowtie2 process for efficient streaming:
 *
 *   1. build_index() - Creates bowtie2 index files using bowtie2-build
 *   2. First align() call - Forks bowtie2 process, sets up pipes, starts reader thread
 *   3. Subsequent align() calls - Writes queries to stdin, reads results from queue
 *   4. Destructor - Closes stdin (signals EOF), waits for bowtie2 to finish
 *
 * THREADING MODEL:
 * ===============
 * To avoid pipe buffer deadlocks (stdout fills while we're still writing to stdin),
 * we use two threads:
 *
 *   Main Thread:
 *   - Calls align() to write query batches to bowtie2 stdin
 *   - Retrieves results from thread-safe result queue
 *
 *   Reader Thread:
 *   - Continuously reads SAM records from bowtie2 stdout
 *   - Parses using SAMReader (reuses HTSlib integration)
 *   - Pushes parsed records to thread-safe result queue
 *
 * SEQUENCE FORMAT:
 * ===============
 * - If input has quality scores → FASTQ format
 * - If input lacks quality scores → FASTA format (with bowtie2 -f flag)
 * - Never fabricates fake quality scores
 *
 * USAGE PATTERN:
 * =============
 *   Bowtie2Aligner aligner(config);
 *   aligner.build_index(subjects);
 *
 *   // Can call align() multiple times - bowtie2 stays running
 *   aligner.align(batch1, output);
 *   aligner.align(batch2, output);
 *   // ... more batches ...
 *
 *   // Destructor or explicit finish() closes stdin, drains remaining output
 *
 * GLOBAL STATE WARNING:
 * ====================
 * This class sets `signal(SIGPIPE, SIG_IGN)` on first use to prevent the process
 * from being killed when writing to a pipe whose reader has closed (e.g., if
 * bowtie2 exits unexpectedly). This is a PROCESS-WIDE setting that persists
 * for the lifetime of the process and affects all code in the process, not just
 * this class.
 *
 * Implications:
 * - All writes to broken pipes will return EPIPE instead of killing the process
 * - Code that relies on SIGPIPE for notification will no longer work as expected
 * - This is set once and never restored (there's no safe way to restore it)
 *
 * This trade-off is necessary because:
 * 1. We cannot predict when bowtie2 might exit unexpectedly
 * 2. A SIGPIPE-killed process provides no opportunity for cleanup
 * 3. Handling EPIPE allows graceful error reporting
 */

#include "Minimap2Aligner.hpp" // For AlignmentSubject (shared struct)
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <filesystem>
#include <memory>
#include <mutex>
#include <queue>
#include <string>
#include <thread>
#include <vector>

namespace miint {

// ============================================================================
// Configuration
// ============================================================================

// Configuration for Bowtie2 alignment
// These map to bowtie2 command-line arguments
struct Bowtie2Config {
	std::string preset = "";      // Preset: --very-fast, --fast, --sensitive, --very-sensitive
	bool local = false;           // --local vs --end-to-end (default)
	int threads = 1;              // -p parameter (bowtie2 internal threading)
	int max_secondary = 0;        // -k parameter (0 = default bowtie2 behavior)
	std::string extra_args = "";  // Additional bowtie2 arguments (space-separated)
	bool quiet = true;            // Suppress stderr output (alignment statistics)
};

// ============================================================================
// Aligner Class
// ============================================================================

class Bowtie2Aligner {
public:
	explicit Bowtie2Aligner(const Bowtie2Config &config);
	~Bowtie2Aligner();

	// Non-copyable (owns temp directory, process handles, threads)
	Bowtie2Aligner(const Bowtie2Aligner &) = delete;
	Bowtie2Aligner &operator=(const Bowtie2Aligner &) = delete;

	// Movable (transfers ownership)
	// Move operations throw if aligner has an active process (thread safety)
	Bowtie2Aligner(Bowtie2Aligner &&other);
	Bowtie2Aligner &operator=(Bowtie2Aligner &&other);

	// ========================================================================
	// Index Building / Loading
	// ========================================================================

	// Build index from multiple reference sequences
	// Creates index files in temp directory using bowtie2-build
	// Must be called before align()
	void build_index(const std::vector<AlignmentSubject> &subjects);

	// Convenience: build index from a single reference
	void build_single_index(const AlignmentSubject &subject);

	// Load a pre-built index from prefix (validates files exist)
	// The prefix should be the path without extension, e.g., "/path/to/index"
	// which expects files like /path/to/index.1.bt2, /path/to/index.rev.1.bt2, etc.
	void load_index(const std::string &index_prefix);

	// Check if a bowtie2 index exists at the given prefix
	// Returns true if the required index files (.1.bt2, .2.bt2, .rev.1.bt2, .rev.2.bt2) exist
	static bool is_index_prefix(const std::string &prefix);

	// ========================================================================
	// Alignment
	// ========================================================================

	// Align query sequences against the built index
	// - First call starts the bowtie2 process and reader thread
	// - Writes queries to bowtie2 stdin (FASTA if no quals, FASTQ if quals present)
	// - Drains any currently available results to output (non-blocking)
	// - Can be called multiple times; bowtie2 process persists between calls
	// - Call finish() after all batches to get remaining results
	void align(const SequenceRecordBatch &queries, SAMRecordBatch &output);

	// Finish alignment: close stdin, drain remaining output, wait for process
	// Must be called after all align() calls to collect remaining results
	// After finish(), call reset() to reuse for another index, or destroy the instance
	void finish(SAMRecordBatch &output);

	// Reset aligner state for reuse with a new index (after finish())
	// This is more efficient than creating a new Bowtie2Aligner instance because:
	// - Preserves validated executable paths (avoids fork/exec for 'which' commands)
	// - Preserves temp directory (if one was created for build_index())
	// Usage: after finish(), call reset(), then load_index() or build_index(), then align()
	void reset();

	// ========================================================================
	// Accessors (for testing)
	// ========================================================================

	std::filesystem::path get_temp_dir() const {
		return temp_dir_;
	}

	bool is_aligner_running() const {
		return aligner_running_;
	}

private:
	// ========================================================================
	// Configuration and State
	// ========================================================================

	Bowtie2Config config_;
	std::string bowtie2_path_;       // Path to bowtie2 binary
	std::string bowtie2_build_path_; // Path to bowtie2-build binary
	std::filesystem::path temp_dir_; // Temp directory for index and intermediate files
	std::string index_prefix_;       // Base path for index files (temp_dir_/index)
	bool index_built_ = false;       // True after build_index() succeeds

	// ========================================================================
	// Persistent Process State
	// ========================================================================

	pid_t bowtie2_pid_ = -1;         // PID of running bowtie2 process
	FILE *stdin_pipe_ = nullptr;     // Pipe to write queries to bowtie2 stdin
	// Note: stdout pipe fd is transferred to reader thread (SAMReader) immediately after fork
	bool aligner_running_ = false;   // True while bowtie2 process is active
	bool use_fasta_ = false;         // True if first batch had no quality scores
	bool is_paired_ = false;         // True if first batch was paired-end

	// ========================================================================
	// Reader Thread and Result Queue
	// ========================================================================

	std::thread reader_thread_;                  // Background thread reading SAM output
	std::atomic<bool> reader_should_stop_{false}; // Signal reader thread to stop

	// Thread-safe queue for passing SAM records from reader to main thread
	std::mutex queue_mutex_;
	std::condition_variable queue_cv_;
	std::queue<SAMRecordBatch> result_queue_;
	bool reader_finished_ = false;               // True when reader thread exits
	std::string reader_error_;                   // Error message from reader thread

	// ========================================================================
	// Binary Discovery
	// ========================================================================

	// Find executable in PATH using 'which' command
	std::string find_executable(const std::string &name);

	// Validate that bowtie2 and bowtie2-build are in PATH
	void validate_executables();

	// ========================================================================
	// Temp Directory Management
	// ========================================================================

	// Create unique temp directory for this instance
	void create_temp_directory();

	// Remove temp directory and all contents
	void cleanup_temp_directory();

	// ========================================================================
	// Index Building Helpers
	// ========================================================================

	// Write reference sequences to FASTA file
	void write_subjects_fasta(const std::vector<AlignmentSubject> &subjects,
	                          const std::filesystem::path &fasta_path);

	// Run bowtie2-build to create index files
	void run_bowtie2_build(const std::filesystem::path &fasta_path);

	// ========================================================================
	// Alignment Process Management
	// ========================================================================

	// Start the bowtie2 process with pipes for stdin/stdout
	// Called on first align() invocation
	void start_aligner_process(bool use_fasta, bool is_paired);

	// Stop the bowtie2 process and reader thread
	// Closes stdin, waits for process to finish, joins reader thread
	void stop_aligner_process();

	// Reader thread function: reads SAM from bowtie2 stdout, pushes to queue
	// Takes ownership of the fd parameter
	void reader_thread_func(int fd);

	// ========================================================================
	// Query Writing Helpers
	// ========================================================================

	// Write queries in FASTQ format (when quality scores present)
	void write_queries_fastq(FILE *fp, const SequenceRecordBatch &queries);

	// Write queries in FASTA format (when no quality scores)
	void write_queries_fasta(FILE *fp, const SequenceRecordBatch &queries);

	// Write paired-end queries in interleaved FASTQ format
	void write_queries_interleaved_fastq(FILE *fp, const SequenceRecordBatch &queries);

	// Write paired-end queries in interleaved FASTA format
	void write_queries_interleaved_fasta(FILE *fp, const SequenceRecordBatch &queries);

	// Check if any query in the batch has quality scores
	static bool has_quality_scores(const SequenceRecordBatch &queries);

	// ========================================================================
	// Result Collection
	// ========================================================================

	// Drain currently available results from queue into output (non-blocking)
	void drain_results(SAMRecordBatch &output);

	// Cleanup aligner process without draining results (for destructor)
	void cleanup_process();
};

} // namespace miint
