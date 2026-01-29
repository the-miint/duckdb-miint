#include "Bowtie2Aligner.hpp"
#include "SAMReader.hpp"
#include <array>
#include <chrono>
#include <csignal>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <sstream>
#include <stdexcept>
#include <fcntl.h>
#include <sys/wait.h>
#include <unistd.h>

#ifdef __linux__
#include <sys/prctl.h>
#endif

namespace miint {

// =============================================================================
// Binary Discovery
// =============================================================================

std::string Bowtie2Aligner::find_executable(const std::string &name) {
	// Use fork/exec to run 'which' - avoids shell injection risks
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		return "";
	}

	pid_t pid = fork();
	if (pid == -1) {
		close(pipefd[0]);
		close(pipefd[1]);
		return "";
	}

	if (pid == 0) {
		// Child process
		close(pipefd[0]); // Close read end
		dup2(pipefd[1], STDOUT_FILENO);
		close(pipefd[1]);

		// Redirect stderr to /dev/null
		int devnull = open("/dev/null", O_WRONLY);
		if (devnull != -1) {
			dup2(devnull, STDERR_FILENO);
			close(devnull);
		}

		execlp("which", "which", name.c_str(), nullptr);
		_exit(127); // exec failed
	}

	// Parent process
	close(pipefd[1]); // Close write end

	std::array<char, 256> buffer;
	std::string result;
	ssize_t n;
	while ((n = read(pipefd[0], buffer.data(), buffer.size())) > 0) {
		result.append(buffer.data(), n);
	}
	close(pipefd[0]);

	int status;
	waitpid(pid, &status, 0);

	if (!WIFEXITED(status) || WEXITSTATUS(status) != 0) {
		return "";
	}

	// Remove trailing newline
	while (!result.empty() && (result.back() == '\n' || result.back() == '\r')) {
		result.pop_back();
	}

	return result;
}

void Bowtie2Aligner::validate_executables() {
	bowtie2_path_ = find_executable("bowtie2");
	bowtie2_build_path_ = find_executable("bowtie2-build");

	if (bowtie2_path_.empty() || bowtie2_build_path_.empty()) {
		std::string missing;
		if (bowtie2_path_.empty()) {
			missing = "bowtie2";
		}
		if (bowtie2_build_path_.empty()) {
			if (!missing.empty()) {
				missing += " and ";
			}
			missing += "bowtie2-build";
		}
		throw std::runtime_error(missing + " not found in PATH. "
		                                   "Please install bowtie2 and ensure it is in your PATH.");
	}
}

// =============================================================================
// Temp Directory Management
// =============================================================================

void Bowtie2Aligner::create_temp_directory() {
	// Generate unique directory name: bowtie2_<pid>_<timestamp>_<random>
	auto now = std::chrono::steady_clock::now();
	auto timestamp = std::chrono::duration_cast<std::chrono::nanoseconds>(now.time_since_epoch()).count();

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_int_distribution<> dis(0, 999999);

	std::ostringstream oss;
	oss << "bowtie2_" << getpid() << "_" << timestamp << "_" << dis(gen);

	temp_dir_ = std::filesystem::temp_directory_path() / oss.str();

	// Create the directory
	std::error_code ec;
	if (!std::filesystem::create_directories(temp_dir_, ec)) {
		if (ec) {
			throw std::runtime_error("Failed to create temp directory: " + ec.message());
		}
	}
}

void Bowtie2Aligner::cleanup_temp_directory() {
	if (!temp_dir_.empty() && std::filesystem::exists(temp_dir_)) {
		std::error_code ec;
		std::filesystem::remove_all(temp_dir_, ec);
		// Ignore errors during cleanup (best effort)
	}
	temp_dir_.clear();
}

// =============================================================================
// Constructor / Destructor
// =============================================================================

Bowtie2Aligner::Bowtie2Aligner(const Bowtie2Config &config) : config_(config) {
	// Validate that bowtie2 binaries are available
	validate_executables();

	// Create unique temp directory for this instance
	create_temp_directory();

	// Set index prefix (will be used when building index)
	index_prefix_ = (temp_dir_ / "index").string();
}

Bowtie2Aligner::~Bowtie2Aligner() {
	// Stop any running alignment process before cleanup
	cleanup_process();
	cleanup_temp_directory();
}

// =============================================================================
// Move Semantics
// =============================================================================

Bowtie2Aligner::Bowtie2Aligner(Bowtie2Aligner &&other)
    : config_(std::move(other.config_)), bowtie2_path_(std::move(other.bowtie2_path_)),
      bowtie2_build_path_(std::move(other.bowtie2_build_path_)), temp_dir_(std::move(other.temp_dir_)),
      index_prefix_(std::move(other.index_prefix_)), index_built_(other.index_built_) {
	// Cannot safely move an aligner with an active process - the reader thread
	// would still be accessing the source object's mutex/queue
	if (other.aligner_running_) {
		throw std::runtime_error("Cannot move Bowtie2Aligner with active alignment process. "
		                         "Call finish() before moving.");
	}

	// Safe to move - no active process
	bowtie2_pid_ = other.bowtie2_pid_;
	stdin_pipe_ = other.stdin_pipe_;
	aligner_running_ = other.aligner_running_;
	use_fasta_ = other.use_fasta_;
	is_paired_ = other.is_paired_;
	reader_finished_ = other.reader_finished_;
	reader_error_ = std::move(other.reader_error_);

	// Transfer ownership - clear other's state so it won't cleanup on destruction
	other.temp_dir_.clear();
	other.index_built_ = false;
	other.bowtie2_pid_ = -1;
	other.stdin_pipe_ = nullptr;
	other.aligner_running_ = false;
	other.use_fasta_ = false;
	other.is_paired_ = false;
	other.reader_finished_ = false;
}

Bowtie2Aligner &Bowtie2Aligner::operator=(Bowtie2Aligner &&other) {
	if (this != &other) {
		// Cannot safely move an aligner with an active process
		if (other.aligner_running_) {
			throw std::runtime_error("Cannot move Bowtie2Aligner with active alignment process. "
			                         "Call finish() before moving.");
		}

		// Stop our current alignment process and cleanup
		cleanup_process();
		cleanup_temp_directory();

		// Move from other
		config_ = std::move(other.config_);
		bowtie2_path_ = std::move(other.bowtie2_path_);
		bowtie2_build_path_ = std::move(other.bowtie2_build_path_);
		temp_dir_ = std::move(other.temp_dir_);
		index_prefix_ = std::move(other.index_prefix_);
		index_built_ = other.index_built_;
		bowtie2_pid_ = other.bowtie2_pid_;
		stdin_pipe_ = other.stdin_pipe_;
		aligner_running_ = other.aligner_running_;
		use_fasta_ = other.use_fasta_;
		is_paired_ = other.is_paired_;
		reader_finished_ = other.reader_finished_;
		reader_error_ = std::move(other.reader_error_);

		// Clear other's state so it won't cleanup on destruction
		other.temp_dir_.clear();
		other.index_built_ = false;
		other.bowtie2_pid_ = -1;
		other.stdin_pipe_ = nullptr;
		other.aligner_running_ = false;
		other.use_fasta_ = false;
		other.is_paired_ = false;
		other.reader_finished_ = false;
	}
	return *this;
}

// =============================================================================
// Index Building
// =============================================================================

void Bowtie2Aligner::write_subjects_fasta(const std::vector<AlignmentSubject> &subjects,
                                          const std::filesystem::path &fasta_path) {
	FILE *fp = fopen(fasta_path.c_str(), "w");
	if (!fp) {
		throw std::runtime_error("Failed to create FASTA file: " + fasta_path.string());
	}

	for (const auto &subject : subjects) {
		// Write FASTA format: >name\nsequence\n
		fprintf(fp, ">%s\n%s\n", subject.read_id.c_str(), subject.sequence.c_str());
	}

	fclose(fp);
}

void Bowtie2Aligner::run_bowtie2_build(const std::filesystem::path &fasta_path) {
	// Use fork/exec to run bowtie2-build - avoids shell injection risks
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		throw std::runtime_error("Failed to create pipe for bowtie2-build: " + std::string(strerror(errno)));
	}

	pid_t pid = fork();
	if (pid == -1) {
		close(pipefd[0]);
		close(pipefd[1]);
		throw std::runtime_error("Failed to fork for bowtie2-build: " + std::string(strerror(errno)));
	}

	if (pid == 0) {
		// Child process
		close(pipefd[0]); // Close read end

		// Redirect both stdout and stderr to the pipe for error capture
		dup2(pipefd[1], STDOUT_FILENO);
		dup2(pipefd[1], STDERR_FILENO);
		close(pipefd[1]);

		// Execute bowtie2-build with arguments passed directly (no shell)
		std::string fasta_str = fasta_path.string();
		execl(bowtie2_build_path_.c_str(), bowtie2_build_path_.c_str(), "--quiet", fasta_str.c_str(),
		      index_prefix_.c_str(), nullptr);

		// If exec fails, write error and exit
		std::string err = "Failed to execute bowtie2-build: " + std::string(strerror(errno));
		write(STDERR_FILENO, err.c_str(), err.size());
		_exit(127);
	}

	// Parent process
	close(pipefd[1]); // Close write end

	// Capture output for error reporting
	std::array<char, 256> buffer;
	std::string output;
	ssize_t n;
	while ((n = read(pipefd[0], buffer.data(), buffer.size())) > 0) {
		output.append(buffer.data(), n);
	}
	close(pipefd[0]);

	int status;
	waitpid(pid, &status, 0);

	if (!WIFEXITED(status)) {
		throw std::runtime_error("bowtie2-build terminated abnormally: " + output);
	}

	int exit_code = WEXITSTATUS(status);
	if (exit_code != 0) {
		throw std::runtime_error("bowtie2-build failed with exit code " + std::to_string(exit_code) + ": " + output);
	}
}

void Bowtie2Aligner::build_index(const std::vector<AlignmentSubject> &subjects) {
	if (subjects.empty()) {
		throw std::runtime_error("Cannot build index from empty subject list");
	}

	// Write subjects to FASTA file
	auto fasta_path = temp_dir_ / "reference.fasta";
	write_subjects_fasta(subjects, fasta_path);

	// Run bowtie2-build
	run_bowtie2_build(fasta_path);

	index_built_ = true;
}

void Bowtie2Aligner::build_single_index(const AlignmentSubject &subject) {
	std::vector<AlignmentSubject> subjects = {subject};
	build_index(subjects);
}

bool Bowtie2Aligner::is_index_prefix(const std::string &prefix) {
	// Bowtie2 creates 6 index files for small genomes (<4GB):
	// .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2
	// For large genomes, it creates .1.bt2l files instead
	// We check for the minimum required files
	std::vector<std::string> required_extensions = {".1.bt2", ".2.bt2", ".rev.1.bt2", ".rev.2.bt2"};

	// Also check for large index format (.bt2l)
	std::vector<std::string> required_extensions_large = {".1.bt2l", ".2.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"};

	// Check standard format first
	bool all_exist = true;
	for (const auto &ext : required_extensions) {
		if (!std::filesystem::exists(prefix + ext)) {
			all_exist = false;
			break;
		}
	}
	if (all_exist) {
		return true;
	}

	// Check large index format
	all_exist = true;
	for (const auto &ext : required_extensions_large) {
		if (!std::filesystem::exists(prefix + ext)) {
			all_exist = false;
			break;
		}
	}
	return all_exist;
}

void Bowtie2Aligner::load_index(const std::string &index_prefix) {
	if (!is_index_prefix(index_prefix)) {
		throw std::runtime_error("No valid bowtie2 index found at prefix: " + index_prefix + ". Expected files like " +
		                         index_prefix + ".1.bt2, " + index_prefix + ".rev.1.bt2, etc.");
	}

	index_prefix_ = index_prefix;
	index_built_ = true;
}

// =============================================================================
// Alignment - Query Writing
// =============================================================================

bool Bowtie2Aligner::has_quality_scores(const SequenceRecordBatch &queries) {
	// Check if any query has non-empty quality scores
	for (const auto &qual : queries.quals1) {
		if (!qual.as_string().empty()) {
			return true;
		}
	}
	return false;
}

void Bowtie2Aligner::write_queries_fasta(FILE *fp, const SequenceRecordBatch &queries) {
	for (size_t i = 0; i < queries.read_ids.size(); ++i) {
		// FASTA format: >read_id\nsequence\n
		// Fix #8: Check fprintf return value
		if (fprintf(fp, ">%s\n%s\n", queries.read_ids[i].c_str(), queries.sequences1[i].c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
	}
	if (fflush(fp) != 0) {
		throw std::runtime_error("Failed to flush bowtie2 stdin: " + std::string(strerror(errno)));
	}
}

void Bowtie2Aligner::write_queries_fastq(FILE *fp, const SequenceRecordBatch &queries) {
	for (size_t i = 0; i < queries.read_ids.size(); ++i) {
		// FASTQ format: @read_id\nsequence\n+\nquality\n
		// Fix #8: Check fprintf return value
		if (fprintf(fp, "@%s\n%s\n+\n%s\n", queries.read_ids[i].c_str(), queries.sequences1[i].c_str(),
		            queries.quals1[i].as_string().c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
	}
	if (fflush(fp) != 0) {
		throw std::runtime_error("Failed to flush bowtie2 stdin: " + std::string(strerror(errno)));
	}
}

void Bowtie2Aligner::write_queries_interleaved_fasta(FILE *fp, const SequenceRecordBatch &queries) {
	// Interleaved format: R1, R2, R1, R2, ...
	// bowtie2 with --interleaved expects pairs of consecutive reads
	for (size_t i = 0; i < queries.read_ids.size(); ++i) {
		// Write R1 (read 1 of pair)
		if (fprintf(fp, ">%s/1\n%s\n", queries.read_ids[i].c_str(), queries.sequences1[i].c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
		// Write R2 (read 2 of pair)
		if (fprintf(fp, ">%s/2\n%s\n", queries.read_ids[i].c_str(), queries.sequences2[i].c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
	}
	if (fflush(fp) != 0) {
		throw std::runtime_error("Failed to flush bowtie2 stdin: " + std::string(strerror(errno)));
	}
}

void Bowtie2Aligner::write_queries_interleaved_fastq(FILE *fp, const SequenceRecordBatch &queries) {
	// Interleaved FASTQ format: R1, R2, R1, R2, ...
	for (size_t i = 0; i < queries.read_ids.size(); ++i) {
		// Write R1 (read 1 of pair)
		if (fprintf(fp, "@%s/1\n%s\n+\n%s\n", queries.read_ids[i].c_str(), queries.sequences1[i].c_str(),
		            queries.quals1[i].as_string().c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
		// Write R2 (read 2 of pair)
		if (fprintf(fp, "@%s/2\n%s\n+\n%s\n", queries.read_ids[i].c_str(), queries.sequences2[i].c_str(),
		            queries.quals2[i].as_string().c_str()) < 0) {
			throw std::runtime_error("Failed to write query to bowtie2 stdin: " + std::string(strerror(errno)));
		}
	}
	if (fflush(fp) != 0) {
		throw std::runtime_error("Failed to flush bowtie2 stdin: " + std::string(strerror(errno)));
	}
}

// =============================================================================
// Alignment - Process Management
// =============================================================================

void Bowtie2Aligner::start_aligner_process(bool use_fasta, bool is_paired) {
	if (aligner_running_) {
		return; // Already running
	}

	use_fasta_ = use_fasta;
	is_paired_ = is_paired;

	// Create pipes for stdin and stdout
	int stdin_pipe[2];  // [0] = read end, [1] = write end
	int stdout_pipe[2]; // [0] = read end, [1] = write end

	if (pipe(stdin_pipe) == -1) {
		throw std::runtime_error("Failed to create stdin pipe: " + std::string(strerror(errno)));
	}
	if (pipe(stdout_pipe) == -1) {
		close(stdin_pipe[0]);
		close(stdin_pipe[1]);
		throw std::runtime_error("Failed to create stdout pipe: " + std::string(strerror(errno)));
	}

	// Fork the process
	bowtie2_pid_ = fork();
	if (bowtie2_pid_ == -1) {
		close(stdin_pipe[0]);
		close(stdin_pipe[1]);
		close(stdout_pipe[0]);
		close(stdout_pipe[1]);
		throw std::runtime_error("Failed to fork bowtie2 process: " + std::string(strerror(errno)));
	}

	if (bowtie2_pid_ == 0) {
		// Child process - becomes bowtie2

		// If parent dies, kill this child too (Linux-specific)
#ifdef __linux__
		prctl(PR_SET_PDEATHSIG, SIGTERM);
#endif

		// Redirect stdin to read end of stdin_pipe
		dup2(stdin_pipe[0], STDIN_FILENO);
		close(stdin_pipe[0]);
		close(stdin_pipe[1]); // Close write end in child

		// Redirect stdout to write end of stdout_pipe
		dup2(stdout_pipe[1], STDOUT_FILENO);
		close(stdout_pipe[0]); // Close read end in child
		close(stdout_pipe[1]);

		// Redirect stderr to /dev/null if quiet mode
		if (config_.quiet) {
			int devnull = open("/dev/null", O_WRONLY);
			if (devnull != -1) {
				dup2(devnull, STDERR_FILENO);
				close(devnull);
			}
		}

		// Build argument list
		// Reserve enough space upfront to avoid reallocation (which would invalidate c_str() pointers)
		std::vector<std::string> args_storage;
		args_storage.reserve(20); // More than enough for all arguments
		std::vector<char *> args;

		args_storage.push_back(bowtie2_path_);
		args.push_back(const_cast<char *>(args_storage.back().c_str()));

		// Index prefix
		args_storage.push_back("-x");
		args.push_back(const_cast<char *>(args_storage.back().c_str()));
		args_storage.push_back(index_prefix_);
		args.push_back(const_cast<char *>(args_storage.back().c_str()));

		// Input format and source
		if (use_fasta) {
			args_storage.push_back("-f"); // FASTA input
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}

		// Paired-end mode uses --interleaved, single-end uses just - for stdin
		if (is_paired) {
			// --interleaved takes the filename as its argument
			args_storage.push_back("--interleaved");
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}
		args_storage.push_back("-"); // Read from stdin
		args.push_back(const_cast<char *>(args_storage.back().c_str()));

		// Preset if specified (bowtie2 presets need -- prefix, e.g. --very-sensitive)
		if (!config_.preset.empty()) {
			args_storage.push_back("--" + config_.preset);
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}

		// Local vs end-to-end
		if (config_.local) {
			args_storage.push_back("--local");
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}

		// Threads
		if (config_.threads > 1) {
			args_storage.push_back("-p");
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
			args_storage.push_back(std::to_string(config_.threads));
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}

		// Max secondary alignments (-k)
		if (config_.max_secondary > 0) {
			args_storage.push_back("-k");
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
			args_storage.push_back(std::to_string(config_.max_secondary));
			args.push_back(const_cast<char *>(args_storage.back().c_str()));
		}

		// Extra arguments (space-separated)
		if (!config_.extra_args.empty()) {
			std::istringstream iss(config_.extra_args);
			std::string arg;
			while (iss >> arg) {
				args_storage.push_back(arg);
				args.push_back(const_cast<char *>(args_storage.back().c_str()));
			}
		}

		args.push_back(nullptr); // Null-terminate

		// Execute bowtie2
		execv(bowtie2_path_.c_str(), args.data());

		// If execv returns, it failed
		_exit(127);
	}

	// Parent process
	// Note: SIGPIPE is set to SIG_IGN at extension initialization (miint_extension.cpp)
	// to ensure thread-safe signal handling across all subprocess management.

	// Close unused pipe ends
	close(stdin_pipe[0]);  // Close read end of stdin pipe
	close(stdout_pipe[1]); // Close write end of stdout pipe

	// Wrap write end of stdin pipe in FILE* for fprintf
	stdin_pipe_ = fdopen(stdin_pipe[1], "w");
	if (!stdin_pipe_) {
		// Fix #1: Kill child process and clean up on fdopen failure
		kill(bowtie2_pid_, SIGTERM);
		waitpid(bowtie2_pid_, nullptr, 0);
		bowtie2_pid_ = -1;
		close(stdin_pipe[1]);
		close(stdout_pipe[0]);
		throw std::runtime_error("Failed to wrap stdin pipe in FILE*");
	}

	// Set line buffering on stdin pipe to ensure data is sent to bowtie2 immediately
	// Without this, data may be stuck in stdio buffer even after fflush() in some edge cases
	setvbuf(stdin_pipe_, nullptr, _IOLBF, 0);

	// Fix #3/#4: Transfer fd ownership to reader thread immediately
	// Pass fd by value to thread, so ownership is clear from the start
	int reader_fd = stdout_pipe[0];
	// stdout pipe fd transferred to reader thread (SAMReader owns it now)

	aligner_running_ = true;
	reader_finished_ = false;
	reader_error_.clear();
	reader_should_stop_.store(false);

	// Start the reader thread with the fd
	reader_thread_ = std::thread(&Bowtie2Aligner::reader_thread_func, this, reader_fd);
}

void Bowtie2Aligner::reader_thread_func(int fd) {
	// Fix #3/#4: fd is passed by value, so ownership is unambiguous
	// This thread now owns the fd and is responsible for closing it (via SAMReader)
	try {
		// Create SAMReader from the stdout pipe file descriptor
		// SAMReader takes ownership of the fd and will close it on destruction
		SAMReader reader(fd, "bowtie2_stdout", false);

		// Read batches of SAM records and push to queue
		const int batch_size = 1000;
		while (!reader_should_stop_.load()) {
			SAMRecordBatch batch = reader.read(batch_size);
			if (batch.empty()) {
				// EOF reached
				break;
			}

			// Push batch to queue
			{
				std::lock_guard<std::mutex> lock(queue_mutex_);
				result_queue_.push(std::move(batch));
			}
			queue_cv_.notify_one();
		}
	} catch (const std::exception &e) {
		// Catch any exception from SAMReader constructor or read() to prevent std::terminate
		// fd is closed by SAMReader destructor (hdopen transfers ownership to hFILE)
		std::lock_guard<std::mutex> lock(queue_mutex_);
		reader_error_ = e.what();
	}

	// Mark reader as finished
	{
		std::lock_guard<std::mutex> lock(queue_mutex_);
		reader_finished_ = true;
	}
	queue_cv_.notify_all();
}

void Bowtie2Aligner::stop_aligner_process() {
	if (!aligner_running_) {
		return;
	}

	// Signal reader thread to stop
	reader_should_stop_.store(true);

	// Close stdin pipe - this signals EOF to bowtie2
	if (stdin_pipe_) {
		fclose(stdin_pipe_);
		stdin_pipe_ = nullptr;
	}

	// Join reader thread (it will exit when bowtie2 closes stdout)
	if (reader_thread_.joinable()) {
		reader_thread_.join();
	}

	// Note: stdout pipe is owned by reader thread (passed to SAMReader), don't close here

	// Wait for bowtie2 process to exit
	if (bowtie2_pid_ > 0) {
		int status;
		waitpid(bowtie2_pid_, &status, 0);
		bowtie2_pid_ = -1;
	}

	aligner_running_ = false;
}

void Bowtie2Aligner::cleanup_process() {
	if (!aligner_running_) {
		return;
	}

	// Close stdin to signal EOF to bowtie2
	// This causes bowtie2 to finish processing and close stdout
	if (stdin_pipe_) {
		fclose(stdin_pipe_);
		stdin_pipe_ = nullptr;
	}

	// Fix #5: Wait for reader thread with timeout to avoid deadlock
	// If bowtie2 hangs, we don't want to block forever
	if (reader_thread_.joinable()) {
		// Give bowtie2 reasonable time to finish (30 seconds)
		auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(30);
		bool joined = false;

		while (std::chrono::steady_clock::now() < deadline) {
			// Check if thread finished
			{
				std::lock_guard<std::mutex> lock(queue_mutex_);
				if (reader_finished_) {
					joined = true;
					break;
				}
			}
			std::this_thread::sleep_for(std::chrono::milliseconds(100));
		}

		if (!joined) {
			// Thread is hung - kill bowtie2 to unblock it
			if (bowtie2_pid_ > 0) {
				kill(bowtie2_pid_, SIGKILL);
				// Reap immediately to avoid zombie process
				waitpid(bowtie2_pid_, nullptr, 0);
				bowtie2_pid_ = -1;
			}
		}

		reader_thread_.join();
	}

	// Now safe to signal stop (though thread should already be done)
	reader_should_stop_.store(true);

	// Note: stdout pipe is owned by reader thread (passed to SAMReader), don't close here

	// Wait for bowtie2 process and check exit status
	if (bowtie2_pid_ > 0) {
		int status;
		if (waitpid(bowtie2_pid_, &status, 0) > 0) {
			if (WIFEXITED(status)) {
				int exit_code = WEXITSTATUS(status);
				if (exit_code != 0) {
					// Store error for later retrieval
					std::lock_guard<std::mutex> lock(queue_mutex_);
					if (reader_error_.empty()) {
						reader_error_ = "bowtie2 exited with code " + std::to_string(exit_code);
					}
				}
			} else if (WIFSIGNALED(status)) {
				int sig = WTERMSIG(status);
				std::lock_guard<std::mutex> lock(queue_mutex_);
				if (reader_error_.empty()) {
					reader_error_ = "bowtie2 killed by signal " + std::to_string(sig);
				}
			}
		}
		bowtie2_pid_ = -1;
	}

	aligner_running_ = false;
}

// =============================================================================
// Alignment - Result Draining Helper
// =============================================================================

void Bowtie2Aligner::drain_results(SAMRecordBatch &output) {
	// Drain all currently available results from queue (non-blocking)
	std::lock_guard<std::mutex> lock(queue_mutex_);

	// Check for reader errors
	if (!reader_error_.empty()) {
		throw std::runtime_error("Reader thread error: " + reader_error_);
	}

	while (!result_queue_.empty()) {
		SAMRecordBatch &batch = result_queue_.front();
		// Append batch to output - core fields
		output.read_ids.insert(output.read_ids.end(), std::make_move_iterator(batch.read_ids.begin()),
		                       std::make_move_iterator(batch.read_ids.end()));
		output.flags.insert(output.flags.end(), batch.flags.begin(), batch.flags.end());
		output.references.insert(output.references.end(), std::make_move_iterator(batch.references.begin()),
		                         std::make_move_iterator(batch.references.end()));
		output.positions.insert(output.positions.end(), batch.positions.begin(), batch.positions.end());
		output.stop_positions.insert(output.stop_positions.end(), batch.stop_positions.begin(),
		                             batch.stop_positions.end());
		output.mapqs.insert(output.mapqs.end(), batch.mapqs.begin(), batch.mapqs.end());
		output.cigars.insert(output.cigars.end(), std::make_move_iterator(batch.cigars.begin()),
		                     std::make_move_iterator(batch.cigars.end()));
		output.mate_references.insert(output.mate_references.end(),
		                              std::make_move_iterator(batch.mate_references.begin()),
		                              std::make_move_iterator(batch.mate_references.end()));
		output.mate_positions.insert(output.mate_positions.end(), batch.mate_positions.begin(),
		                             batch.mate_positions.end());
		output.template_lengths.insert(output.template_lengths.end(), batch.template_lengths.begin(),
		                               batch.template_lengths.end());
		// Tag fields
		output.tag_as_values.insert(output.tag_as_values.end(), batch.tag_as_values.begin(), batch.tag_as_values.end());
		output.tag_xs_values.insert(output.tag_xs_values.end(), batch.tag_xs_values.begin(), batch.tag_xs_values.end());
		output.tag_ys_values.insert(output.tag_ys_values.end(), batch.tag_ys_values.begin(), batch.tag_ys_values.end());
		output.tag_xn_values.insert(output.tag_xn_values.end(), batch.tag_xn_values.begin(), batch.tag_xn_values.end());
		output.tag_xm_values.insert(output.tag_xm_values.end(), batch.tag_xm_values.begin(), batch.tag_xm_values.end());
		output.tag_xo_values.insert(output.tag_xo_values.end(), batch.tag_xo_values.begin(), batch.tag_xo_values.end());
		output.tag_xg_values.insert(output.tag_xg_values.end(), batch.tag_xg_values.begin(), batch.tag_xg_values.end());
		output.tag_nm_values.insert(output.tag_nm_values.end(), batch.tag_nm_values.begin(), batch.tag_nm_values.end());
		output.tag_yt_values.insert(output.tag_yt_values.end(), std::make_move_iterator(batch.tag_yt_values.begin()),
		                            std::make_move_iterator(batch.tag_yt_values.end()));
		output.tag_md_values.insert(output.tag_md_values.end(), std::make_move_iterator(batch.tag_md_values.begin()),
		                            std::make_move_iterator(batch.tag_md_values.end()));
		output.tag_sa_values.insert(output.tag_sa_values.end(), std::make_move_iterator(batch.tag_sa_values.begin()),
		                            std::make_move_iterator(batch.tag_sa_values.end()));
		// Optional seq/qual fields
		output.sequences.insert(output.sequences.end(), std::make_move_iterator(batch.sequences.begin()),
		                        std::make_move_iterator(batch.sequences.end()));
		output.quals.insert(output.quals.end(), std::make_move_iterator(batch.quals.begin()),
		                    std::make_move_iterator(batch.quals.end()));
		result_queue_.pop();
	}
}

// =============================================================================
// Alignment - Main Interface
// =============================================================================

void Bowtie2Aligner::align(const SequenceRecordBatch &queries, SAMRecordBatch &output) {
	if (queries.empty()) {
		return;
	}

	if (!index_built_) {
		throw std::runtime_error("No index built. Call build_index() first.");
	}

	// Start aligner process on first call
	if (!aligner_running_) {
		bool use_fasta = !has_quality_scores(queries);
		bool is_paired = queries.is_paired;
		start_aligner_process(use_fasta, is_paired);
	} else {
		// Check for paired/unpaired mismatch - cannot mix after process started
		if (queries.is_paired != is_paired_) {
			throw std::runtime_error("Cannot mix paired and unpaired batches in the same alignment session");
		}
	}

	// Write queries to bowtie2 stdin
	if (is_paired_) {
		// Paired-end: write interleaved format
		if (use_fasta_) {
			write_queries_interleaved_fasta(stdin_pipe_, queries);
		} else {
			write_queries_interleaved_fastq(stdin_pipe_, queries);
		}
	} else {
		// Single-end
		if (use_fasta_) {
			write_queries_fasta(stdin_pipe_, queries);
		} else {
			write_queries_fastq(stdin_pipe_, queries);
		}
	}

	// Drain any results that are currently available (non-blocking).
	// For large datasets, bowtie2 streams output continuously.
	// For small test batches, results may not be available until finish().
	drain_results(output);
}

void Bowtie2Aligner::finish(SAMRecordBatch &output) {
	// Stop the process and wait for completion
	cleanup_process();

	// Drain remaining results
	drain_results(output);
}

void Bowtie2Aligner::reset() {
	// Ensure any running process is stopped
	if (aligner_running_) {
		cleanup_process();
	}

	// Reset process state (preserves executable paths and temp directory)
	index_built_ = false;
	index_prefix_.clear();
	aligner_running_ = false;
	use_fasta_ = false;
	is_paired_ = false;
	reader_finished_ = false;
	reader_error_.clear();
	reader_should_stop_.store(false);

	// Clear result queue
	{
		std::lock_guard<std::mutex> lock(queue_mutex_);
		while (!result_queue_.empty()) {
			result_queue_.pop();
		}
	}
}

} // namespace miint
