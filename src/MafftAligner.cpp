#include "include/MafftAligner.hpp"
#include <cerrno>
#include <climits>
#include <cstdio>
#include <cstring>
#include <stdexcept>
#include <sys/wait.h>
#include <unistd.h>

extern "C" {
extern int splittbfast_library(int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv,
                               int (*callback)(int, int, char *));
#define GUI_LENGTHOVER 2
}

namespace miint {

struct MafftAligner::Impl {
	Impl() {
	}
};

MafftAligner::MafftAligner() {
	impl_ = std::make_unique<Impl>();
}

MafftAligner::~MafftAligner() = default;
MafftAligner::MafftAligner(MafftAligner &&) noexcept = default;
MafftAligner &MafftAligner::operator=(MafftAligner &&) noexcept = default;

static char detect_sequence_type(const std::vector<std::string> &sequences) {
	int aa_only_count = 0;
	int total_count = 0;
	for (const auto &seq : sequences) {
		for (char c : seq) {
			char upper = static_cast<char>(toupper(c));
			total_count++;
			if (upper == 'E' || upper == 'F' || upper == 'I' || upper == 'L' || upper == 'P' || upper == 'Q') {
				aa_only_count++;
			}
		}
	}
	if (total_count > 0 && (double)aa_only_count / total_count > 0.05) {
		return 'p';
	}
	return 'd';
}

static void restore_case(const std::string &original, std::string &aligned) {
	size_t orig_idx = 0;
	for (size_t i = 0; i < aligned.size(); i++) {
		if (aligned[i] != '-') {
			if (orig_idx < original.size()) {
				aligned[i] = original[orig_idx++];
			}
		}
	}
}

// Write exactly len bytes to fd, retrying on short writes.
static bool write_all(int fd, const void *buf, size_t len) {
	const char *p = static_cast<const char *>(buf);
	while (len > 0) {
		ssize_t w = write(fd, p, len);
		if (w <= 0) {
			return false;
		}
		p += w;
		len -= w;
	}
	return true;
}

// Read exactly len bytes from fd, retrying on short reads.
static bool read_all(int fd, void *buf, size_t len) {
	char *p = static_cast<char *>(buf);
	while (len > 0) {
		ssize_t r = read(fd, p, len);
		if (r <= 0) {
			return false;
		}
		p += r;
		len -= r;
	}
	return true;
}

static bool write_string(int fd, const std::string &s) {
	uint32_t len = static_cast<uint32_t>(s.size());
	if (!write_all(fd, &len, sizeof(len))) {
		return false;
	}
	return len == 0 || write_all(fd, s.data(), len);
}

static std::string read_string(int fd) {
	uint32_t len = 0;
	if (!read_all(fd, &len, sizeof(len))) {
		throw std::runtime_error("MafftAligner: failed to read from child process");
	}
	std::string s(len, '\0');
	if (len > 0 && !read_all(fd, &s[0], len)) {
		throw std::runtime_error("MafftAligner: incomplete read from child process");
	}
	return s;
}

MafftAlignResult MafftAligner::align(const std::vector<std::string> &names, const std::vector<std::string> &comments,
                                     const std::vector<std::string> &sequences) {
	if (sequences.size() < 2) {
		throw std::invalid_argument("MafftAligner: at least 2 sequences required");
	}
	if (names.size() != sequences.size() || comments.size() != sequences.size()) {
		throw std::invalid_argument("MafftAligner: names, comments, and sequences must have the same size");
	}
	for (size_t i = 0; i < sequences.size(); i++) {
		if (sequences[i].empty()) {
			throw std::invalid_argument("MafftAligner: sequence " + std::to_string(i) + " is empty");
		}
	}

	int n = static_cast<int>(sequences.size());
	char seqtype = detect_sequence_type(sequences);

	size_t max_len = 0;
	for (const auto &s : sequences) {
		if (s.size() > max_len) {
			max_len = s.size();
		}
	}
	size_t lgui_sz = max_len * 3 + 1000;
	if (lgui_sz > static_cast<size_t>(INT_MAX)) {
		throw std::invalid_argument("MafftAligner: sequences too long (max aligned length would overflow)");
	}
	int lgui = static_cast<int>(lgui_sz);

	// Fork to isolate MAFFT's ~150 global variables.
	// splittbfast's internal memory management (recursive splitseq_mq with reallocs)
	// does not support repeated calls in the same process.
	int pipefd[2];
	if (pipe(pipefd) == -1) {
		throw std::runtime_error("MafftAligner: pipe() failed: " + std::string(strerror(errno)));
	}

	pid_t child = fork();
	if (child == -1) {
		close(pipefd[0]);
		close(pipefd[1]);
		throw std::runtime_error("MafftAligner: fork() failed: " + std::string(strerror(errno)));
	}

	if (child == 0) {
		// === CHILD PROCESS ===
		close(pipefd[0]);

		// Suppress MAFFT's stderr output
		FILE *devnull = fopen("/dev/null", "w");
		if (devnull) {
			dup2(fileno(devnull), STDERR_FILENO);
			fclose(devnull);
		}

		char **seqgui = static_cast<char **>(calloc(n, sizeof(char *)));
		char **namegui = static_cast<char **>(calloc(n, sizeof(char *)));
		for (int i = 0; i < n; i++) {
			seqgui[i] = static_cast<char *>(calloc(lgui + 1, sizeof(char)));
			namegui[i] = static_cast<char *>(calloc(1000, sizeof(char)));
			strcpy(seqgui[i], sequences[i].c_str());
			std::string full_name = names[i];
			if (!comments[i].empty()) {
				full_name += " " + comments[i];
			}
			strncpy(namegui[i], full_name.c_str(), 999);
			namegui[i][999] = '\0';
		}

		int argc_m = 12;
		char **argv_m = static_cast<char **>(calloc(argc_m, sizeof(char *)));
		for (int i = 0; i < argc_m; i++) {
			argv_m[i] = static_cast<char *>(calloc(100, sizeof(char)));
		}
		strcpy(argv_m[0], "splittbfast");
		strcpy(argv_m[1], seqtype == 'd' ? "-D" : "-P");
		strcpy(argv_m[2], "-f");
		strcpy(argv_m[3], "-1.53");
		strcpy(argv_m[4], "-Q");
		strcpy(argv_m[5], "100");
		strcpy(argv_m[6], "-h");
		strcpy(argv_m[7], "0");
		strcpy(argv_m[8], "-p");
		strcpy(argv_m[9], "50");
		strcpy(argv_m[10], "-s");
		strcpy(argv_m[11], "-1");

		int result = splittbfast_library(n, lgui, namegui, seqgui, argc_m, argv_m, nullptr);

		write_all(pipefd[1], &result, sizeof(result));
		if (result == 0) {
			for (int i = 0; i < n; i++) {
				write_string(pipefd[1], std::string(seqgui[i]));
			}
		}

		close(pipefd[1]);
		_exit(result);
	}

	// === PARENT PROCESS ===
	close(pipefd[1]);

	int result = -1;
	if (!read_all(pipefd[0], &result, sizeof(result))) {
		close(pipefd[0]);
		int status;
		waitpid(child, &status, 0);
		if (WIFSIGNALED(status)) {
			throw std::runtime_error("MafftAligner: child killed by signal " + std::to_string(WTERMSIG(status)));
		}
		throw std::runtime_error("MafftAligner: child process died unexpectedly");
	}

	MafftAlignResult align_result;
	if (result == 0) {
		align_result.names = names;
		align_result.comments = comments;

		for (int i = 0; i < n; i++) {
			std::string aligned = read_string(pipefd[0]);
			restore_case(sequences[i], aligned);
			align_result.sequences.push_back(std::move(aligned));
			align_result.original_lengths.push_back(static_cast<int>(sequences[i].size()));
		}
		align_result.aligned_length =
		    align_result.sequences.empty() ? 0 : static_cast<int>(align_result.sequences[0].size());
	}

	close(pipefd[0]);

	int status;
	waitpid(child, &status, 0);

	if (result != 0) {
		throw std::runtime_error("MAFFT alignment failed with error code " + std::to_string(result));
	}

	return align_result;
}

} // namespace miint
