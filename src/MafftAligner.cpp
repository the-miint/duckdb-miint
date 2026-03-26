#include "include/MafftAligner.hpp"
#include <climits>
#include <cstring>
#include <mutex>
#include <stdexcept>

extern "C" {
extern int splittbfast_library(int ngui, int lgui, char **namegui, char **seqgui, int argc, char **argv,
                               int (*callback)(int, int, char *));
}

namespace miint {

// Process-wide mutex: MAFFT uses ~150 global variables and is not reentrant.
static std::mutex g_mafft_mutex;

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

// Detect whether sequences are DNA or protein.
// Scans for amino-acid-only characters (E, F, I, L, P, Q) not in the DNA alphabet.
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

// Restore original case after MAFFT lowercases DNA or uppercases protein.
// Same algorithm as restoreu.c:fillorichar()
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
		// MAFFT's internal code calls ErrorExit() (which invokes exit(1)) for
		// sequences shorter than 6 characters. Validate here to throw a C++
		// exception instead of killing the process.
		if (sequences[i].size() < 6) {
			throw std::invalid_argument("MafftAligner: sequence " + std::to_string(i) + " too short (" +
			                            std::to_string(sequences[i].size()) + " chars, minimum 6)");
		}
	}

	int n = static_cast<int>(sequences.size());

	size_t max_len = 0;
	for (const auto &s : sequences) {
		if (s.size() > max_len) {
			max_len = s.size();
		}
	}

	// lgui = buffer size per sequence (must hold aligned result with gaps)
	size_t lgui_sz = max_len * 3 + 1000;
	if (lgui_sz > static_cast<size_t>(INT_MAX)) {
		throw std::invalid_argument("MafftAligner: sequences too long (max aligned length would overflow)");
	}
	int lgui = static_cast<int>(lgui_sz);

	char seqtype = detect_sequence_type(sequences);

	// Allocate C arrays for MAFFT
	char **seqgui = static_cast<char **>(calloc(n, sizeof(char *)));
	char **namegui = static_cast<char **>(calloc(n, sizeof(char *)));
	if (!seqgui || !namegui) {
		free(seqgui);
		free(namegui);
		throw std::runtime_error("MafftAligner: allocation failed");
	}
	for (int i = 0; i < n; i++) {
		seqgui[i] = static_cast<char *>(calloc(lgui + 1, sizeof(char)));
		namegui[i] = static_cast<char *>(calloc(1000, sizeof(char)));
		strcpy(seqgui[i], sequences[i].c_str());
		// copydatafromgui() expects plain name (it prepends '=')
		std::string full_name = names[i];
		if (!comments[i].empty()) {
			full_name += " " + comments[i];
		}
		strncpy(namegui[i], full_name.c_str(), 999);
		namegui[i][999] = '\0';
	}

	// Construct argv for splittbfast --parttree mode
	int argc_m = 12;
	char **argv_m = static_cast<char **>(calloc(argc_m, sizeof(char *)));
	for (int i = 0; i < argc_m; i++) {
		argv_m[i] = static_cast<char *>(calloc(100, sizeof(char)));
	}
	strcpy(argv_m[0], "splittbfast");
	strcpy(argv_m[1], seqtype == 'd' ? "-D" : "-P");
	strcpy(argv_m[2], "-f");
	strcpy(argv_m[3], "-1.53"); // gap open penalty
	strcpy(argv_m[4], "-Q");
	strcpy(argv_m[5], "100"); // spfactor
	strcpy(argv_m[6], "-h");
	strcpy(argv_m[7], "0"); // aof
	strcpy(argv_m[8], "-p");
	strcpy(argv_m[9], "50"); // partsize (PICKSIZE)
	strcpy(argv_m[10], "-s");
	strcpy(argv_m[11], "-1"); // groupsize = njob+1 (full alignment)

	// Call MAFFT (mutex-protected — MAFFT uses global state)
	int result;
	{
		std::lock_guard<std::mutex> lock(g_mafft_mutex);
		result = splittbfast_library(n, lgui, namegui, seqgui, argc_m, argv_m, nullptr);
	}

	// Build result before cleanup (seqgui holds aligned sequences)
	MafftAlignResult align_result;
	if (result == 0) {
		align_result.names = names;
		align_result.comments = comments;
		align_result.aligned_length = static_cast<int>(strlen(seqgui[0]));

		for (int i = 0; i < n; i++) {
			std::string aligned(seqgui[i]);
			restore_case(sequences[i], aligned);
			align_result.sequences.push_back(std::move(aligned));
			align_result.original_lengths.push_back(static_cast<int>(sequences[i].size()));
		}
	}

	// Cleanup C arrays
	for (int i = 0; i < n; i++) {
		free(seqgui[i]);
		free(namegui[i]);
	}
	free(seqgui);
	free(namegui);
	for (int i = 0; i < argc_m; i++) {
		free(argv_m[i]);
	}
	free(argv_m);

	if (result != 0) {
		throw std::runtime_error("MAFFT alignment failed with error code " + std::to_string(result));
	}

	return align_result;
}

} // namespace miint
