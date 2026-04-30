#include "include/MafftAligner.hpp"
#include "mafft_api.h"
#include <climits>
#include <cstring>
#include <mutex>
#include <stdexcept>

namespace miint {

// Process-wide mutex: MAFFT uses ~150 global variables and is not reentrant.
static std::mutex g_mafft_mutex;

struct MafftAligner::Impl {};

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
                                     const std::vector<std::string> &sequences, int n_threads) {
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

	// Build C-array views over the input. The API copies these internally,
	// so it is safe for the underlying std::strings to outlive only the call.
	std::vector<std::string> full_names(n);
	std::vector<const char *> name_ptrs(n);
	std::vector<const char *> seq_ptrs(n);
	for (int i = 0; i < n; i++) {
		full_names[i] = names[i];
		if (!comments[i].empty()) {
			full_names[i] += " " + comments[i];
		}
		name_ptrs[i] = full_names[i].c_str();
		seq_ptrs[i] = sequences[i].c_str();
	}

	mafft_config_t cfg;
	mafft_config_init(&cfg);
	cfg.strategy = MAFFT_STRATEGY_AUTO;
	cfg.seqtype = (detect_sequence_type(sequences) == 'd') ? MAFFT_SEQ_DNA : MAFFT_SEQ_PROTEIN;
	cfg.n_threads = n_threads > 0 ? n_threads : 1;

	mafft_output_t *out = nullptr;
	int rc;
	std::string err_msg;
	std::string mafft_log;
	{
		std::lock_guard<std::mutex> lock(g_mafft_mutex);
		mafft_ctx_t *ctx = mafft_create(&cfg);
		if (!ctx) {
			throw std::runtime_error("MafftAligner: mafft_create failed");
		}
		rc = mafft_align(ctx, name_ptrs.data(), seq_ptrs.data(), n, &out, nullptr);
		if (rc != MAFFT_OK) {
			const char *e = mafft_last_error(ctx);
			err_msg = e ? e : mafft_strerror(rc);
			const char *log = mafft_ctx_log(ctx);
			if (log && *log) {
				mafft_log = log;
			}
		}
		mafft_destroy(ctx);
	}

	if (rc != MAFFT_OK) {
		if (out) {
			mafft_output_free(out);
		}
		std::string full = "MAFFT alignment failed: " + err_msg + " (code " + std::to_string(rc) + ")";
		if (!mafft_log.empty()) {
			full += "\n--- captured MAFFT log ---\n" + mafft_log + "--- end MAFFT log ---";
		} else {
			full += "\n(captured MAFFT log was empty)";
		}
		throw std::runtime_error(full);
	}

	MafftAlignResult result;
	result.names = names;
	result.comments = comments;
	result.aligned_length = out->aligned_len;
	result.sequences.reserve(n);
	result.original_lengths.reserve(n);
	for (int i = 0; i < n; i++) {
		std::string aligned(out->seqs[i]);
		restore_case(sequences[i], aligned);
		result.sequences.push_back(std::move(aligned));
		result.original_lengths.push_back(static_cast<int>(sequences[i].size()));
	}
	mafft_output_free(out);
	return result;
}

} // namespace miint
