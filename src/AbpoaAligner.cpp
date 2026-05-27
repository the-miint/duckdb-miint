#include "AbpoaAligner.hpp"
#include <stdexcept>
#include <cstring>

extern "C" {
#include "abpoa.h"
extern unsigned char ab_nt4_table[256];
extern const char ab_nt256_table[256];
}

namespace miint {

AbpoaAligner::AbpoaAligner() {
}

AbpoaAligner::~AbpoaAligner() {
}

static int ParseAlignMode(const std::string &mode) {
	if (mode == "global") {
		return ABPOA_GLOBAL_MODE;
	}
	if (mode == "local") {
		return ABPOA_LOCAL_MODE;
	}
	if (mode == "extension") {
		return ABPOA_EXTEND_MODE;
	}
	throw std::invalid_argument("Invalid align_mode '" + mode + "': must be 'global', 'local', or 'extension'");
}

static int ParseConsAlgorithm(const std::string &alg) {
	if (alg == "heaviest_bundling") {
		return ABPOA_HB;
	}
	if (alg == "majority_voting") {
		return ABPOA_MF;
	}
	throw std::invalid_argument("Invalid algorithm '" + alg + "': must be 'heaviest_bundling' or 'majority_voting'");
}

static void ConfigureParams(abpoa_para_t *abpt, const AbpoaAlignParams &params, bool want_msa, bool want_cons) {
	abpt->match = params.match;
	abpt->mismatch = params.mismatch;
	abpt->gap_open1 = params.gap_open1;
	abpt->gap_open2 = params.gap_open2;
	abpt->gap_ext1 = params.gap_ext1;
	abpt->gap_ext2 = params.gap_ext2;
	abpt->align_mode = ParseAlignMode(params.align_mode);
	abpt->progressive_poa = params.progressive ? 1 : 0;
	abpt->disable_seeding = params.disable_seeding ? 1 : 0;
	abpt->amb_strand = params.amb_strand ? 1 : 0;
	abpt->k = params.k;
	abpt->w = params.w;
	abpt->min_w = params.min_w;
	abpt->wb = params.bandwidth;
	abpt->wf = params.bandwidth_frac;
	abpt->max_n_cons = params.max_num_cons;
	abpt->min_freq = params.min_freq;
	abpt->cons_algrm = ParseConsAlgorithm(params.algorithm);
	abpt->out_msa = want_msa ? 1 : 0;
	abpt->out_cons = want_cons ? 1 : 0;
	abpt->out_gfa = 0;
	abpt->use_read_ids = (want_msa || params.max_num_cons > 1) ? 1 : 0;
	abpoa_post_set_para(abpt);
}

struct AbpoaGuard {
	abpoa_t *ab = nullptr;
	abpoa_para_t *abpt = nullptr;
	uint8_t **seqs = nullptr;
	int *seq_lens = nullptr;
	char **seq_names = nullptr;
	int n_seq = 0;

	~AbpoaGuard() {
		if (seqs) {
			for (int i = 0; i < n_seq; i++) {
				free(seqs[i]);
			}
			free(seqs);
		}
		if (seq_lens) {
			free(seq_lens);
		}
		if (seq_names) {
			for (int i = 0; i < n_seq; i++) {
				free(seq_names[i]);
			}
			free(seq_names);
		}
		if (ab) {
			abpoa_free(ab);
		}
		if (abpt) {
			abpoa_free_para(abpt);
		}
	}
};

static void ValidateInputs(const std::vector<std::string> &names, const std::vector<std::string> &sequences) {
	if (names.size() != sequences.size()) {
		throw std::invalid_argument("names and sequences must have the same size");
	}
	if (sequences.size() < 2) {
		throw std::invalid_argument("at least 2 sequences required");
	}
	for (size_t i = 0; i < sequences.size(); i++) {
		if (sequences[i].empty()) {
			throw std::invalid_argument("sequence at index " + std::to_string(i) + " is empty");
		}
	}
}

static void PrepareSequences(AbpoaGuard &guard, const std::vector<std::string> &names,
                             const std::vector<std::string> &sequences) {
	int n = static_cast<int>(sequences.size());
	guard.seqs = (uint8_t **)calloc(n, sizeof(uint8_t *));
	guard.seq_lens = (int *)calloc(n, sizeof(int));
	guard.seq_names = (char **)calloc(n, sizeof(char *));
	if (!guard.seqs || !guard.seq_lens || !guard.seq_names) {
		throw std::runtime_error("failed to allocate sequence arrays");
	}

	for (int i = 0; i < n; i++) {
		guard.n_seq = i + 1;
		int len = static_cast<int>(sequences[i].size());
		guard.seq_lens[i] = len;
		guard.seqs[i] = (uint8_t *)malloc(len * sizeof(uint8_t));
		if (!guard.seqs[i]) {
			throw std::runtime_error("failed to allocate sequence buffer");
		}
		for (int j = 0; j < len; j++) {
			guard.seqs[i][j] = ab_nt4_table[(unsigned char)sequences[i][j]];
		}
		guard.seq_names[i] = strdup(names[i].c_str());
	}
}

AbpoaMsaResult AbpoaAligner::align(const std::vector<std::string> &names, const std::vector<std::string> &sequences,
                                    const AbpoaAlignParams &params) {
	ValidateInputs(names, sequences);

	AbpoaGuard guard;
	guard.abpt = abpoa_init_para();
	ConfigureParams(guard.abpt, params, true, false);
	guard.ab = abpoa_init();
	PrepareSequences(guard, names, sequences);

	(void)abpoa_msa(guard.ab, guard.abpt, guard.n_seq, guard.seq_names, guard.seq_lens, guard.seqs, nullptr, nullptr);

	abpoa_cons_t *abc = guard.ab->abc;
	if (!abc || abc->msa_len <= 0) {
		throw std::runtime_error("abPOA produced no MSA output");
	}

	AbpoaMsaResult result;
	result.aligned_length = abc->msa_len;
	int n = guard.n_seq;
	result.names.reserve(n);
	result.aligned_sequences.reserve(n);
	result.original_lengths.reserve(n);

	for (int i = 0; i < n; i++) {
		result.names.push_back(names[i]);
		std::string aligned(abc->msa_len, ' ');
		for (int j = 0; j < abc->msa_len; j++) {
			uint8_t code = abc->msa_base[i][j];
			aligned[j] = ab_nt256_table[code];
		}
		result.aligned_sequences.push_back(std::move(aligned));
		result.original_lengths.push_back(static_cast<int>(sequences[i].size()));
	}

	return result;
}

AbpoaConsensusResult AbpoaAligner::consensus(const std::vector<std::string> &names,
                                              const std::vector<std::string> &sequences,
                                              const AbpoaAlignParams &params) {
	ValidateInputs(names, sequences);

	AbpoaGuard guard;
	guard.abpt = abpoa_init_para();
	ConfigureParams(guard.abpt, params, false, true);
	guard.ab = abpoa_init();
	PrepareSequences(guard, names, sequences);

	(void)abpoa_msa(guard.ab, guard.abpt, guard.n_seq, guard.seq_names, guard.seq_lens, guard.seqs, nullptr, nullptr);

	abpoa_cons_t *abc = guard.ab->abc;
	if (!abc || abc->n_cons <= 0) {
		throw std::runtime_error("abPOA produced no consensus output");
	}

	AbpoaConsensusResult result;
	for (int i = 0; i < abc->n_cons; i++) {
		AbpoaConsensusEntry entry;
		entry.consensus_id = i;
		entry.length = abc->cons_len[i];
		entry.num_reads = abc->clu_n_seq ? abc->clu_n_seq[i] : guard.n_seq;

		std::string seq(entry.length, ' ');
		for (int j = 0; j < entry.length; j++) {
			uint8_t code = abc->cons_base[i][j];
			seq[j] = ab_nt256_table[code];
		}
		entry.sequence = std::move(seq);
		result.entries.push_back(std::move(entry));
	}

	return result;
}

} // namespace miint
