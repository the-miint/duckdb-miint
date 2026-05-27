#pragma once

#include <string>
#include <vector>

namespace miint {

struct AbpoaAlignParams {
	int match = 2;
	int mismatch = 4;
	int gap_open1 = 4;
	int gap_open2 = 24;
	int gap_ext1 = 2;
	int gap_ext2 = 1;
	std::string align_mode = "global";
	bool progressive = true;
	bool disable_seeding = false;
	bool amb_strand = false;
	int k = 19;
	int w = 10;
	int min_w = 500;
	int bandwidth = 10;
	float bandwidth_frac = 0.01f;
	int max_num_cons = 1;
	float min_freq = 0.25f;
	std::string algorithm = "heaviest_bundling";
};

struct AbpoaMsaResult {
	std::vector<std::string> names;
	std::vector<std::string> aligned_sequences;
	std::vector<int> original_lengths;
	int aligned_length;
};

struct AbpoaConsensusEntry {
	int consensus_id;
	std::string sequence;
	int length;
	int num_reads;
};

struct AbpoaConsensusResult {
	std::vector<AbpoaConsensusEntry> entries;
};

class AbpoaAligner {
public:
	AbpoaAligner();
	~AbpoaAligner();
	AbpoaAligner(const AbpoaAligner &) = delete;
	AbpoaAligner &operator=(const AbpoaAligner &) = delete;

	AbpoaMsaResult align(const std::vector<std::string> &names, const std::vector<std::string> &sequences,
	                      const AbpoaAlignParams &params = {});

	AbpoaConsensusResult consensus(const std::vector<std::string> &names,
	                               const std::vector<std::string> &sequences,
	                               const AbpoaAlignParams &params = {});
};

} // namespace miint
