#pragma once
#include <memory>
#include <string>
#include <vector>

namespace miint {

struct MafftAlignResult {
	std::vector<std::string> names;
	std::vector<std::string> comments;
	std::vector<std::string> sequences;
	std::vector<int> original_lengths;
	int aligned_length;
};

class MafftAligner {
public:
	MafftAligner();
	~MafftAligner();

	MafftAligner(const MafftAligner &) = delete;
	MafftAligner &operator=(const MafftAligner &) = delete;
	MafftAligner(MafftAligner &&) noexcept;
	MafftAligner &operator=(MafftAligner &&) noexcept;

	// Align sequences using MAFFT (auto strategy: FFT-NS-2 for nseq < 200k,
	// PartTree above). n_threads is forwarded to MAFFT's internal pthread
	// pool; use 1 to disable. Thread safety: protected by internal mutex —
	// only one alignment may run at a time process-wide because MAFFT uses
	// ~150 globals.
	MafftAlignResult align(const std::vector<std::string> &names, const std::vector<std::string> &comments,
	                       const std::vector<std::string> &sequences, int n_threads = 1);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
