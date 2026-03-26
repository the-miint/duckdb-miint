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

	// Align sequences using MAFFT PartTree algorithm.
	// Thread safety: protected by internal mutex (single alignment at a time, process-wide).
	MafftAlignResult align(const std::vector<std::string> &names, const std::vector<std::string> &comments,
	                       const std::vector<std::string> &sequences);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

} // namespace miint
