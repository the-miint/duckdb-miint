#pragma once

#include <cstdint>
#include <string>
#include <vector>

struct smr_context;
struct smr_index;

namespace miint {

struct SortMeRNAConfig {
	int32_t num_threads = 0;
	int32_t match = 2;
	int32_t mismatch = -3;
	int32_t gap_open = 5;
	int32_t gap_ext = 2;
	int32_t score_N = -3;
	double evalue = -1.0;
	uint32_t seed_win_len = 18;
	uint32_t num_alignments = 1;
	bool best = true;
	bool paired = false;
	bool forward_only = false;
	bool reverse_only = false;
	bool full_search = false;
};

// Minimal SoA query container. Phase 3 adapters map SequenceRecordBatch → this.
// For paired-end, sequences2 must be the same length as sequences.
struct SortMeRNAQueryBatch {
	std::vector<std::string> read_ids;
	std::vector<std::string> sequences;
	std::vector<std::string> sequences2;

	size_t size() const {
		return read_ids.size();
	}
	bool empty() const {
		return read_ids.empty();
	}
	bool is_paired() const {
		return !sequences2.empty();
	}
};

// SoA result container. Parallel to SmrOutput but with C++ strings and the
// caller-supplied read_id written back by position (no library strdup kept).
// For paired-end, each input row produces two result rows with segment_idx
// 0 (fwd) and 1 (rev).
struct SortMeRNAResultBatch {
	std::vector<std::string> read_ids;
	std::vector<int32_t> aligned;       // 1 or 0
	std::vector<int32_t> strands;       // 1=fwd, 0=rc, -1=unaligned
	std::vector<std::string> ref_names; // "" if unaligned
	std::vector<int32_t> ref_starts;    // 1-based, 0 if unaligned
	std::vector<int32_t> ref_ends;      // 1-based, 0 if unaligned
	std::vector<std::string> cigars;    // "" if unaligned
	std::vector<int32_t> scores;        // -1 if unaligned
	std::vector<double> e_values;
	std::vector<double> identities;       // 0-100
	std::vector<double> coverages;        // 0-100
	std::vector<int32_t> edit_distances;  // -1 if unaligned
	std::vector<int32_t> segment_indices; // 0=single/fwd, 1=rev

	size_t size() const {
		return read_ids.size();
	}
	bool empty() const {
		return read_ids.empty();
	}
	void clear();
	void reserve(size_t n);
};

class SortMeRNAAligner {
public:
	SortMeRNAAligner(const SortMeRNAConfig &cfg, const std::vector<std::string> &ref_paths);
	~SortMeRNAAligner();

	SortMeRNAAligner(const SortMeRNAAligner &) = delete;
	SortMeRNAAligner &operator=(const SortMeRNAAligner &) = delete;

	// Align a sub-batch of queries and append results to `out`. Thread-safety:
	// sortmerna's g_run_mutex serialises concurrent callers; the wrapper adds
	// no additional synchronisation.
	void align(const SortMeRNAQueryBatch &queries, SortMeRNAResultBatch &out);

	static std::string version();

private:
	smr_context *ctx_ = nullptr;
	smr_index *idx_ = nullptr;
	bool paired_ = false;
};

} // namespace miint
