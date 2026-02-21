#pragma once

#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <minimap2/minimap.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Subject sequence for indexing (subjects cannot be paired-end)
struct AlignmentSubject {
	std::string read_id;
	std::string sequence;

	// Length is computed from sequence
	size_t length() const {
		return sequence.size();
	}
};

// Configuration for alignment
struct Minimap2Config {
	std::string preset = "sr"; // Preset: sr, map-ont, map-pb, asm5, etc.
	int max_secondary = 5;     // Max secondary alignments per query
	bool eqx = true;           // Use =/X instead of M in CIGAR
	int k = 0;                 // k-mer size (0 = use preset default)
	int w = 0;                 // minimizer window (0 = use preset default)
};

// Custom deleter for minimap2 index
struct Minimap2IndexDeleter {
	void operator()(mm_idx_t *idx) const;
};

// Custom deleter for minimap2 thread buffer
struct Minimap2TbufDeleter {
	void operator()(mm_tbuf_t *tbuf) const;
};

// RAII wrapper types
using Minimap2IndexPtr = std::unique_ptr<mm_idx_t, Minimap2IndexDeleter>;
using Minimap2TbufPtr = std::unique_ptr<mm_tbuf_t, Minimap2TbufDeleter>;

// Shared, immutable minimap2 index for multi-thread-per-shard alignment.
// Multiple Minimap2Aligner instances can reference the same SharedMinimap2Index
// concurrently (each aligner has its own mm_tbuf_t).
class SharedMinimap2Index {
public:
	SharedMinimap2Index(const std::string &index_path, const Minimap2Config &config);
	~SharedMinimap2Index();

	// Non-copyable
	SharedMinimap2Index(const SharedMinimap2Index &) = delete;
	SharedMinimap2Index &operator=(const SharedMinimap2Index &) = delete;

	const mm_idx_t *index() const;
	const mm_mapopt_t &mapopt() const;
	const std::vector<std::string> &subject_names() const;

private:
	Minimap2IndexPtr index_;
	mm_mapopt_t mopt_;
	std::vector<std::string> subject_names_;
};

// Main aligner class.
// NOT thread-safe: each thread must have its own Minimap2Aligner instance.
// Multiple instances may share a SharedMinimap2Index concurrently, but the
// aligner itself (including its mm_tbuf_t) must not be used from multiple threads.
class Minimap2Aligner {
public:
	explicit Minimap2Aligner(const Minimap2Config &config);
	~Minimap2Aligner();

	// Non-copyable
	Minimap2Aligner(const Minimap2Aligner &) = delete;
	Minimap2Aligner &operator=(const Minimap2Aligner &) = delete;

	// Movable
	Minimap2Aligner(Minimap2Aligner &&) noexcept;
	Minimap2Aligner &operator=(Minimap2Aligner &&) noexcept;

	// Build index from multiple subjects
	void build_index(const std::vector<AlignmentSubject> &subjects);

	// Build index from a single subject (for per_subject_database mode)
	void build_single_index(const AlignmentSubject &subject);

	// Load index from .mmi file
	void load_index(const std::string &index_path);

	// Save current index to .mmi file
	void save_index(const std::string &output_path) const;

	// Check if file is a valid minimap2 index
	static bool is_index_file(const std::string &path);

	// Static helpers for option/index initialization (shared with SharedMinimap2Index)
	static void InitOptions(const Minimap2Config &config, mm_idxopt_t &iopt, mm_mapopt_t &mopt);
	static void LoadIndexFromFile(const std::string &path, const mm_idxopt_t &iopt, mm_idx_t *&out_idx,
	                              std::vector<std::string> &out_names);

	// Attach a shared index (clears any owned index)
	void attach_shared_index(std::shared_ptr<SharedMinimap2Index> shared_idx);
	// Detach the shared index (does not destroy it; other aligners may still reference it)
	void detach_shared_index();

	// Align queries against current index, append results to batch
	// Uses SequenceRecordBatch which matches read_fastx output schema
	void align(const SequenceRecordBatch &queries, SAMRecordBatch &output);

private:
	Minimap2Config config_;
	std::unique_ptr<mm_idxopt_t> iopt_;
	std::unique_ptr<mm_mapopt_t> mopt_;
	Minimap2IndexPtr index_;
	std::vector<std::string> subject_names_;            // For reference name lookup
	Minimap2TbufPtr tbuf_;                              // Reusable thread buffer
	std::shared_ptr<SharedMinimap2Index> shared_index_; // Shared index (mutually exclusive with index_)

	// Accessors that transparently pick shared or owned state
	const mm_idx_t *active_index() const;
	const mm_mapopt_t *active_mapopt() const;
	const std::vector<std::string> &active_subject_names() const;

	// Internal alignment functions
	void align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output);
	void align_paired(const std::string &read_id, const std::string &sequence1, const std::string &sequence2,
	                  SAMRecordBatch &output);

	// Convert minimap2 result to SAM fields
	void reg_to_sam(const void *reg_ptr, const std::string &read_id, const std::string &query_seq,
	                SAMRecordBatch &batch, int segment_idx, bool mate_mapped, bool mate_rev, int32_t mate_rid,
	                int32_t mate_pos, int32_t tlen);

	// Generate CIGAR string from mm_extra_t
	std::string cigar_string(const void *reg_ptr) const;

	// Calculate stop position from CIGAR
	int64_t calculate_stop_position(int64_t start_pos, const void *reg_ptr) const;

	// Calculate SAM flags
	uint16_t calculate_flags(const void *reg_ptr, bool is_paired, int segment_idx, bool mate_mapped, bool mate_rev,
	                         bool is_unmapped) const;

	// Get reference name by ID
	const std::string &get_reference_name(int32_t rid) const;
};

} // namespace miint
