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

// Main aligner class
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

	// Align queries against current index, append results to batch
	// Uses SequenceRecordBatch which matches read_fastx output schema
	void align(const SequenceRecordBatch &queries, SAMRecordBatch &output);

private:
	Minimap2Config config_;
	std::unique_ptr<mm_idxopt_t> iopt_;
	std::unique_ptr<mm_mapopt_t> mopt_;
	Minimap2IndexPtr index_;
	std::vector<std::string> subject_names_; // For reference name lookup
	Minimap2TbufPtr tbuf_;                   // Reusable thread buffer

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
