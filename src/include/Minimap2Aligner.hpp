#pragma once

#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <minimap2/minimap.h>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Alignment stats computed during CIGAR string generation
struct AlignmentStats {
	int64_t mismatches = 0;  // XM: number of mismatches
	int64_t gap_opens = 0;   // XO: number of gap opens
	int64_t gap_extends = 0; // XG: number of gap extensions
};

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
	std::string preset = "sr";       // Preset: sr, map-ont, map-pb, asm5, etc.
	int max_secondary = 5;           // Max secondary alignments per query
	bool eqx = true;                 // Use =/X instead of M in CIGAR
	int k = 0;                       // k-mer size (0 = use preset default)
	int w = 0;                       // minimizer window (0 = use preset default)
	float min_chain_coverage = 0.0f; // min best-chain span (qe-qs)/qlen to attempt DP (0.0 = disabled)

	// High-occurrence minimizer filter, minimap2's -f. Parsed at bind into minimap2's own three
	// variables so InitOptions is a direct transcription; occ_filter_set == false leaves whatever
	// the preset chose. occ_mid_frac defaults to minimap2's own default so that the absolute form
	// ("-f 100000") leaves it untouched exactly as the CLI does.
	bool occ_filter_set = false;
	float occ_mid_frac = 2e-4f; // -> mm_mapopt_t::mid_occ_frac; only read when occ_mid <= 0
	int32_t occ_mid = 0;        // -> mm_mapopt_t::mid_occ; <= 0 makes mm_mapopt_update derive it
	// -> mm_mapopt_t::max_occ. -1, not 0, means "no second value was given": minimap2 accepts
	// `-f N,0` and 0 there is meaningful (it disables the re-chain pass), so the two cannot share
	// a sentinel.
	int32_t occ_max = -1;

	// Emit one row per query that produced no alignment, instead of nothing at all (#185)
	bool include_unmapped = false;
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
	// Load from .mmi file
	SharedMinimap2Index(const std::string &index_path, const Minimap2Config &config);
	// Take ownership of a pre-built index. Takes Minimap2IndexPtr (not a raw
	// mm_idx_t*) so ownership transfers atomically with argument binding: a
	// caller building this via std::make_shared can pass std::move(idx) and be
	// certain the index is freed even if make_shared's own allocation throws
	// before this constructor ever runs — a raw pointer argument would leave
	// nothing owning the index during that window.
	SharedMinimap2Index(Minimap2IndexPtr idx, const mm_mapopt_t &mopt, std::vector<std::string> subject_names);
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

// Streams a possibly multi-part .mmi file one part at a time.
//
// minimap2's own answer to "reference larger than RAM" is the multi-part index
// (built with `minimap2 -d out.mmi -I <batch_size> ref.fa`): mm_idx_reader_read
// returns one part per call, and the CLI loops until mm_idx_reader_eof, aligning
// every read against each part in turn. Peak memory is one part, not the whole
// index. This class is that loop, wrapped so each part comes out as an
// independently owned SharedMinimap2Index that can be dropped before the next
// part is loaded.
//
// NOT thread-safe: callers (align_minimap2's part-transition logic) serialize
// access via their own lock, since only one thread should ever be advancing the
// underlying mm_idx_reader_t at a time.
class Minimap2IndexReader {
public:
	Minimap2IndexReader(const std::string &index_path, const Minimap2Config &config);
	~Minimap2IndexReader();

	Minimap2IndexReader(const Minimap2IndexReader &) = delete;
	Minimap2IndexReader &operator=(const Minimap2IndexReader &) = delete;

	// Reads and returns the next part, or nullptr once the file is exhausted.
	// mm_mapopt_update is applied against THIS part's index (mid_occ is derived
	// from the loaded index's minimizer distribution, so reusing an earlier
	// part's value would apply the wrong high-occurrence filter).
	std::shared_ptr<SharedMinimap2Index> ReadNextPart();

	// True if there is no next part. Confirms this by actually attempting to
	// read it (same reasoning as LoadIndexFromFile: mm_idx_reader_eof's
	// file-position heuristic can report "not eof" for a single-part file that
	// has trailing bytes for an unrelated reason, e.g. a padded transfer or an
	// appended sidecar). Unlike a naive peek-and-cache, the confirming read is
	// immediately destroyed and the file position rewound rather than retained:
	// keeping the loaded part alive here would sit it resident for the entire
	// time the PREVIOUS part is being aligned against, doubling peak memory for
	// every multi-part index. ReadNextPart() re-reads it for real, for keeps,
	// only once the caller actually asks for it.
	bool AtEof();

private:
	mm_idx_reader_t *reader_ = nullptr;
	// Built once at construction (preset/k/w parsed and validated here, not
	// per part) and copied into a fresh mm_mapopt_t for each ReadNextPart call,
	// which then runs only mm_mapopt_update against that part's own minimizer
	// distribution — the one piece of mopt that must be re-derived per part.
	// mm_idxopt_t is NOT cached the same way: mm_idx_reader_open copies it into
	// the reader's own state at open time and never consults the caller's copy
	// again, so keeping it as a member here would just be dead weight.
	mm_mapopt_t mopt_template_;
	Minimap2Config config_;
	// Set by AtEof() once it has confirmed whether a next part exists, so a
	// second call doesn't re-probe. Deliberately does NOT cache the part
	// itself — see AtEof().
	bool has_peeked_ = false;
	bool next_part_exists_ = false;

	// The actual read, shared by AtEof()'s confirmation read and ReadNextPart()'s
	// direct read when nothing is pending from a previous peek.
	std::shared_ptr<SharedMinimap2Index> ReadNextPartUncached();
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

	// Build a SharedMinimap2Index from subjects (for multi-threaded standard mode)
	static std::shared_ptr<SharedMinimap2Index> BuildSharedIndex(const std::vector<AlignmentSubject> &subjects,
	                                                             const Minimap2Config &config);

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
	char *md_buf_ = nullptr;                            // Pooled buffer for mm_gen_MD (Opt 5)
	int md_max_len_ = 0;                                // Current capacity of md_buf_

	// Build raw mm_idx_t* from subjects (shared by build_index and BuildSharedIndex)
	static mm_idx_t *BuildRawIndex(const std::vector<AlignmentSubject> &subjects, const mm_idxopt_t &iopt,
	                               std::vector<std::string> &out_names);

	// Accessors that transparently pick shared or owned state
	const mm_idx_t *active_index() const;
	const mm_mapopt_t *active_mapopt() const;
	const std::vector<std::string> &active_subject_names() const;

	// Internal alignment functions
	void align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output);
	void align_paired(const std::string &read_id, const std::string &sequence1, const std::string &sequence2,
	                  SAMRecordBatch &output);

	// Convert minimap2 result to SAM fields
	void reg_to_sam(const mm_reg1_t *reg, const std::string &read_id, const std::string &query_seq,
	                SAMRecordBatch &batch, int segment_idx, bool mate_mapped, bool mate_rev, int32_t mate_rid,
	                int32_t mate_pos, int32_t tlen);

	// Append a synthetic row for a query (or paired segment) that produced no alignment, so that
	// "no alignment" is representable in the result set rather than inferable only from a missing
	// row (#185). segment_idx is -1 for single-end. Only called when config_.include_unmapped.
	// mate_rid < 0 means the mate did not map; when it did, its coordinates are carried onto this
	// row as RNEXT/PNEXT, which SAM requires of a paired record whose mate is mapped.
	void append_unmapped(const std::string &read_id, int segment_idx, bool mate_mapped, bool mate_rev, int32_t mate_rid,
	                     int32_t mate_pos, SAMRecordBatch &batch);

	// Generate CIGAR string from mm_extra_t, including soft/hard clips.
	// When stats_out is non-null, computes XM/XO/XG during the same CIGAR walk.
	std::string cigar_string(const mm_reg1_t *reg, int32_t query_len, uint16_t sam_flags,
	                         AlignmentStats *stats_out = nullptr) const;

	// Calculate stop position from CIGAR
	int64_t calculate_stop_position(int64_t start_pos, const mm_reg1_t *reg) const;

	// Calculate SAM flags
	uint16_t calculate_flags(const mm_reg1_t *reg, bool is_paired, int segment_idx, bool mate_mapped, bool mate_rev,
	                         bool is_unmapped) const;

	// Get reference name by ID
	const std::string &get_reference_name(int32_t rid) const;
};

} // namespace miint
