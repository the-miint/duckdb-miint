#include "Minimap2Aligner.hpp"
#include <minimap2/minimap.h>
#include <minimap2/mmpriv.h>
#include <cstdlib>
#include <sstream>
#include <stdexcept>

namespace miint {

// Custom deleters
void Minimap2IndexDeleter::operator()(mm_idx_t *idx) const {
	if (idx) {
		mm_idx_destroy(idx);
	}
}

void Minimap2TbufDeleter::operator()(mm_tbuf_t *tbuf) const {
	if (tbuf) {
		mm_tbuf_destroy(tbuf);
	}
}

// Helper struct to hold alignment stats computed from CIGAR
struct AlignmentStats {
	int64_t mismatches = 0;  // XM: number of mismatches
	int64_t gap_opens = 0;   // XO: number of gap opens
	int64_t gap_extends = 0; // XG: number of gap extensions
	int64_t edit_distance = 0; // NM: total edit distance
};

// Compute alignment statistics from CIGAR
static AlignmentStats compute_alignment_stats(const mm_reg1_t *reg) {
	AlignmentStats stats;

	if (!reg->p || reg->p->n_cigar == 0) {
		return stats;
	}

	for (uint32_t i = 0; i < reg->p->n_cigar; i++) {
		uint32_t op = reg->p->cigar[i] & 0xf;
		uint32_t len = reg->p->cigar[i] >> 4;

		switch (op) {
		case MM_CIGAR_X_MISMATCH: // X: sequence mismatch
			stats.mismatches += len;
			stats.edit_distance += len;
			break;
		case MM_CIGAR_INS: // I: insertion to reference
			stats.gap_opens += 1;
			stats.gap_extends += (len > 1) ? (len - 1) : 0;
			stats.edit_distance += len;
			break;
		case MM_CIGAR_DEL: // D: deletion from reference
			stats.gap_opens += 1;
			stats.gap_extends += (len > 1) ? (len - 1) : 0;
			stats.edit_distance += len;
			break;
		case MM_CIGAR_EQ_MATCH: // =: sequence match
		case MM_CIGAR_MATCH:    // M: alignment match (can be match or mismatch)
		case MM_CIGAR_N_SKIP:   // N: skipped region (intron)
		case MM_CIGAR_SOFTCLIP: // S: soft clipping
		case MM_CIGAR_HARDCLIP: // H: hard clipping
		case MM_CIGAR_PADDING:  // P: padding
		default:
			break;
		}
	}

	return stats;
}

// Constructor
Minimap2Aligner::Minimap2Aligner(const Minimap2Config &config)
    : config_(config), iopt_(std::make_unique<mm_idxopt_t>()), mopt_(std::make_unique<mm_mapopt_t>()),
      tbuf_(mm_tbuf_init()) {

	// Initialize with default values first (required before mm_set_opt)
	mm_idxopt_init(iopt_.get());
	mm_mapopt_init(mopt_.get());

	// Apply preset-specific options
	if (mm_set_opt(config_.preset.c_str(), iopt_.get(), mopt_.get()) != 0) {
		throw std::runtime_error("Unknown minimap2 preset: " + config_.preset);
	}

	// Override k and w if specified
	if (config_.k > 0) {
		iopt_->k = static_cast<short>(config_.k);
	}
	if (config_.w > 0) {
		iopt_->w = static_cast<short>(config_.w);
	}

	// Enable CIGAR output
	mopt_->flag |= MM_F_CIGAR;

	// Enable EQX mode (=/X instead of M)
	if (config_.eqx) {
		mopt_->flag |= MM_F_EQX;
	}

	// Enable MD tag output
	mopt_->flag |= MM_F_OUT_MD;

	// Set max secondary alignments
	// best_n controls how many chains go to DP alignment
	mopt_->best_n = config_.max_secondary + 1;
}

// Destructor
Minimap2Aligner::~Minimap2Aligner() = default;

// Move constructor
Minimap2Aligner::Minimap2Aligner(Minimap2Aligner &&other) noexcept
    : config_(std::move(other.config_)), iopt_(std::move(other.iopt_)), mopt_(std::move(other.mopt_)),
      index_(std::move(other.index_)), subject_names_(std::move(other.subject_names_)),
      tbuf_(std::move(other.tbuf_)) {
}

// Move assignment
Minimap2Aligner &Minimap2Aligner::operator=(Minimap2Aligner &&other) noexcept {
	if (this != &other) {
		config_ = std::move(other.config_);
		iopt_ = std::move(other.iopt_);
		mopt_ = std::move(other.mopt_);
		index_ = std::move(other.index_);
		subject_names_ = std::move(other.subject_names_);
		tbuf_ = std::move(other.tbuf_);
	}
	return *this;
}

void Minimap2Aligner::build_index(const std::vector<AlignmentSubject> &subjects) {
	if (subjects.empty()) {
		throw std::runtime_error("Cannot build index from empty subject list");
	}

	// Prepare sequence and name arrays for mm_idx_str
	std::vector<const char *> seqs;
	std::vector<const char *> names;
	subject_names_.clear();

	seqs.reserve(subjects.size());
	names.reserve(subjects.size());
	subject_names_.reserve(subjects.size());

	for (const auto &subject : subjects) {
		seqs.push_back(subject.sequence.c_str());
		names.push_back(subject.read_id.c_str());
		subject_names_.push_back(subject.read_id);
	}

	// Build index using mm_idx_str
	mm_idx_t *idx = mm_idx_str(iopt_->w, iopt_->k,
	                           0, // is_hpc
	                           iopt_->bucket_bits, static_cast<int>(subjects.size()), seqs.data(), names.data());

	if (!idx) {
		throw std::runtime_error("Failed to build minimap2 index");
	}

	index_.reset(idx);

	// Update mapping options based on index
	mm_mapopt_update(mopt_.get(), index_.get());
}

void Minimap2Aligner::build_single_index(const AlignmentSubject &subject) {
	std::vector<AlignmentSubject> subjects = {subject};
	build_index(subjects);
}

void Minimap2Aligner::align(const SequenceRecordBatch &queries, SAMRecordBatch &output) {
	if (queries.empty()) {
		return;
	}

	if (!index_) {
		throw std::runtime_error("No index built. Call build_index() first.");
	}

	// Process each query
	for (size_t i = 0; i < queries.size(); i++) {
		if (queries.is_paired) {
			align_paired(queries.read_ids[i], queries.sequences1[i], queries.sequences2[i], output);
		} else {
			align_single(queries.read_ids[i], queries.sequences1[i], output);
		}
	}
}

void Minimap2Aligner::align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output) {
	int n_regs = 0;
	mm_reg1_t *regs =
	    mm_map(index_.get(), static_cast<int>(sequence.length()), sequence.c_str(), &n_regs, tbuf_.get(), mopt_.get(),
	           read_id.c_str());

	int secondary_count = 0;

	// Process alignments (minimap2 returns them sorted by score, best first)
	for (int j = 0; j < n_regs; j++) {
		mm_reg1_t *reg = &regs[j];

		// Bounds check for reference ID before using it
		if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= subject_names_.size()) {
			continue; // Skip alignments with invalid reference ID
		}

		// Check if this is primary or secondary
		bool is_primary = (reg->parent == reg->id);
		bool is_secondary = !is_primary;

		// Limit secondary alignments
		if (is_secondary) {
			if (secondary_count >= config_.max_secondary) {
				continue;
			}
			secondary_count++;
		}

		// Convert to SAM record
		reg_to_sam(reg, read_id, sequence, output,
		           -1,    // segment_idx: -1 for single-end
		           false, // mate_mapped
		           false, // mate_rev
		           -1,    // mate_rid
		           0,     // mate_pos
		           0      // tlen
		);
	}

	// Free results
	for (int j = 0; j < n_regs; j++) {
		free(regs[j].p);
	}
	free(regs);
}

void Minimap2Aligner::align_paired(const std::string &read_id, const std::string &sequence1,
                                   const std::string &sequence2, SAMRecordBatch &output) {
	// Setup for paired-end
	int qlens[2] = {static_cast<int>(sequence1.length()), static_cast<int>(sequence2.length())};
	const char *seqs[2] = {sequence1.c_str(), sequence2.c_str()};
	int n_regs[2] = {0, 0};
	mm_reg1_t *regs[2] = {nullptr, nullptr};

	// Enable fragment mode for paired-end
	mm_mapopt_t mopt_copy = *mopt_;
	mopt_copy.flag |= MM_F_FRAG_MODE;

	mm_map_frag(index_.get(), 2, qlens, seqs, n_regs, regs, tbuf_.get(), &mopt_copy, read_id.c_str());

	// Find primary alignments for each segment (with bounds checking)
	mm_reg1_t *primary[2] = {nullptr, nullptr};
	for (int seg = 0; seg < 2; seg++) {
		for (int j = 0; j < n_regs[seg]; j++) {
			mm_reg1_t *reg = &regs[seg][j];
			// Bounds check for reference ID
			if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= subject_names_.size()) {
				continue;
			}
			if (reg->parent == reg->id) { // Primary
				primary[seg] = reg;
				break;
			}
		}
	}

	// Determine mate information
	bool mate_mapped[2] = {primary[1] != nullptr, primary[0] != nullptr};
	bool mate_rev[2] = {primary[1] ? (primary[1]->rev != 0) : false, primary[0] ? (primary[0]->rev != 0) : false};
	int32_t mate_rid[2] = {primary[1] ? primary[1]->rid : -1, primary[0] ? primary[0]->rid : -1};
	int32_t mate_pos[2] = {primary[1] ? (primary[1]->rs + 1) : 0, // 1-based
	                       primary[0] ? (primary[0]->rs + 1) : 0};

	// Calculate template length (only for proper pairs on same reference)
	int32_t tlen = 0;
	if (primary[0] && primary[1] && primary[0]->rid == primary[1]->rid) {
		// tlen = rightmost end - leftmost start, with sign indicating direction
		int32_t start0 = primary[0]->rs;
		int32_t end0 = primary[0]->re;
		int32_t start1 = primary[1]->rs;
		int32_t end1 = primary[1]->re;

		int32_t leftmost = std::min(start0, start1);
		int32_t rightmost = std::max(end0, end1);
		tlen = rightmost - leftmost;

		// Sign: positive for read1 if it's leftmost
		if (start0 > start1) {
			tlen = -tlen;
		}
	}

	// Process alignments for each segment
	for (int seg = 0; seg < 2; seg++) {
		const std::string &query_seq = (seg == 0) ? sequence1 : sequence2;
		int n_output = 0;

		for (int j = 0; j < n_regs[seg]; j++) {
			mm_reg1_t *reg = &regs[seg][j];

			// Bounds check for reference ID
			if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= subject_names_.size()) {
				continue;
			}

			bool is_primary = (reg->parent == reg->id);
			bool is_secondary = !is_primary;

			if (is_secondary) {
				int secondary_count = n_output - 1;
				if (secondary_count >= config_.max_secondary) {
					continue;
				}
			}

			// For paired, tlen has opposite sign for read2
			int32_t this_tlen = (seg == 0) ? tlen : -tlen;

			reg_to_sam(reg, read_id, query_seq, output, seg, mate_mapped[seg], mate_rev[seg], mate_rid[seg],
			           mate_pos[seg], this_tlen);

			n_output++;
		}
	}

	// Free results
	for (int seg = 0; seg < 2; seg++) {
		for (int j = 0; j < n_regs[seg]; j++) {
			free(regs[seg][j].p);
		}
		free(regs[seg]);
	}
}

void Minimap2Aligner::reg_to_sam(const void *reg_ptr, const std::string &read_id, const std::string &query_seq,
                                 SAMRecordBatch &batch, int segment_idx, bool mate_mapped, bool mate_rev,
                                 int32_t mate_rid, int32_t mate_pos, int32_t tlen) {
	const mm_reg1_t *reg = static_cast<const mm_reg1_t *>(reg_ptr);

	bool is_paired = (segment_idx >= 0);
	bool is_unmapped = (reg->rid < 0);

	batch.read_ids.push_back(read_id);
	batch.flags.push_back(calculate_flags(reg, is_paired, segment_idx, mate_mapped, mate_rev, is_unmapped));

	if (is_unmapped) {
		batch.references.push_back("*");
		batch.positions.push_back(0);
		batch.stop_positions.push_back(0);
		batch.mapqs.push_back(0);
		batch.cigars.push_back("*");
	} else {
		batch.references.push_back(get_reference_name(reg->rid));
		batch.positions.push_back(reg->rs + 1); // Convert to 1-based
		batch.stop_positions.push_back(calculate_stop_position(reg->rs + 1, reg));
		batch.mapqs.push_back(static_cast<uint8_t>(reg->mapq));
		batch.cigars.push_back(cigar_string(reg));
	}

	// Mate reference
	if (is_paired && mate_mapped && mate_rid >= 0) {
		const std::string &mate_ref = get_reference_name(mate_rid);
		if (!is_unmapped && mate_ref == batch.references.back()) {
			batch.mate_references.push_back("=");
		} else {
			batch.mate_references.push_back(mate_ref);
		}
		batch.mate_positions.push_back(mate_pos);
	} else {
		batch.mate_references.push_back("*");
		batch.mate_positions.push_back(0);
	}

	batch.template_lengths.push_back(tlen);

	// Compute alignment statistics from CIGAR
	AlignmentStats stats = compute_alignment_stats(reg);

	// Tags
	batch.tag_as_values.push_back(reg->score);
	batch.tag_xs_values.push_back(reg->subsc > 0 ? reg->subsc : -1);
	batch.tag_ys_values.push_back(-1); // Not available from minimap2
	batch.tag_xn_values.push_back(-1); // Not available from minimap2
	batch.tag_xm_values.push_back(is_unmapped ? -1 : stats.mismatches);  // XM: mismatches
	batch.tag_xo_values.push_back(is_unmapped ? -1 : stats.gap_opens);   // XO: gap opens
	batch.tag_xg_values.push_back(is_unmapped ? -1 : stats.gap_extends); // XG: gap extensions
	batch.tag_nm_values.push_back(is_unmapped ? -1 : stats.edit_distance); // NM: edit distance

	// YT tag (pair type)
	std::string yt;
	if (!is_paired) {
		yt = "UU"; // Unpaired
	} else if (mate_mapped && !is_unmapped && reg->proper_frag) {
		yt = "CP"; // Concordant pair
	} else if (mate_mapped && !is_unmapped) {
		yt = "DP"; // Discordant pair
	} else {
		yt = "UP"; // Unpaired (one unmapped)
	}
	batch.tag_yt_values.push_back(yt);

	// MD tag - generate if available
	std::string md_tag;
	if (reg->p && !is_unmapped) {
		char *md_buf = nullptr;
		int md_max_len = 0;
		int md_len = mm_gen_MD(nullptr, &md_buf, &md_max_len, index_.get(), reg, query_seq.c_str());
		if (md_len > 0 && md_buf) {
			md_tag = std::string(md_buf, md_len);
		}
		free(md_buf);
	}
	batch.tag_md_values.push_back(md_tag);

	// SA tag (supplementary alignments) - not implemented yet
	batch.tag_sa_values.push_back("");
}

std::string Minimap2Aligner::cigar_string(const void *reg_ptr) const {
	const mm_reg1_t *reg = static_cast<const mm_reg1_t *>(reg_ptr);

	if (!reg->p || reg->p->n_cigar == 0) {
		return "*";
	}

	std::ostringstream oss;
	for (uint32_t i = 0; i < reg->p->n_cigar; i++) {
		uint32_t op = reg->p->cigar[i] & 0xf;
		uint32_t len = reg->p->cigar[i] >> 4;
		oss << len << MM_CIGAR_STR[op];
	}

	return oss.str();
}

int64_t Minimap2Aligner::calculate_stop_position(int64_t start_pos, const void *reg_ptr) const {
	const mm_reg1_t *reg = static_cast<const mm_reg1_t *>(reg_ptr);

	// minimap2 uses 0-based half-open coordinates: [rs, re)
	// For example, alignment from position 0 to 49 is represented as rs=0, re=50
	// SAM format uses 1-based inclusive coordinates: [POS, stop_position]
	// Converting: if [0, 50) in 0-based half-open = [1, 50] in 1-based inclusive
	// So 0-based exclusive end (re) equals 1-based inclusive end directly
	return reg->re;
}

uint16_t Minimap2Aligner::calculate_flags(const void *reg_ptr, bool is_paired, int segment_idx, bool mate_mapped,
                                          bool mate_rev, bool is_unmapped) const {
	const mm_reg1_t *reg = static_cast<const mm_reg1_t *>(reg_ptr);

	uint16_t flags = 0;

	if (is_paired) {
		flags |= 0x1; // Paired

		if (segment_idx == 0) {
			flags |= 0x40; // First in pair
		} else {
			flags |= 0x80; // Second in pair
		}

		if (reg->proper_frag && mate_mapped && !is_unmapped) {
			flags |= 0x2; // Proper pair
		}

		if (!mate_mapped) {
			flags |= 0x8; // Mate unmapped
		}

		if (mate_rev) {
			flags |= 0x20; // Mate reverse strand
		}
	}

	if (is_unmapped) {
		flags |= 0x4; // Unmapped
	} else {
		if (reg->rev) {
			flags |= 0x10; // Reverse strand
		}
	}

	// Check for secondary/supplementary
	if (reg->parent != reg->id) {
		if (reg->sam_pri == 0) {
			flags |= 0x100; // Secondary alignment
		}
	}

	// Supplementary alignment flag
	if (reg->split == 2) { // split type indicates supplementary
		flags |= 0x800;    // Supplementary
	}

	return flags;
}

const std::string &Minimap2Aligner::get_reference_name(int32_t rid) const {
	if (rid < 0 || static_cast<size_t>(rid) >= subject_names_.size()) {
		static const std::string unknown = "*";
		return unknown;
	}
	return subject_names_[rid];
}

} // namespace miint
