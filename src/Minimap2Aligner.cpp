#include "Minimap2Aligner.hpp"
#include <minimap2/minimap.h>
#include <minimap2/mmpriv.h>
#include <cstdlib>
#include <stdexcept>

// When MIINT_USE_JEMALLOC is defined, minimap2 is compiled with malloc/free
// redirected to duckdb_je_* (jemalloc). Our C++ code must free minimap2-allocated
// memory through the same allocator. On builds without jemalloc (musl, macOS),
// minimap2 uses system malloc and we use system free.
#ifdef MIINT_USE_JEMALLOC
#include "duckdb_je_decl.h"
#define MM_FREE(ptr) duckdb_je_free(ptr)
#else
#define MM_FREE(ptr) free(ptr)
#endif

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

// SharedMinimap2Index implementation
SharedMinimap2Index::SharedMinimap2Index(const std::string &index_path, const Minimap2Config &config) : mopt_() {
	mm_idxopt_t iopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt_);
	mm_idx_t *raw_idx = nullptr;
	Minimap2Aligner::LoadIndexFromFile(index_path, iopt, raw_idx, subject_names_);
	index_.reset(raw_idx);
	mm_mapopt_update(&mopt_, index_.get());
}

SharedMinimap2Index::SharedMinimap2Index(mm_idx_t *idx, const mm_mapopt_t &mopt, std::vector<std::string> subject_names)
    : index_(idx), mopt_(mopt), subject_names_(std::move(subject_names)) {
}

SharedMinimap2Index::~SharedMinimap2Index() = default;

const mm_idx_t *SharedMinimap2Index::index() const {
	return index_.get();
}

const mm_mapopt_t &SharedMinimap2Index::mapopt() const {
	return mopt_;
}

const std::vector<std::string> &SharedMinimap2Index::subject_names() const {
	return subject_names_;
}

// Static helper: initialize indexing and mapping options from config
void Minimap2Aligner::InitOptions(const Minimap2Config &config, mm_idxopt_t &iopt, mm_mapopt_t &mopt) {
	// Initialize with default values first (required before mm_set_opt)
	mm_idxopt_init(&iopt);
	mm_mapopt_init(&mopt);

	// Apply preset-specific options
	if (mm_set_opt(config.preset.c_str(), &iopt, &mopt) != 0) {
		throw std::runtime_error("Unknown minimap2 preset: " + config.preset);
	}

	// Validate preset set valid k and w values
	if (iopt.k <= 0 || iopt.k > 28) {
		throw std::runtime_error("Preset '" + config.preset + "' set invalid k-mer size: " + std::to_string(iopt.k));
	}
	if (iopt.w <= 0 || iopt.w >= 256) {
		throw std::runtime_error("Preset '" + config.preset + "' set invalid window size: " + std::to_string(iopt.w));
	}

	// Override k and w if specified
	if (config.k > 0) {
		if (config.k > 28) {
			throw std::runtime_error("k-mer size must be <= 28 (got " + std::to_string(config.k) + ")");
		}
		iopt.k = static_cast<short>(config.k);
	}
	if (config.w > 0) {
		if (config.w >= 256) {
			throw std::runtime_error("Window size must be < 256 (got " + std::to_string(config.w) + ")");
		}
		iopt.w = static_cast<short>(config.w);
	}

	// Enable CIGAR output
	mopt.flag |= MM_F_CIGAR;

	// Enable EQX mode (=/X instead of M)
	if (config.eqx) {
		mopt.flag |= MM_F_EQX;
	}

	// Enable MD tag output
	mopt.flag |= MM_F_OUT_MD;

	// Set max secondary alignments
	// best_n controls how many chains go to DP alignment and how many
	// secondary alignments are retained. Match minimap2 command-line
	// behavior: -N sets best_n directly (default 5).
	mopt.best_n = config.max_secondary;

	// Coverage pre-filter threshold
	mopt.min_chain_coverage = config.min_chain_coverage;

	// Map-time scoring/chaining overrides. Applied after mm_set_opt(preset) so a caller
	// can tune alignment against a prebuilt index. mm_mapopt_update() (run after the index
	// loads) only adjusts occurrence thresholds, not these fields, so the overrides persist.
	// -1 / -1.0f means "unset" (keep the preset default); ParseMinimap2ConfigParams has
	// already rejected negative user values, so >= 0 reliably distinguishes a set override.
	if (config.match_score >= 0) {
		mopt.a = config.match_score;
	}
	if (config.mismatch_penalty >= 0) {
		mopt.b = config.mismatch_penalty;
	}
	if (config.gap_open >= 0) {
		mopt.q = config.gap_open;
	}
	if (config.gap_extend >= 0) {
		mopt.e = config.gap_extend;
	}
	if (config.gap_open2 >= 0) {
		mopt.q2 = config.gap_open2;
	}
	if (config.gap_extend2 >= 0) {
		mopt.e2 = config.gap_extend2;
	}
	if (config.bandwidth >= 0) {
		mopt.bw = config.bandwidth;
	}
	if (config.zdrop >= 0) {
		mopt.zdrop = config.zdrop;
	}
	if (config.zdrop_inv >= 0) {
		mopt.zdrop_inv = config.zdrop_inv;
	}
	if (config.min_chain_score >= 0) {
		mopt.min_chain_score = config.min_chain_score;
	}
	if (config.min_count >= 0) {
		mopt.min_cnt = config.min_count;
	}
	if (config.max_gap >= 0) {
		mopt.max_gap = config.max_gap;
	}
	// min_dp_max is the minimal peak DP score to keep a hit (minimap2 -s). The preset sets it
	// (sr -> 40); like minimap2's CLI we do NOT auto-recompute it from match_score/min_chain_score,
	// so a caller that drops the match score should set min_dp_max too if the default would filter
	// their hits.
	if (config.min_dp_max >= 0) {
		mopt.min_dp_max = config.min_dp_max;
	}
	if (config.pri_ratio >= 0.0f) {
		mopt.pri_ratio = config.pri_ratio;
	}
	if (config.mask_level >= 0.0f) {
		mopt.mask_level = config.mask_level;
	}
}

// Static helper: load index from .mmi file
void Minimap2Aligner::LoadIndexFromFile(const std::string &path, const mm_idxopt_t &iopt, mm_idx_t *&out_idx,
                                        std::vector<std::string> &out_names) {
	mm_idx_reader_t *reader = mm_idx_reader_open(path.c_str(), &iopt, nullptr);
	if (!reader) {
		throw std::runtime_error("Cannot open index file: " + path);
	}

	mm_idx_t *idx = mm_idx_reader_read(reader, 1);
	mm_idx_reader_close(reader);

	if (!idx) {
		throw std::runtime_error("Failed to load index from: " + path);
	}

	// Extract reference names from loaded index
	out_names.clear();
	out_names.reserve(idx->n_seq);
	for (uint32_t i = 0; i < idx->n_seq; i++) {
		if (!idx->seq[i].name) {
			mm_idx_destroy(idx);
			throw std::runtime_error("Index contains unnamed sequence at position " + std::to_string(i) +
			                         " in file: " + path);
		}
		out_names.push_back(std::string(idx->seq[i].name));
	}

	out_idx = idx;
}

// Constructor
Minimap2Aligner::Minimap2Aligner(const Minimap2Config &config)
    : config_(config), iopt_(std::make_unique<mm_idxopt_t>()), mopt_(std::make_unique<mm_mapopt_t>()),
      tbuf_(mm_tbuf_init()) {
	InitOptions(config_, *iopt_, *mopt_);
}

// Destructor
Minimap2Aligner::~Minimap2Aligner() {
	MM_FREE(md_buf_);
}

// Move constructor
Minimap2Aligner::Minimap2Aligner(Minimap2Aligner &&other) noexcept
    : config_(std::move(other.config_)), iopt_(std::move(other.iopt_)), mopt_(std::move(other.mopt_)),
      index_(std::move(other.index_)), subject_names_(std::move(other.subject_names_)), tbuf_(std::move(other.tbuf_)),
      shared_index_(std::move(other.shared_index_)), md_buf_(other.md_buf_), md_max_len_(other.md_max_len_) {
	other.md_buf_ = nullptr;
	other.md_max_len_ = 0;
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
		shared_index_ = std::move(other.shared_index_);
		MM_FREE(md_buf_);
		md_buf_ = other.md_buf_;
		md_max_len_ = other.md_max_len_;
		other.md_buf_ = nullptr;
		other.md_max_len_ = 0;
	}
	return *this;
}

void Minimap2Aligner::attach_shared_index(std::shared_ptr<SharedMinimap2Index> shared_idx) {
	// Clear owned index (mutually exclusive)
	index_.reset();
	subject_names_.clear();
	shared_index_ = std::move(shared_idx);
}

void Minimap2Aligner::detach_shared_index() {
	shared_index_.reset();
}

mm_idx_t *Minimap2Aligner::BuildRawIndex(const std::vector<AlignmentSubject> &subjects, const mm_idxopt_t &iopt,
                                         std::vector<std::string> &out_names) {
	if (subjects.empty()) {
		throw std::runtime_error("Cannot build index from empty subject list");
	}

	std::vector<const char *> seqs;
	std::vector<const char *> names;

	seqs.reserve(subjects.size());
	names.reserve(subjects.size());
	out_names.clear();
	out_names.reserve(subjects.size());

	for (const auto &subject : subjects) {
		if (subject.sequence.empty()) {
			throw std::runtime_error("Cannot build index: sequence '" + subject.read_id + "' is empty");
		}
		seqs.push_back(subject.sequence.c_str());
		names.push_back(subject.read_id.c_str());
		out_names.push_back(subject.read_id);
	}

	mm_idx_t *idx = mm_idx_str(iopt.w, iopt.k, iopt.flag & 1, iopt.bucket_bits, static_cast<int>(subjects.size()),
	                           seqs.data(), names.data());

	if (!idx) {
		throw std::runtime_error("Failed to build minimap2 index");
	}

	return idx;
}

std::shared_ptr<SharedMinimap2Index> Minimap2Aligner::BuildSharedIndex(const std::vector<AlignmentSubject> &subjects,
                                                                       const Minimap2Config &config) {
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	InitOptions(config, iopt, mopt);

	std::vector<std::string> subject_names;
	mm_idx_t *idx = BuildRawIndex(subjects, iopt, subject_names);
	mm_mapopt_update(&mopt, idx);

	return std::make_shared<SharedMinimap2Index>(idx, mopt, std::move(subject_names));
}

const mm_idx_t *Minimap2Aligner::active_index() const {
	if (shared_index_) {
		return shared_index_->index();
	}
	return index_.get();
}

const mm_mapopt_t *Minimap2Aligner::active_mapopt() const {
	if (shared_index_) {
		return &shared_index_->mapopt();
	}
	return mopt_.get();
}

const std::vector<std::string> &Minimap2Aligner::active_subject_names() const {
	if (shared_index_) {
		return shared_index_->subject_names();
	}
	return subject_names_;
}

void Minimap2Aligner::build_index(const std::vector<AlignmentSubject> &subjects) {
	index_.reset(BuildRawIndex(subjects, *iopt_, subject_names_));
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

	if (!active_index()) {
		throw std::runtime_error("No index built. Call build_index() or attach_shared_index() first.");
	}

	// Heuristic: ~3 alignments per query (primary + secondaries)
	output.reserve(output.size() + queries.size() * 3);

	// Process each query
	for (size_t i = 0; i < queries.size(); i++) {
		// Check if this query is actually paired (has non-empty sequence2)
		// is_paired flag just means the column exists, not that all rows have data
		bool actually_paired = queries.is_paired && i < queries.sequences2.size() && !queries.sequences2[i].empty();

		if (actually_paired) {
			align_paired(queries.read_ids[i], queries.sequences1[i], queries.sequences2[i], output);
		} else {
			align_single(queries.read_ids[i], queries.sequences1[i], output);
		}
	}
}

void Minimap2Aligner::align_single(const std::string &read_id, const std::string &sequence, SAMRecordBatch &output) {
	// Skip empty query sequences (minimap2 requires len > 0)
	if (sequence.empty()) {
		return; // No alignments for empty query
	}

	int n_regs = 0;
	mm_reg1_t *regs = mm_map(active_index(), static_cast<int>(sequence.length()), sequence.c_str(), &n_regs,
	                         tbuf_.get(), active_mapopt(), read_id.c_str());

	int secondary_count = 0;

	// Process alignments (minimap2 returns them sorted by score, best first)
	for (int j = 0; j < n_regs; j++) {
		mm_reg1_t *reg = &regs[j];

		// Bounds check for reference ID before using it
		if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= active_subject_names().size()) {
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

	// Free results (must use MM_FREE — minimap2 may use jemalloc allocator)
	for (int j = 0; j < n_regs; j++) {
		MM_FREE(regs[j].p);
	}
	MM_FREE(regs);
}

void Minimap2Aligner::align_paired(const std::string &read_id, const std::string &sequence1,
                                   const std::string &sequence2, SAMRecordBatch &output) {
	// Skip if both sequences are empty (minimap2 requires len > 0)
	if (sequence1.empty() && sequence2.empty()) {
		return; // No alignments for empty queries
	}

	// Setup for paired-end
	int qlens[2] = {static_cast<int>(sequence1.length()), static_cast<int>(sequence2.length())};
	const char *seqs[2] = {sequence1.c_str(), sequence2.c_str()};
	int n_regs[2] = {0, 0};
	mm_reg1_t *regs[2] = {nullptr, nullptr};

	// Enable fragment mode for paired-end
	mm_mapopt_t mopt_copy = *active_mapopt();
	mopt_copy.flag |= MM_F_FRAG_MODE;

	mm_map_frag(active_index(), 2, qlens, seqs, n_regs, regs, tbuf_.get(), &mopt_copy, read_id.c_str());

	// Find primary alignments for each segment (with bounds checking)
	mm_reg1_t *primary[2] = {nullptr, nullptr};
	for (int seg = 0; seg < 2; seg++) {
		for (int j = 0; j < n_regs[seg]; j++) {
			mm_reg1_t *reg = &regs[seg][j];
			// Bounds check for reference ID
			if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= active_subject_names().size()) {
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
		int secondary_count = 0;

		for (int j = 0; j < n_regs[seg]; j++) {
			mm_reg1_t *reg = &regs[seg][j];

			// Bounds check for reference ID
			if (reg->rid < 0 || static_cast<size_t>(reg->rid) >= active_subject_names().size()) {
				continue;
			}

			bool is_primary = (reg->parent == reg->id);
			bool is_secondary = !is_primary;

			if (is_secondary) {
				if (secondary_count >= config_.max_secondary) {
					continue;
				}
				secondary_count++;
			}

			// For paired, tlen has opposite sign for read2
			int32_t this_tlen = (seg == 0) ? tlen : -tlen;

			reg_to_sam(reg, read_id, query_seq, output, seg, mate_mapped[seg], mate_rev[seg], mate_rid[seg],
			           mate_pos[seg], this_tlen);
		}
	}

	// Free results (must use MM_FREE — minimap2 may use jemalloc allocator)
	for (int seg = 0; seg < 2; seg++) {
		for (int j = 0; j < n_regs[seg]; j++) {
			MM_FREE(regs[seg][j].p);
		}
		MM_FREE(regs[seg]);
	}
}

void Minimap2Aligner::reg_to_sam(const mm_reg1_t *reg, const std::string &read_id, const std::string &query_seq,
                                 SAMRecordBatch &batch, int segment_idx, bool mate_mapped, bool mate_rev,
                                 int32_t mate_rid, int32_t mate_pos, int32_t tlen) {
	bool is_paired = (segment_idx >= 0);
	bool is_unmapped = (reg->rid < 0);

	uint16_t flags = calculate_flags(reg, is_paired, segment_idx, mate_mapped, mate_rev, is_unmapped);
	batch.read_ids.push_back(read_id);
	batch.flags.push_back(flags);

	// Compute CIGAR string and alignment stats in a single pass (Opt 2+3)
	AlignmentStats stats;
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
		batch.cigars.push_back(cigar_string(reg, static_cast<int32_t>(query_seq.length()), flags, &stats));
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

	// Tags
	batch.tag_as_values.push_back(is_unmapped ? -1 : (reg->p ? reg->p->dp_score : -1));
	batch.tag_xs_values.push_back(reg->subsc > 0 ? reg->subsc : -1);
	batch.tag_ys_values.push_back(-1);                                   // Not available from minimap2
	batch.tag_xn_values.push_back(-1);                                   // Not available from minimap2
	batch.tag_xm_values.push_back(is_unmapped ? -1 : stats.mismatches);  // XM: mismatches
	batch.tag_xo_values.push_back(is_unmapped ? -1 : stats.gap_opens);   // XO: gap opens
	batch.tag_xg_values.push_back(is_unmapped ? -1 : stats.gap_extends); // XG: gap extensions

	// NM: use minimap2's O(1) formula instead of CIGAR walk (Opt 7)
	int64_t nm = reg->blen - reg->mlen + (reg->p ? reg->p->n_ambi : 0);
	batch.tag_nm_values.push_back(is_unmapped ? -1 : nm);

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

	// MD tag - use arena allocator for internal temp buffers (Opt 1),
	// and pooled md_buf_ member to avoid per-call alloc/free (Opt 5)
	std::string md_tag;
	if (reg->p && !is_unmapped) {
		int md_len =
		    mm_gen_MD(mm_tbuf_get_km(tbuf_.get()), &md_buf_, &md_max_len_, active_index(), reg, query_seq.c_str());
		if (md_len > 0 && md_buf_) {
			md_tag = std::string(md_buf_, md_len);
		}
		// Do NOT free md_buf_ here — reused across calls, freed in destructor
	}
	batch.tag_md_values.push_back(md_tag);

	// SA tag (supplementary alignments) - not implemented yet
	batch.tag_sa_values.push_back("");
}

std::string Minimap2Aligner::cigar_string(const mm_reg1_t *reg, int32_t query_len, uint16_t sam_flags,
                                          AlignmentStats *stats_out) const {
	if (!reg->p || reg->p->n_cigar == 0) {
		return "*";
	}

	// Bounds validation: qs/qe must be within [0, query_len]
	if (reg->qs < 0 || reg->qe < 0 || reg->qs > query_len || reg->qe > query_len || reg->qs > reg->qe) {
		return "*";
	}

	// Compute clip lengths from query start/end, accounting for strand
	// (matches minimap2 format.c:490-491)
	int clip_front = reg->rev ? (query_len - reg->qe) : reg->qs;
	int clip_back = reg->rev ? reg->qs : (query_len - reg->qe);

	// Supplementary (0x800) → hard clip; primary/secondary → soft clip.
	char clip_char = (sam_flags & 0x800) ? 'H' : 'S';

	// Pre-size: each CIGAR op is at most 10 digits + 1 char, plus 2 clips
	std::string result;
	result.reserve((reg->p->n_cigar + 2) * 11);

	auto append_int = [&result](uint32_t val) {
		char buf[10];
		int len = 0;
		do {
			buf[len++] = '0' + (val % 10);
			val /= 10;
		} while (val);
		for (int i = len - 1; i >= 0; --i) {
			result.push_back(buf[i]);
		}
	};

	if (clip_front > 0) {
		append_int(clip_front);
		result.push_back(clip_char);
	}

	for (uint32_t i = 0; i < reg->p->n_cigar; i++) {
		uint32_t op = reg->p->cigar[i] & 0xf;
		uint32_t len = reg->p->cigar[i] >> 4;
		append_int(len);
		result.push_back(MM_CIGAR_STR[op]);
		if (stats_out) {
			switch (op) {
			case MM_CIGAR_X_MISMATCH:
				stats_out->mismatches += len;
				break;
			case MM_CIGAR_INS:
			case MM_CIGAR_DEL:
				stats_out->gap_opens += 1;
				stats_out->gap_extends += (len > 1) ? (len - 1) : 0;
				break;
			default:
				break;
			}
		}
	}

	if (clip_back > 0) {
		append_int(clip_back);
		result.push_back(clip_char);
	}

	return result;
}

int64_t Minimap2Aligner::calculate_stop_position(int64_t start_pos, const mm_reg1_t *reg) const {
	// minimap2 uses 0-based half-open coordinates: [rs, re)
	// stop_position uses 1-based half-open coordinates: position + ref_length
	// Convert: 0-based exclusive end → 1-based exclusive end = re + 1
	// Example: rs=0, re=50 → position=1, stop_position=51, length=50
	return reg->re + 1;
}

uint16_t Minimap2Aligner::calculate_flags(const mm_reg1_t *reg, bool is_paired, int segment_idx, bool mate_mapped,
                                          bool mate_rev, bool is_unmapped) const {
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

	// Secondary/supplementary detection (matches minimap2 format.c:546-547)
	if (reg->parent != reg->id) {
		flags |= 0x100; // Secondary alignment
	} else if (!reg->sam_pri) {
		flags |= 0x800; // Supplementary alignment
	}

	return flags;
}

const std::string &Minimap2Aligner::get_reference_name(int32_t rid) const {
	auto &names = active_subject_names();
	if (rid < 0 || static_cast<size_t>(rid) >= names.size()) {
		static const std::string unknown = "*";
		return unknown;
	}
	return names[rid];
}

void Minimap2Aligner::load_index(const std::string &index_path) {
	mm_idx_t *idx = nullptr;
	LoadIndexFromFile(index_path, *iopt_, idx, subject_names_);

	// Store index and update mapping options
	index_.reset(idx);
	mm_mapopt_update(mopt_.get(), index_.get());

	// Clear shared index (owned and shared are mutually exclusive)
	shared_index_.reset();
}

void Minimap2Aligner::save_index(const std::string &output_path) const {
	// Validate index exists
	if (!index_) {
		throw std::runtime_error("No index to save. Build or load an index first.");
	}

	// Open output file
	FILE *fp = fopen(output_path.c_str(), "wb");
	if (!fp) {
		throw std::runtime_error("Cannot create index file: " + output_path);
	}

	// Write index using minimap2 API
	mm_idx_dump(fp, index_.get());

	// Check for write errors after dump
	if (ferror(fp)) {
		fclose(fp);
		throw std::runtime_error("Write error while saving index to: " + output_path);
	}

	// Close and check for errors (fclose flushes, so check return value)
	if (fclose(fp) != 0) {
		throw std::runtime_error("Error closing index file: " + output_path);
	}
}

bool Minimap2Aligner::is_index_file(const std::string &path) {
	// Use minimap2's built-in check function
	// Returns file size if valid index, 0 if not
	int64_t file_size = mm_idx_is_idx(path.c_str());
	return file_size > 0;
}

} // namespace miint
