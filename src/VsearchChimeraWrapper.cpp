// VsearchChimeraWrapper: bridges vsearch's chimera detection library API
// to the miint DuckDB extension's table function interface.
//
// vsearch (BSD-2-Clause) provides the UCHIME algorithm with SIMD-optimized
// Needleman-Wunsch alignment, DUST masking, k-mer candidate search, and
// parent selection. This wrapper manages vsearch's global state lifecycle
// and converts results to miint's UchimeResult format.

#include "VsearchChimeraWrapper.hpp"

// vsearch headers (Pimpl pattern — only included in this .cpp file).
// vsearch.h is the aggregate header: config.h, system headers, global opt_*
// variables, and core module headers (db.h, searchcore.h, etc.).
// Additional module headers included explicitly for their API declarations.
#include "vsearch.h"
#include "chimera.h"
#include "dbindex.h"
#include "mask.h"

#include <cstring>

namespace miint {

struct VsearchChimeraWrapper::Impl {
	chimera_info_s *ci = nullptr;
	bool ci_initialized = false;

	Impl() : ci(chimera_info_alloc()) {
	}
	~Impl() {
		if (ci) {
			chimera_info_free(ci);
		}
	}
	Impl(const Impl &) = delete;
	Impl &operator=(const Impl &) = delete;
};

VsearchChimeraWrapper::VsearchChimeraWrapper(const UchimeParams &params)
    : impl_(std::make_unique<Impl>()), params_(params) {
}

VsearchChimeraWrapper::~VsearchChimeraWrapper() {
	if (impl_ && impl_->ci_initialized) {
		chimera_detect_cleanup(impl_->ci);
	}
}

VsearchChimeraWrapper::VsearchChimeraWrapper(VsearchChimeraWrapper &&) noexcept = default;
VsearchChimeraWrapper &VsearchChimeraWrapper::operator=(VsearchChimeraWrapper &&) noexcept = default;

void VsearchChimeraWrapper::init_vsearch_globals() {
	// Set vsearch global options to match our UchimeParams.
	// These are the opt_* extern variables defined in vsearch.cc.
	opt_minh = params_.minh;
	opt_xn = params_.xn;
	opt_dn = params_.dn;
	opt_mindiv = params_.mindiv;
	opt_mindiffs = params_.mindiffs;
	opt_abskew = params_.abskew;

	// Alignment scoring defaults (matching vsearch's chimera detection)
	opt_match = 2;
	opt_mismatch = -4;
	opt_gap_open_query_interior = 20;
	opt_gap_open_target_interior = 20;
	opt_gap_open_query_left = 20;
	opt_gap_open_target_left = 20;
	opt_gap_open_query_right = 20;
	opt_gap_open_target_right = 20;
	opt_gap_extension_query_interior = 2;
	opt_gap_extension_target_interior = 2;
	opt_gap_extension_query_left = 1;
	opt_gap_extension_target_left = 1;
	opt_gap_extension_query_right = 1;
	opt_gap_extension_target_right = 1;

	// Masking and search defaults
	opt_dbmask = MASK_DUST;
	opt_qmask = MASK_DUST;
	opt_hardmask = 0;
	opt_wordlength = 8;
	opt_notrunclabels = 0;
	opt_xsize = 0;
	opt_xee = 0;
	opt_xlength = 0;
	opt_fasta_score = 0;

	// Chimera-specific (NULL = not in that mode)
	opt_uchime_ref = nullptr;
	opt_uchime_denovo = nullptr;
	opt_uchime2_denovo = nullptr;
	opt_uchime3_denovo = nullptr;
	opt_chimeras_denovo = nullptr;
	opt_chimeras = nullptr;
	opt_nonchimeras = nullptr;
	opt_borderline = nullptr;
	opt_uchimealns = nullptr;
	opt_uchimeout = nullptr;
	opt_uchimeout5 = 0;
	opt_tabbedout = nullptr;
	opt_alnout = nullptr;
}

void VsearchChimeraWrapper::set_reference(const std::vector<std::string> &labels,
                                          const std::vector<std::string> &sequences) {
	ref_labels_ = labels;
	ref_sequences_ = sequences;
	ref_abundances_.assign(sequences.size(), 0);

	init_vsearch_globals();

	// Clear any previous database
	db_free();
	dbindex_free();

	// Load sequences into vsearch's global DB via db_add()
	for (size_t i = 0; i < labels.size(); i++) {
		// Normalize U→T for RNA sequences
		std::string seq = sequences[i];
		for (auto &c : seq) {
			if (c == 'U') {
				c = 'T';
			} else if (c == 'u') {
				c = 't';
			}
		}

		db_add(false, // not fastq
		       labels[i].c_str(), seq.c_str(),
		       nullptr, // no quality
		       labels[i].size(), seq.size(),
		       1); // abundance (default 1 for uchime_ref)
	}

	// Apply DUST masking to all sequences
	dust_all();

	// Build k-mer index
	dbindex_prepare(1, opt_dbmask);
	dbindex_addallsequences(opt_dbmask);

	// Initialize per-thread chimera working state
	if (impl_->ci_initialized) {
		chimera_detect_cleanup(impl_->ci);
		chimera_info_free(impl_->ci);
		impl_->ci = chimera_info_alloc();
	}
	chimera_detect_init(impl_->ci);
	impl_->ci_initialized = true;

	initialized_ = true;
}

void VsearchChimeraWrapper::add_to_reference(const std::string &label, const std::string &sequence, int64_t abundance) {
	ref_labels_.push_back(label);
	ref_sequences_.push_back(sequence);
	ref_abundances_.push_back(abundance);

	std::string seq = sequence;
	for (auto &c : seq) {
		if (c == 'U') {
			c = 'T';
		} else if (c == 'u') {
			c = 't';
		}
	}

	db_add(false, label.c_str(), seq.c_str(), nullptr, label.size(), seq.size(), abundance);
}

void VsearchChimeraWrapper::index_sequence(uint64_t seqno) {
	dbindex_addsequence(static_cast<unsigned int>(seqno), opt_qmask);
}

UchimeResult VsearchChimeraWrapper::convert_result(const void *vsearch_result_ptr) {
	const auto *r = static_cast<const chimera_result_s *>(vsearch_result_ptr);

	UchimeResult result;
	result.score = r->score;
	result.query_label = r->query_label;
	result.flag = std::string(1, r->flag);

	if (r->flag == 'N') {
		// Non-chimeric: NULL parents and identities (vsearch convention)
		return result;
	}

	result.parent_a_label = r->parent_a_label;
	result.parent_b_label = r->parent_b_label;
	result.closest_parent_label = r->closest_parent_label;
	result.id_query_model = r->id_query_model;
	result.id_query_a = r->id_query_a;
	result.id_query_b = r->id_query_b;
	result.id_a_b = r->id_a_b;
	result.id_query_top = r->id_query_top;
	result.left_yes = r->left_yes;
	result.left_no = r->left_no;
	result.left_abstain = r->left_abstain;
	result.right_yes = r->right_yes;
	result.right_no = r->right_no;
	result.right_abstain = r->right_abstain;
	result.divergence = r->divergence;

	return result;
}

UchimeResult VsearchChimeraWrapper::detect(const std::string &query_label, const std::string &query_sequence) {
	// Normalize U→T
	std::string seq = query_sequence;
	for (auto &c : seq) {
		if (c == 'U') {
			c = 'T';
		} else if (c == 'u') {
			c = 't';
		}
	}

	chimera_result_s vsearch_result {};
	chimera_detect_single(impl_->ci, seq.c_str(), query_label.c_str(), static_cast<int>(seq.size()),
	                      1, // abundance = 1 for uchime_ref
	                      &vsearch_result);

	return convert_result(&vsearch_result);
}

UchimeResult VsearchChimeraWrapper::detect_denovo(const std::string &query_label, const std::string &query_sequence,
                                                  int64_t query_abundance) {
	std::string seq = query_sequence;
	for (auto &c : seq) {
		if (c == 'U') {
			c = 'T';
		} else if (c == 'u') {
			c = 't';
		}
	}

	chimera_result_s vsearch_result {};
	chimera_detect_single(impl_->ci, seq.c_str(), query_label.c_str(), static_cast<int>(seq.size()),
	                      static_cast<int>(query_abundance), &vsearch_result);

	return convert_result(&vsearch_result);
}

} // namespace miint
