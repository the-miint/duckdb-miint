// VsearchChimeraWrapper: bridges vsearch's chimera detection library API
// to the miint DuckDB extension's table function interface.
//
// vsearch (BSD-2-Clause) provides the UCHIME algorithm with SIMD-optimized
// Needleman-Wunsch alignment, DUST masking, k-mer candidate search, and
// parent selection. This wrapper manages vsearch's global state lifecycle
// and converts results to miint's UchimeResult format.
//
// Session lifecycle (from vsearch_api.h):
//   1.  vsearch_init_defaults()            — acquires session mutex
//   2.  set opt_* overrides
//   3.  vsearch_apply_defaults_fixups()
//   4.  db_init / db_add / dust_all        — load reference
//   5.  dbindex_prepare / addallsequences   — build k-mer index
//   6.  chimera_session_init()              — session-level chimera setup
//   7.  chimera_detect_thread_init(ci)      — per-thread working state (N times)
//   8.  chimera_detect_single(ci, ...)      — detection (thread-safe)
//   9.  chimera_detect_thread_cleanup(ci)   — per-thread teardown (N times)
//   10. chimera_session_cleanup()           — session-level chimera teardown
//   11. dbindex_free / db_free              — release database
//   12. vsearch_session_end()               — release session mutex

#include "VsearchChimeraWrapper.hpp"
#include "vsearch_utils.hpp"

#include "vsearch_api.h"

#include <cassert>

namespace miint {

// ============================================================================
// DetectHandle — per-thread chimera detection state
// ============================================================================

struct VsearchChimeraWrapper::DetectHandle::Impl {
	chimera_info_s *ci = nullptr;
	bool initialized = false;

	// ci is allocated lazily by create_detect_handle, NOT in the constructor,
	// because chimera_info_alloc must be called after vsearch_init_defaults.
	Impl() = default;

	~Impl() {
		if (initialized) {
			chimera_detect_thread_cleanup(ci);
		}
		if (ci) {
			chimera_info_free(ci);
		}
	}
	Impl(const Impl &) = delete;
	Impl &operator=(const Impl &) = delete;
};

VsearchChimeraWrapper::DetectHandle::DetectHandle() : impl_(std::make_unique<Impl>()) {
}

VsearchChimeraWrapper::DetectHandle::~DetectHandle() = default;
VsearchChimeraWrapper::DetectHandle::DetectHandle(DetectHandle &&) noexcept = default;
VsearchChimeraWrapper::DetectHandle &VsearchChimeraWrapper::DetectHandle::operator=(DetectHandle &&) noexcept = default;

UchimeResult VsearchChimeraWrapper::DetectHandle::detect(const std::string &query_label,
                                                         const std::string &query_sequence) {
	assert(impl_ && impl_->initialized && "DetectHandle::detect called on uninitialized handle");

	validate_label_length(query_label, "detect");

	auto seq = normalize_rna(query_sequence);

	chimera_result_s vsearch_result {};
	chimera_detect_single(impl_->ci, seq.c_str(), query_label.c_str(), static_cast<int>(seq.size()),
	                      1, // abundance = 1 for uchime_ref
	                      &vsearch_result);

	return convert_result(&vsearch_result);
}

// ============================================================================
// VsearchState — all vsearch global state in one place for correct teardown
// ============================================================================

struct VsearchChimeraWrapper::VsearchState {
	chimera_info_s *ci = nullptr;
	bool thread_initialized = false;
	bool session_initialized = false;
	bool db_initialized = false;
	bool session_mutex_held = false;

	VsearchState() = default;

	~VsearchState() {
		// Teardown in reverse order of initialization (steps 9-12).
		if (thread_initialized) {
			chimera_detect_thread_cleanup(ci);
		}
		if (ci) {
			chimera_info_free(ci);
		}
		if (session_initialized) {
			chimera_session_cleanup();
		}
		if (db_initialized) {
			dbindex_free();
			db_free();
		}
		if (session_mutex_held) {
			vsearch_session_end();
		}
	}

	VsearchState(const VsearchState &) = delete;
	VsearchState &operator=(const VsearchState &) = delete;
};

// ============================================================================
// VsearchChimeraWrapper — public API
// ============================================================================

VsearchChimeraWrapper::VsearchChimeraWrapper(const UchimeParams &params) : params_(params) {
	// No vsearch state allocated here. chimera_info_alloc() requires
	// vsearch_init_defaults() to have been called first, so we defer
	// all allocation to set_reference() / prepare_denovo().
}

VsearchChimeraWrapper::~VsearchChimeraWrapper() {
	// VsearchState destructor handles all cleanup in correct order.
	// If state_ is null (never initialized, or moved-from), this is a no-op.
}

VsearchChimeraWrapper::VsearchChimeraWrapper(VsearchChimeraWrapper &&other) noexcept
    : params_(other.params_), ref_labels_(std::move(other.ref_labels_)),
      ref_sequences_(std::move(other.ref_sequences_)), ref_abundances_(std::move(other.ref_abundances_)),
      state_(std::move(other.state_)) {
	// other.state_ is now nullptr — its destructor is a no-op.
}

VsearchChimeraWrapper &VsearchChimeraWrapper::operator=(VsearchChimeraWrapper &&other) noexcept {
	if (this != &other) {
		// Release our current session (if any). VsearchState destructor
		// handles all cleanup including session mutex release.
		state_.reset();

		params_ = other.params_;
		ref_labels_ = std::move(other.ref_labels_);
		ref_sequences_ = std::move(other.ref_sequences_);
		ref_abundances_ = std::move(other.ref_abundances_);
		state_ = std::move(other.state_);
	}
	return *this;
}

void VsearchChimeraWrapper::init_common(bool denovo) {
	// Allocate state RAII object. If anything below throws, the destructor
	// releases the session mutex and any partially-initialized resources.
	auto state = std::make_unique<VsearchState>();

	// Step 1: acquire session mutex
	vsearch_init_defaults();
	state->session_mutex_held = true;

	// Step 2-3: set parameters and resolve sentinels
	opt_wordlength = 8;
	opt_minh = params_.minh;
	opt_xn = params_.xn;
	opt_dn = params_.dn;
	opt_mindiv = params_.mindiv;
	opt_mindiffs = params_.mindiffs;
	// abskew: set explicitly in both modes to avoid depending on vsearch defaults.
	opt_abskew = denovo ? params_.abskew : 0.0;
	vsearch_apply_defaults_fixups();

	// Step 4-5: load DB (caller finishes index setup after this)
	db_init();
	state->db_initialized = true;

	// Commit state to the wrapper. Caller (set_reference/prepare_denovo)
	// continues filling in DB, index, and chimera init on this state.
	state_ = std::move(state);
}

void VsearchChimeraWrapper::load_sequences_into_db(const std::vector<std::string> &labels,
                                                   const std::vector<std::string> &sequences,
                                                   const std::vector<int64_t> &abundances) {
	for (size_t i = 0; i < labels.size(); i++) {
		validate_label_length(labels[i], "chimera detection");
		auto seq = normalize_rna(sequences[i]);
		int64_t abund = i < abundances.size() ? abundances[i] : 1;
		db_add(false, labels[i].c_str(), seq.c_str(), nullptr, labels[i].size(), seq.size(), abund);
	}
}

void VsearchChimeraWrapper::set_reference(const std::vector<std::string> &labels,
                                          const std::vector<std::string> &sequences) {
	ref_labels_ = labels;
	ref_sequences_ = sequences;
	ref_abundances_.assign(sequences.size(), 1);

	// Release previous session if re-initializing. VsearchState destructor
	// handles all cleanup including session mutex release.
	state_.reset();

	init_common(false); // Reference mode: abskew=0

	load_sequences_into_db(labels, sequences, ref_abundances_);
	dust_all();
	dbindex_prepare(1, opt_dbmask);
	dbindex_addallsequences(opt_dbmask);

	// Step 6-7: chimera session + per-thread init
	chimera_session_init();
	state_->session_initialized = true;
	state_->ci = chimera_info_alloc();
	chimera_detect_thread_init(state_->ci);
	state_->thread_initialized = true;
}

void VsearchChimeraWrapper::prepare_denovo(const std::vector<std::string> &labels,
                                           const std::vector<std::string> &sequences,
                                           const std::vector<int64_t> &abundances) {
	ref_labels_ = labels;
	ref_sequences_ = sequences;
	ref_abundances_ = abundances;

	state_.reset();

	init_common(true); // De novo mode: abskew from params

	load_sequences_into_db(labels, sequences, abundances);
	dust_all();
	dbindex_prepare(1, opt_dbmask);
	// De novo: do NOT index all sequences — caller indexes incrementally.

	chimera_session_init();
	state_->session_initialized = true;
	state_->ci = chimera_info_alloc();
	chimera_detect_thread_init(state_->ci);
	state_->thread_initialized = true;
}

VsearchChimeraWrapper::DetectHandle VsearchChimeraWrapper::create_detect_handle() {
	assert(state_ && state_->session_initialized && "create_detect_handle called before set_reference/prepare_denovo");

	DetectHandle handle;
	handle.impl_->ci = chimera_info_alloc();
	chimera_detect_thread_init(handle.impl_->ci);
	handle.impl_->initialized = true;
	return handle;
}

void VsearchChimeraWrapper::index_sequence(uint64_t seqno) {
	dbindex_addsequence(static_cast<unsigned int>(seqno), opt_dbmask);
}

UchimeResult VsearchChimeraWrapper::convert_result(const void *vsearch_result_ptr) {
	const auto *r = static_cast<const chimera_result_s *>(vsearch_result_ptr);

	UchimeResult result;
	result.score = r->score;
	result.query_label = r->query_label;
	result.flag = std::string(1, r->flag);

	if (r->flag == 'N') {
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

UchimeResult VsearchChimeraWrapper::detect_denovo(const std::string &query_label, const std::string &query_sequence,
                                                  int64_t query_abundance) {
	assert(state_ && state_->thread_initialized && "detect_denovo called before prepare_denovo");

	validate_label_length(query_label, "detect_denovo");

	auto seq = normalize_rna(query_sequence);

	chimera_result_s vsearch_result {};
	chimera_detect_single(state_->ci, seq.c_str(), query_label.c_str(), static_cast<int>(seq.size()),
	                      static_cast<int>(query_abundance), &vsearch_result);

	return convert_result(&vsearch_result);
}

} // namespace miint
