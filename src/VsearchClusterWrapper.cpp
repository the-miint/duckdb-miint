// VsearchClusterWrapper: bridges vsearch's greedy clustering library API
// to the miint DuckDB extension.
//
// vsearch (BSD-2-Clause) provides incremental greedy clustering with
// SIMD-optimized Needleman-Wunsch alignment and k-mer candidate search.
// This wrapper manages vsearch's global state lifecycle and converts
// results to miint's ClusterResult format.
//
// Session lifecycle:
//   1.  vsearch_init_defaults()            — acquires session mutex
//   2.  set opt_* overrides (id, strand, wordlength)
//   3.  vsearch_apply_defaults_fixups()
//   4.  db_init / db_add / dust_all        — load sequences
//   5.  dbindex_prepare(1, mask)           — allocate index (empty)
//   6.  cluster_session_alloc + init       — allocate clustering state
//   7.  cluster_assign_single(cs, i, &r)   — sequential assignment
//   8.  cluster_session_cleanup + free     — teardown clustering
//   9.  dbindex_free / db_free             — release database
//   10. vsearch_session_end()              — release session mutex
//
// Centroids are indexed incrementally by vsearch during step 7 — do NOT
// call dbindex_addallsequences().

#include "VsearchClusterWrapper.hpp"
#include "vsearch_utils.hpp"

#include "vsearch_api.h"

namespace miint {

// ============================================================================
// VsearchState — all vsearch global state in one place for correct teardown
// ============================================================================

struct VsearchClusterWrapper::VsearchState {
	bool db_initialized = false;
	bool index_initialized = false;
	bool session_mutex_held = false;
	int sequence_count = 0;

	// Stored labels for mapping seqno -> read_id in results.
	// vsearch's db_getheader() returns the label, but we store our own copy
	// to avoid depending on vsearch's internal string lifecycle after db_free.
	std::vector<std::string> labels;

	VsearchState() = default;

	~VsearchState() {
		if (index_initialized) {
			dbindex_free();
		}
		if (db_initialized) {
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
// VsearchClusterWrapper — public API
// ============================================================================

VsearchClusterWrapper::VsearchClusterWrapper(const ClusterParams &params) : params_(params) {
}

VsearchClusterWrapper::~VsearchClusterWrapper() = default;

VsearchClusterWrapper::VsearchClusterWrapper(VsearchClusterWrapper &&other) noexcept
    : params_(other.params_), state_(std::move(other.state_)) {
}

VsearchClusterWrapper &VsearchClusterWrapper::operator=(VsearchClusterWrapper &&other) noexcept {
	if (this != &other) {
		state_.reset();
		params_ = other.params_;
		state_ = std::move(other.state_);
	}
	return *this;
}

void VsearchClusterWrapper::set_sequences(const std::vector<std::string> &labels,
                                          const std::vector<std::string> &sequences) {
	if (labels.empty()) {
		throw std::runtime_error("Input sequence set is empty");
	}
	if (labels.size() != sequences.size()) {
		throw std::runtime_error("Labels and sequences must have the same size");
	}

	// Release previous session if re-initializing.
	state_.reset();

	auto state = std::make_unique<VsearchState>();

	// Step 1: acquire session mutex
	vsearch_init_defaults();
	state->session_mutex_held = true;

	// Step 2-3: set parameters and resolve sentinels.
	opt_wordlength = 8;
	opt_id = params_.id;
	opt_strand = params_.strand;
	opt_maxaccepts = 1;
	opt_maxrejects = 32;
	// opt_threads must be set before fixups: fixups resolve 0 → arch_get_cores().
	if (params_.threads > 0) {
		opt_threads = params_.threads;
	}
	vsearch_apply_defaults_fixups();

	// Step 4: load DB
	db_init();
	state->db_initialized = true;

	state->labels.reserve(labels.size());
	for (size_t i = 0; i < labels.size(); i++) {
		validate_label_length(labels[i], "set_sequences");
		auto seq = normalize_rna(sequences[i]);
		db_add(false, labels[i].c_str(), seq.c_str(), nullptr, labels[i].size(), seq.size(), 1);
		state->labels.push_back(labels[i]);
	}

	// DUST masking
	dust_all();

	// Step 5: prepare k-mer index (empty — centroids indexed incrementally)
	dbindex_prepare(1, opt_dbmask);
	state->index_initialized = true;

	state->sequence_count = static_cast<int>(labels.size());
	state_ = std::move(state);
}

std::vector<ClusterResult> VsearchClusterWrapper::cluster_all() {
	if (!state_ || !state_->db_initialized) {
		throw std::runtime_error("cluster_all called before set_sequences");
	}

	// Step 6: allocate clustering session
	struct cluster_session_s *cs = cluster_session_alloc();
	cluster_session_init(cs);

	// Step 7: assign sequences to clusters using batch API
	// (internally parallelizes search across opt_threads)
	std::vector<cluster_result_s> raw_results(state_->sequence_count);
	cluster_assign_batch(cs, 0, state_->sequence_count, raw_results.data());

	std::vector<ClusterResult> results;
	results.reserve(state_->sequence_count);
	for (int i = 0; i < state_->sequence_count; i++) {
		auto &raw = raw_results[i];
		ClusterResult cr;
		cr.read_id = state_->labels[i];
		cr.is_centroid = raw.is_centroid;
		cr.cluster_id = raw.cluster_id;
		cr.centroid_id = raw.centroid_label;
		cr.identity = raw.identity;
		cr.cigar = raw.cigar;
		cr.cigar_truncated = raw.cigar_truncated;
		results.push_back(std::move(cr));
	}

	// Step 8: cleanup clustering session
	cluster_session_cleanup(cs);
	cluster_session_free(cs);

	return results;
}

} // namespace miint
