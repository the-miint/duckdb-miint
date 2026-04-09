// VsearchSearchWrapper: bridges vsearch's global search library API
// to the miint DuckDB extension's table function interface.
//
// vsearch (BSD-2-Clause) provides usearch_global with SIMD-optimized
// Needleman-Wunsch alignment, DUST masking, and k-mer candidate search.
// This wrapper manages vsearch's global state lifecycle and converts
// results to miint's SearchResult format.
//
// Session lifecycle (from vsearch_api.h):
//   1.  vsearch_init_defaults()            — acquires session mutex
//   2.  set opt_* overrides (id, maxaccepts, maxrejects)
//   3.  vsearch_apply_defaults_fixups()
//   4.  db_init / db_add / dust_all        — load reference DB
//   5.  dbindex_prepare / addallsequences   — build k-mer index
//   6.  search_info_alloc + search_init     — per-thread search state (N times)
//   7.  search_single(si, ...)              — per-query search (thread-safe)
//   8.  search_cleanup + search_info_free   — per-thread teardown (N times)
//   9.  dbindex_free / db_free              — release database
//   10. vsearch_session_end()               — release session mutex
//
// Session serialization: vsearch uses ~200 global opt_* variables.
// vsearch_init_defaults() acquires a process-wide session mutex that
// blocks until vsearch_session_end() releases it. This means only one
// vsearch-backed operation (search, chimera, cluster) can be active at
// a time across all DuckDB queries. A second search_sequences_vsearch() call
// will block until the first completes.

#include "VsearchSearchWrapper.hpp"
#include "vsearch_utils.hpp"

#include "vsearch_api.h"

namespace miint {

// ============================================================================
// SearchHandle — per-thread search state
// ============================================================================

struct VsearchSearchWrapper::SearchHandle::Impl {
	searchinfo_s *si = nullptr;
	bool initialized = false;
	int maxaccepts = 1;
	std::vector<search_result_s> result_buf;

	Impl(int max_accepts) : maxaccepts(max_accepts), result_buf(max_accepts) {
	}

	~Impl() {
		if (initialized) {
			search_cleanup(si);
		}
		if (si) {
			search_info_free(si);
		}
	}
	Impl(const Impl &) = delete;
	Impl &operator=(const Impl &) = delete;
};

VsearchSearchWrapper::SearchHandle::SearchHandle(int maxaccepts) : impl_(std::make_unique<Impl>(maxaccepts)) {
}

VsearchSearchWrapper::SearchHandle::~SearchHandle() = default;
VsearchSearchWrapper::SearchHandle::SearchHandle(SearchHandle &&) noexcept = default;
VsearchSearchWrapper::SearchHandle &VsearchSearchWrapper::SearchHandle::operator=(SearchHandle &&) noexcept = default;

std::vector<SearchResult> VsearchSearchWrapper::SearchHandle::search(const std::string &query_label,
                                                                     const std::string &query_sequence,
                                                                     int query_size) {
	if (!impl_ || !impl_->initialized) {
		throw std::runtime_error("SearchHandle::search called on uninitialized handle");
	}

	validate_label_length(query_label, "search");

	if (query_sequence.empty()) {
		return {};
	}

	auto seq = normalize_rna(query_sequence);

	int result_count = 0;
	search_single(impl_->si, seq.c_str(), query_label.c_str(), static_cast<int>(seq.size()), query_size,
	              impl_->result_buf.data(), impl_->maxaccepts, &result_count);

	std::vector<SearchResult> results;
	results.reserve(result_count);
	for (int i = 0; i < result_count; i++) {
		auto &r = impl_->result_buf[i];
		SearchResult sr;
		sr.query_id = query_label;
		sr.target_id = r.target_label;
		sr.identity = r.id;
		sr.matches = r.matches;
		sr.mismatches = r.mismatches;
		sr.gaps = r.gaps;
		sr.alignment_length = r.alignment_length;
		sr.query_length = r.query_length;
		sr.target_length = r.target_length;
		sr.accepted = r.accepted;
		results.push_back(std::move(sr));
	}

	return results;
}

// ============================================================================
// VsearchState — all vsearch global state in one place for correct teardown
// ============================================================================

struct VsearchSearchWrapper::VsearchState {
	bool db_initialized = false;
	bool index_initialized = false;
	bool session_mutex_held = false;

	VsearchState() = default;

	~VsearchState() {
		// Teardown in reverse order of initialization.
		// index_initialized is separate from db_initialized because
		// dbindex_free() on stale globals (from a previous session that
		// was already freed) causes double-free. Only call dbindex_free()
		// if dbindex_prepare() completed in THIS session.
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
// VsearchSearchWrapper — public API
// ============================================================================

VsearchSearchWrapper::VsearchSearchWrapper(const SearchParams &params) : params_(params) {
}

VsearchSearchWrapper::~VsearchSearchWrapper() = default;

VsearchSearchWrapper::VsearchSearchWrapper(VsearchSearchWrapper &&other) noexcept
    : params_(other.params_), state_(std::move(other.state_)) {
}

VsearchSearchWrapper &VsearchSearchWrapper::operator=(VsearchSearchWrapper &&other) noexcept {
	if (this != &other) {
		state_.reset();
		params_ = other.params_;
		state_ = std::move(other.state_);
	}
	return *this;
}

void VsearchSearchWrapper::set_database(const std::vector<std::string> &labels,
                                        const std::vector<std::string> &sequences) {
	if (labels.empty()) {
		throw std::runtime_error("Reference database is empty");
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
	// opt_wordlength MUST be set before fixups — fixups resolve
	// opt_minwordmatches using minwordmatches_defaults[opt_wordlength].
	opt_wordlength = 8;
	opt_id = params_.id;
	opt_maxaccepts = params_.maxaccepts;
	opt_maxrejects = params_.maxrejects;
	vsearch_apply_defaults_fixups();

	// Step 4: load DB
	db_init();
	state->db_initialized = true;

	for (size_t i = 0; i < labels.size(); i++) {
		validate_label_length(labels[i], "set_database");
		auto seq = normalize_rna(sequences[i]);
		db_add(false, labels[i].c_str(), seq.c_str(), nullptr, labels[i].size(), seq.size(), 1);
	}

	// Step 5: DUST masking + k-mer index
	dust_all();
	dbindex_prepare(1, opt_dbmask);
	dbindex_addallsequences(opt_dbmask);
	state->index_initialized = true;

	state_ = std::move(state);
}

VsearchSearchWrapper::SearchHandle VsearchSearchWrapper::create_search_handle() {
	if (!state_ || !state_->db_initialized) {
		throw std::runtime_error("create_search_handle called before set_database");
	}

	SearchHandle handle(params_.maxaccepts);
	handle.impl_->si = search_info_alloc();
	search_init(handle.impl_->si);
	handle.impl_->initialized = true;
	return handle;
}

void VsearchSearchWrapper::search_batch(const std::vector<std::string> &query_labels,
                                        const std::vector<std::string> &query_sequences,
                                        std::vector<SearchResult> &output) {
	if (!state_ || !state_->db_initialized) {
		throw std::runtime_error("search_batch called before set_database");
	}
	if (query_labels.empty()) {
		return;
	}

	VsearchBatchArgs args(query_labels, query_sequences);
	int max_results = params_.maxaccepts;

	std::vector<search_result_s> raw_results(args.count * max_results);
	std::vector<int> result_counts(args.count);

	::search_batch(args.seq_ptrs.data(), args.head_ptrs.data(), args.lens.data(), args.sizes.data(), args.count,
	               raw_results.data(), max_results, result_counts.data());

	for (int qi = 0; qi < args.count; qi++) {
		for (int hi = 0; hi < result_counts[qi]; hi++) {
			auto &r = raw_results[qi * max_results + hi];
			SearchResult sr;
			sr.query_id = query_labels[qi];
			sr.target_id = r.target_label;
			sr.identity = r.id;
			sr.matches = r.matches;
			sr.mismatches = r.mismatches;
			sr.gaps = r.gaps;
			sr.alignment_length = r.alignment_length;
			sr.query_length = r.query_length;
			sr.target_length = r.target_length;
			sr.accepted = r.accepted;
			output.push_back(std::move(sr));
		}
	}
}

} // namespace miint
