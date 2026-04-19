#include "SortMeRNAAligner.hpp"

#include "smr_api.h"

#include <cstring>
#include <memory>
#include <stdexcept>
#include <vector>

namespace miint {

namespace {

void apply_config(const SortMeRNAConfig &in, smr_config_t &out) {
	smr_config_init(&out);
	if (in.num_threads > 0)
		out.num_threads = in.num_threads;
	out.match = in.match;
	out.mismatch = in.mismatch;
	out.gap_open = in.gap_open;
	out.gap_ext = in.gap_ext;
	out.score_N = in.score_N;
	out.evalue = in.evalue;
	out.seed_win_len = in.seed_win_len;
	out.num_alignments = in.num_alignments;
	out.best = in.best ? 1 : 0;
	out.paired = in.paired ? 1 : 0;
	out.forward_only = in.forward_only ? 1 : 0;
	out.reverse_only = in.reverse_only ? 1 : 0;
	out.full_search = in.full_search ? 1 : 0;
}

struct SmrOutputDeleter {
	void operator()(smr_output_t *p) const noexcept {
		if (p)
			smr_output_free(p);
	}
};

using SmrOutputPtr = std::unique_ptr<smr_output_t, SmrOutputDeleter>;

} // namespace

SortMeRNAAligner::SortMeRNAAligner(const SortMeRNAConfig &cfg, const std::vector<std::string> &ref_paths)
    : paired_(cfg.paired) {
	if (ref_paths.empty()) {
		throw std::invalid_argument("SortMeRNAAligner: ref_paths is empty");
	}

	smr_config_t c_cfg;
	apply_config(cfg, c_cfg);

	ctx_ = smr_ctx_create(&c_cfg);
	if (!ctx_) {
		throw std::runtime_error("SortMeRNAAligner: smr_ctx_create returned NULL");
	}

	std::vector<const char *> ref_ptrs;
	ref_ptrs.reserve(ref_paths.size());
	for (const auto &p : ref_paths) {
		ref_ptrs.push_back(p.c_str());
	}

	idx_ = smr_index_load(ctx_, ref_ptrs.data(), static_cast<int32_t>(ref_ptrs.size()));
	if (!idx_) {
		int rc = smr_last_error_code(ctx_);
		const char *last = smr_last_error(ctx_);
		std::string msg = "smr_index_load: ";
		msg += smr_strerror(rc ? rc : SMR_ERR_IO);
		if (last && *last) {
			msg += ": ";
			msg += last;
		}
		smr_ctx_destroy(ctx_);
		ctx_ = nullptr;
		throw std::runtime_error(msg);
	}
}

SortMeRNAAligner::~SortMeRNAAligner() {
	if (idx_)
		smr_index_free(idx_);
	if (ctx_)
		smr_ctx_destroy(ctx_);
}

void SortMeRNAAligner::align(const SortMeRNAQueryBatch &queries, SortMeRNAResultBatch &out) {
	if (queries.empty()) {
		throw std::invalid_argument("SortMeRNAAligner::align: empty query batch");
	}
	if (queries.read_ids.size() != queries.sequences.size()) {
		throw std::invalid_argument("SortMeRNAAligner::align: read_ids.size() != sequences.size()");
	}
	if (!queries.sequences2.empty() && queries.sequences2.size() != queries.sequences.size()) {
		throw std::invalid_argument("SortMeRNAAligner::align: sequences2 size does not match sequences");
	}
	if (queries.is_paired() != paired_) {
		throw std::invalid_argument("SortMeRNAAligner::align: query batch paired-ness does not match aligner config");
	}

	const size_t n_in = queries.size();
	const int32_t segments = paired_ ? 2 : 1;
	const size_t n_seqs = n_in * static_cast<size_t>(segments);

	// Synthetic IDs "0".."n_seqs-1" avoid the upstream uniqueness contract.
	// Keep the std::strings alive for the duration of the smr_run call.
	std::vector<std::string> syn_ids;
	syn_ids.reserve(n_seqs);
	for (size_t i = 0; i < n_seqs; ++i) {
		syn_ids.emplace_back(std::to_string(i));
	}

	std::vector<smr_seq_t> smr_seqs(n_seqs);
	if (paired_) {
		for (size_t i = 0; i < n_in; ++i) {
			smr_seqs[2 * i].id = syn_ids[2 * i].c_str();
			smr_seqs[2 * i].sequence = queries.sequences[i].c_str();
			smr_seqs[2 * i].quality = nullptr;
			smr_seqs[2 * i + 1].id = syn_ids[2 * i + 1].c_str();
			smr_seqs[2 * i + 1].sequence = queries.sequences2[i].c_str();
			smr_seqs[2 * i + 1].quality = nullptr;
		}
	} else {
		for (size_t i = 0; i < n_in; ++i) {
			smr_seqs[i].id = syn_ids[i].c_str();
			smr_seqs[i].sequence = queries.sequences[i].c_str();
			smr_seqs[i].quality = nullptr;
		}
	}

	smr_output_t *raw = nullptr;
	smr_stats_t stats = {};
	int rc = smr_run_seqs_with_index(idx_, smr_seqs.data(), static_cast<int32_t>(n_seqs), &raw, &stats);
	if (rc != SMR_OK) {
		const char *last = smr_last_error(ctx_);
		std::string msg = "smr_run_seqs_with_index: ";
		msg += smr_strerror(rc);
		if (last && *last) {
			msg += ": ";
			msg += last;
		}
		throw std::runtime_error(msg);
	}
	// Ownership of raw transfers to us on SMR_OK; guard it immediately so any
	// exception from the reserve/push_back loop below still frees it.
	SmrOutputPtr raw_out(raw);

	// The library contract is one output row per input seq (aligned=0 for
	// reads without a hit). If that ever changes upstream, the positional
	// ID check below would index past syn_ids — fail loudly instead.
	if (raw_out->num_reads != static_cast<uint64_t>(n_seqs)) {
		throw std::runtime_error("SortMeRNAAligner::align: library returned " + std::to_string(raw_out->num_reads) +
		                         " rows for " + std::to_string(n_seqs) + " input seqs (expected 1:1)");
	}

	out.reserve(out.size() + raw_out->num_reads);
	for (size_t i = 0; i < raw_out->num_reads; ++i) {
		// Verify the library returned rows in the same order we submitted them.
		// We rely on this positional mapping for paired-end input_row math and
		// for mapping output back to caller-supplied read_ids. The upstream
		// traversal guarantee isn't in the public API contract, so fail loud if
		// it ever changes.
		const char *out_id = raw_out->read_ids[i];
		if (!out_id || syn_ids[i] != out_id) {
			throw std::runtime_error("SortMeRNAAligner::align: library returned rows out of order "
			                         "(row " +
			                         std::to_string(i) + " id='" + (out_id ? out_id : "(null)") + "' expected '" +
			                         syn_ids[i] + "')");
		}

		const size_t input_row = paired_ ? (i / 2) : i;
		const int32_t segment_idx = paired_ ? static_cast<int32_t>(i % 2) : 0;

		out.read_ids.push_back(queries.read_ids[input_row]);
		out.aligned.push_back(raw_out->aligned[i]);
		out.strands.push_back(raw_out->strand[i]);
		out.ref_names.emplace_back(raw_out->ref_name[i] ? raw_out->ref_name[i] : "");
		out.ref_starts.push_back(raw_out->ref_start[i]);
		out.ref_ends.push_back(raw_out->ref_end[i]);
		out.cigars.emplace_back(raw_out->cigar[i] ? raw_out->cigar[i] : "");
		out.scores.push_back(raw_out->score[i]);
		out.e_values.push_back(raw_out->e_value[i]);
		out.identities.push_back(raw_out->identity[i]);
		out.coverages.push_back(raw_out->coverage[i]);
		out.edit_distances.push_back(raw_out->edit_distance[i]);
		out.segment_indices.push_back(segment_idx);
	}
}

void SortMeRNAResultBatch::clear() {
	read_ids.clear();
	aligned.clear();
	strands.clear();
	ref_names.clear();
	ref_starts.clear();
	ref_ends.clear();
	cigars.clear();
	scores.clear();
	e_values.clear();
	identities.clear();
	coverages.clear();
	edit_distances.clear();
	segment_indices.clear();
}

void SortMeRNAResultBatch::reserve(size_t n) {
	read_ids.reserve(n);
	aligned.reserve(n);
	strands.reserve(n);
	ref_names.reserve(n);
	ref_starts.reserve(n);
	ref_ends.reserve(n);
	cigars.reserve(n);
	scores.reserve(n);
	e_values.reserve(n);
	identities.reserve(n);
	coverages.reserve(n);
	edit_distances.reserve(n);
	segment_indices.reserve(n);
}

std::string SortMeRNAAligner::version() {
	const char *v = smr_version();
	return v ? std::string(v) : std::string("unknown");
}

} // namespace miint
