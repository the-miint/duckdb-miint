#include "progressive_pcoa_core.hpp"

#include <algorithm>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "ordination.h" // skbb_pcoa_fsvd_fp32

#include "procrustes_core.hpp"   // FitProcrustes / ApplyToReference / ApplyToOther
#include "unifrac_omp_scope.hpp" // OmpThreadScope — process-wide libssu/skbb serialization

namespace miint::progressive {

namespace {

using ::miint::procrustes::ApplyToOther;
using ::miint::procrustes::ApplyToReference;
using ::miint::procrustes::Disparity;
using ::miint::procrustes::FitProcrustes;
using ::miint::procrustes::ProcrustesFit;

// PCoA on a dense m×m fp32 distance block → m×d row-major double coordinates
// (sample-major, one row per block id in block.ids order). When `eig_out` /
// `prop_out` are non-null they receive the d eigenvalues / proportions.
//
// Delegates to libskbb's randomized FSVD — the exact call and sample-major
// layout that unifrac_pcoa_function.cpp's RunPcoaOnMatrix uses (the ordination.h
// header comment says "(n_eighs × n_dims)" but the implementation writes
// sample-major n×n_eighs). The skbb call is wrapped in the process-wide
// OmpThreadScope exactly like every other libssu/skbb caller in the codebase:
// the scope both pins the OpenMP thread count and serializes against skbb calls
// from OTHER concurrent DuckDB connections (which would otherwise race on the
// global omp thread-count and libssu RNG state). This core running one batch at a
// time does not make the mutex unnecessary — a second, unrelated query can be in
// skbb simultaneously. (A future inner-parallel-batches path would need a
// mutex-free skbb entry point + the vendored report_status thread_local patch
// before this scope could be dropped.)
std::vector<double> PcoaBlock(const DistanceBlock &block, uint32_t d, int seed, int n_threads,
                              std::vector<double> *eig_out, std::vector<double> *prop_out) {
	const uint32_t m = static_cast<uint32_t>(block.ids.size());
	std::vector<float> eig(d), samples(static_cast<size_t>(m) * d), prop(d);
	{
		miint::unifrac::OmpThreadScope omp_scope(n_threads);
		skbb_pcoa_fsvd_fp32(m, block.matrix.data(), d, seed, eig.data(), samples.data(), prop.data());
	}
	std::vector<double> coords(static_cast<size_t>(m) * d);
	for (size_t i = 0; i < coords.size(); ++i) {
		coords[i] = static_cast<double>(samples[i]);
	}
	if (eig_out) {
		eig_out->assign(eig.begin(), eig.end());
	}
	if (prop_out) {
		prop_out->assign(prop.begin(), prop.end());
	}
	return coords;
}

// Verify a provider returned a block over EXACTLY the requested samples: same
// count, every requested id present, none extra or duplicated. This is the core's
// fail-loud guard against a BlockProvider silently dropping a sample — e.g. a
// UniFrac sub-biom that omits a sample with no surviving features, or a distance
// slice missing a row. Without it, a dropped batch sample would simply never be
// emitted (short result, no error); a dropped anchor would mis-key the fit.
void ValidateBlockMatchesRequest(const DistanceBlock &block, const std::vector<std::string> &requested,
                                 const char *which) {
	if (block.ids.size() != requested.size()) {
		throw std::invalid_argument("progressive_pcoa: block provider returned " + std::to_string(block.ids.size()) +
		                            " samples for a " + std::to_string(requested.size()) + "-sample " + which +
		                            " request");
	}
	std::unordered_set<std::string> req(requested.begin(), requested.end());
	std::unordered_set<std::string> seen;
	seen.reserve(block.ids.size());
	for (const auto &id : block.ids) {
		if (!req.count(id)) {
			throw std::invalid_argument("progressive_pcoa: block provider returned sample '" + id + "' not in the " +
			                            which + " request");
		}
		if (!seen.insert(id).second) {
			throw std::invalid_argument("progressive_pcoa: block provider returned duplicate sample '" + id +
			                            "' in the " + which + " block");
		}
	}
}

} // namespace

ProgressivePcoaResult RunProgressivePcoa(const std::vector<std::string> &anchor_ids,
                                         const std::vector<std::string> &remaining_ids, uint32_t n_dims,
                                         uint32_t batch_size, int seed, int n_threads, const BlockProvider &get_block,
                                         const WavePrefetch &prefetch, uint32_t wave_batches) {
	if (n_dims < 1) {
		throw std::invalid_argument("progressive_pcoa: n_dims must be >= 1");
	}
	if (batch_size == 0) {
		throw std::invalid_argument("progressive_pcoa: batch_size must be >= 1");
	}
	if (n_threads < 1) {
		throw std::invalid_argument("progressive_pcoa: n_threads must be >= 1");
	}
	const uint32_t d = n_dims;
	const uint32_t a = static_cast<uint32_t>(anchor_ids.size());
	// Need a >= d + 1 for both the reference PCoA (loses one dim to centering) and
	// the procrustes fit (d + 1 points determine a d-dimensional transform). Phrase
	// the guard as `a <= d` so it stays correct even if d is near UINT32_MAX (d + 1
	// would wrap): a > d is exactly a >= d + 1 for integers.
	if (a <= d) {
		throw std::invalid_argument("progressive_pcoa: need at least n_dims + 1 anchors to define the reference frame");
	}

	// Anchors must be distinct and disjoint from the streamed samples: the fit
	// pairs anchor rows 1:1 across the reference and each batch, so a duplicate or
	// an anchor that also appears as a batch sample would silently corrupt the
	// alignment (and double-emit that sample).
	std::unordered_set<std::string> anchor_set;
	anchor_set.reserve(a);
	for (const auto &id : anchor_ids) {
		if (!anchor_set.insert(id).second) {
			throw std::invalid_argument("progressive_pcoa: duplicate anchor id '" + id + "'");
		}
	}
	for (const auto &rid : remaining_ids) {
		if (anchor_set.count(rid)) {
			throw std::invalid_argument("progressive_pcoa: sample '" + rid + "' is both an anchor and a batch sample");
		}
	}

	// ── Phase 0: reference PCoA on the anchors — defines the common frame ────────
	DistanceBlock ref_block = get_block(anchor_ids);
	ValidateBlockMatchesRequest(ref_block, anchor_ids, "anchor");

	ProgressivePcoaResult result;
	result.d = d;
	// Reference coordinates in ref_block.ids order. This is also where the reported
	// eigenvalues / proportions come from (the anchor ordination — a documented
	// caveat: they describe the anchors, not the full data).
	const std::vector<double> ref_coords =
	    PcoaBlock(ref_block, d, seed, n_threads, &result.eigvals, &result.proportion_explained);

	// Emit the reference anchor coordinates once, standardized into the shared
	// frame. ApplyToReference uses only the fit's reference translate/norm, which
	// depend solely on ref_coords, so a self-fit produces the identical
	// standardization every batch fit would — and works even with zero batches.
	// Every batch below aligns its samples into this same standardized frame.
	const ProcrustesFit ref_fit = FitProcrustes(ref_coords.data(), ref_coords.data(), a, d);
	{
		std::vector<double> ref_std(static_cast<size_t>(a) * d);
		ApplyToReference(ref_fit, ref_coords.data(), a, ref_std.data());
		result.coords.reserve(static_cast<size_t>(a + remaining_ids.size()) * d);
		for (uint32_t r = 0; r < a; ++r) {
			for (uint32_t axis = 0; axis < d; ++axis) {
				result.coords.push_back({ref_block.ids[r], static_cast<int32_t>(axis), ref_std[r * d + axis]});
			}
		}
	}

	// ── Phase 1: stream the remaining samples in contiguous batches ──────────────
	// Batches are grouped into waves: before fetching a wave's blocks we announce
	// all of its requests, so a provider can serve them from one pass over its
	// source instead of one pass per block (see WavePrefetch). The batch bodies
	// below are unchanged by this — a wave is purely an I/O grouping, and W has no
	// effect on the coordinates produced.
	const size_t total = remaining_ids.size();

	// Block request = this batch's non-anchor samples followed by all anchors. Built
	// in one place so the request announced by prefetch is byte-identical to the one
	// get_block is then called with — a provider is entitled to key its cache on it.
	const auto build_request = [&](size_t start, size_t end) {
		std::vector<std::string> req;
		req.reserve((end - start) + a);
		req.insert(req.end(), remaining_ids.begin() + start, remaining_ids.begin() + end);
		req.insert(req.end(), anchor_ids.begin(), anchor_ids.end());
		return req;
	};

	// One batch: fetch its block, ordinate it, fit it onto the reference frame
	// through the anchor overlap, and emit its non-anchor coordinates.
	const auto run_one_batch = [&](size_t start, size_t end) {
		const std::vector<std::string> req = build_request(start, end);
		const DistanceBlock blk = get_block(req);
		ValidateBlockMatchesRequest(blk, req, "batch");
		const uint32_t m = static_cast<uint32_t>(blk.ids.size());
		const std::vector<double> batch_coords = PcoaBlock(blk, d, seed, n_threads, nullptr, nullptr);

		// id → batch row index.
		std::unordered_map<std::string, uint32_t> batch_row;
		batch_row.reserve(m);
		for (uint32_t r = 0; r < m; ++r) {
			batch_row.emplace(blk.ids[r], r);
		}

		// Anchor rows of this batch (raw), ordered to pair 1:1 with ref_coords by
		// ref_block.ids order — so batch_anchor row r corresponds to reference row r.
		std::vector<double> batch_anchor(static_cast<size_t>(a) * d);
		for (uint32_t r = 0; r < a; ++r) {
			auto it = batch_row.find(ref_block.ids[r]);
			if (it == batch_row.end()) {
				throw std::invalid_argument("progressive_pcoa: anchor '" + ref_block.ids[r] +
				                            "' missing from batch block");
			}
			const double *src = batch_coords.data() + static_cast<size_t>(it->second) * d;
			std::copy(src, src + d, batch_anchor.begin() + static_cast<size_t>(r) * d);
		}

		// Fit the batch onto the reference through the anchor overlap.
		const ProcrustesFit fit = FitProcrustes(ref_coords.data(), batch_anchor.data(), a, d);

		// Score that fit and report it. This is the run's own accuracy evidence
		// (see BatchDiagnostic): the anchors are the only samples this batch and
		// the reference frame have in common, so the disparity over them is what
		// says whether the batch's placement can be trusted. Costs O(a·d) against
		// a block PCoA, i.e. nothing.
		{
			std::vector<double> ref_std(static_cast<size_t>(a) * d);
			std::vector<double> anchor_fit(static_cast<size_t>(a) * d);
			ApplyToReference(fit, ref_coords.data(), a, ref_std.data());
			ApplyToOther(fit, batch_anchor.data(), a, anchor_fit.data());
			BatchDiagnostic diag;
			diag.batch = static_cast<int32_t>(result.batches.size());
			diag.n_samples = static_cast<uint32_t>(end - start);
			diag.anchor_m2 = Disparity(ref_std.data(), anchor_fit.data(), a, d);
			result.batches.push_back(diag);
		}

		// This batch's non-anchor rows (raw), then mapped into the reference frame.
		std::vector<std::string> nb_ids;
		std::vector<double> nb_raw;
		nb_ids.reserve(end - start);
		nb_raw.reserve((end - start) * static_cast<size_t>(d));
		for (uint32_t r = 0; r < m; ++r) {
			if (anchor_set.count(blk.ids[r])) {
				continue; // anchors already emitted via the reference frame
			}
			nb_ids.push_back(blk.ids[r]);
			const double *src = batch_coords.data() + static_cast<size_t>(r) * d;
			nb_raw.insert(nb_raw.end(), src, src + d);
		}
		const uint32_t nb = static_cast<uint32_t>(nb_ids.size());
		std::vector<double> nb_fit(static_cast<size_t>(nb) * d);
		ApplyToOther(fit, nb_raw.data(), nb, nb_fit.data());
		const int32_t batch_index = result.batches.back().batch;
		for (uint32_t r = 0; r < nb; ++r) {
			for (uint32_t axis = 0; axis < d; ++axis) {
				result.coords.push_back({nb_ids[r], static_cast<int32_t>(axis), nb_fit[r * d + axis], batch_index});
			}
		}
	};

	std::vector<std::pair<size_t, size_t>> batch_ranges;
	for (size_t start = 0; start < total; start += batch_size) {
		batch_ranges.emplace_back(start, std::min(total, start + static_cast<size_t>(batch_size)));
	}
	const size_t wave_width = std::max<size_t>(1, wave_batches);

	for (size_t wave_start = 0; wave_start < batch_ranges.size(); wave_start += wave_width) {
		const size_t wave_end = std::min(batch_ranges.size(), wave_start + wave_width);
		// Announced for EVERY wave, including a trailing one-batch wave: a uniform
		// contract lets a provider serve all batch blocks from its wave cache and
		// keep a direct path only for the one-off anchor (reference) block. The
		// alternative — skipping the announcement when there is nothing to amortize
		// — would force every provider to implement both paths for batch blocks.
		if (prefetch) {
			std::vector<std::vector<std::string>> wave_requests;
			wave_requests.reserve(wave_end - wave_start);
			for (size_t b = wave_start; b < wave_end; ++b) {
				wave_requests.push_back(build_request(batch_ranges[b].first, batch_ranges[b].second));
			}
			prefetch(wave_requests);
		}
		for (size_t b = wave_start; b < wave_end; ++b) {
			run_one_batch(batch_ranges[b].first, batch_ranges[b].second);
		}
	}

	return result;
}

} // namespace miint::progressive
