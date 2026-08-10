#include "progressive_pcoa_core.hpp"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <iterator>
#include <stdexcept>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "ordination.h" // skbb_pcoa_fsvd_fp32

#include "procrustes_core.hpp"   // FitProcrustes / ApplyToReference / ApplyToOther
#include "unifrac_omp_scope.hpp" // ComputeCallScope — per-call OpenMP pin + non-negative seed

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
// sample-major n×n_eighs).
//
// Each call takes a ComputeCallScope: no process-wide lock, because skbb's
// ordination is re-entrant once seeded per call — which is exactly what lets a
// wave's batches ordinate concurrently. `omp_threads` is that call's OpenMP
// width, so the serial path spends the caller's whole thread budget on one
// ordination while a parallel wave gives each worker one thread.
//
// The scope also resolves the seed: a caller-supplied seed >= 0 is passed through
// (so a seeded run is reproducible and identical at any worker count), and an
// unseeded run gets a fresh non-negative seed per call instead of reaching for
// skbb's process-global generator, which concurrent workers would race on.
std::vector<double> PcoaBlock(const DistanceBlock &block, uint32_t d, int seed, int omp_threads,
                              std::vector<double> *eig_out, std::vector<double> *prop_out) {
	const uint32_t m = static_cast<uint32_t>(block.ids.size());
	std::vector<float> eig(d), samples(static_cast<size_t>(m) * d), prop(d);
	{
		miint::unifrac::ComputeCallScope skbb(omp_threads, seed);
		skbb_pcoa_fsvd_fp32(m, block.matrix.data(), d, skbb.seed(), eig.data(), samples.data(), prop.data());
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

uint32_t ChooseWaveWidth(size_t n_anchors, uint32_t batch_size, uint32_t n_workers, uint64_t budget_bytes,
                         size_t n_batches) {
	const double a = static_cast<double>(n_anchors);
	const double k = static_cast<double>(batch_size);
	const double m = a + k;
	const double block_bytes = 5.0 * m * m;
	const double per_batch = block_bytes + 20.0 * (a * k + 0.5 * k * k);
	const double worker_reserve = 2.0 * static_cast<double>(n_workers) * block_bytes;
	const double cap = static_cast<double>(std::max<size_t>(n_batches, 1));
	if (per_batch <= 0.0) {
		return static_cast<uint32_t>(std::min<double>(cap, 4294967295.0)); // no per-batch cost to budget
	}
	// A wave never has more workers than batches, so charging for all n_workers
	// would refuse waves that are affordable precisely because they are narrow —
	// the regime large anchor sets live in, where one block is hundreds of MB.
	// Solve the two regimes instead: below n_workers each batch also carries its
	// worker's blocks; above it the workers' share is a fixed cost. The two agree
	// at W = n_workers, so the result stays monotone in the budget.
	const double budget = static_cast<double>(budget_bytes);
	double fits = budget / (per_batch + 2.0 * block_bytes);
	if (fits > static_cast<double>(n_workers)) {
		fits = (budget - worker_reserve) / per_batch;
	}
	if (!(fits >= 1.0)) {
		return 1; // one block at a time — correct, just one scan per batch
	}
	return static_cast<uint32_t>(std::min(fits, cap));
}

ProgressivePcoaResult RunProgressivePcoa(const std::vector<std::string> &anchor_ids,
                                         const std::vector<std::string> &remaining_ids, uint32_t n_dims,
                                         uint32_t batch_size, int seed, int n_threads, const BlockProvider &get_block,
                                         const WavePrefetch &prefetch, uint32_t wave_batches, uint32_t batch_workers) {
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

	// One batch's contribution to the result. Returned rather than appended so a
	// wave's batches can be produced in any order and merged in batch order, which
	// is what makes the parallel path bit-identical to the serial one.
	struct BatchOutput {
		BatchDiagnostic diag;
		std::vector<ProgressiveCoord> coords;
	};

	// One batch: fetch its block, ordinate it, fit it onto the reference frame
	// through the anchor overlap, and return its non-anchor coordinates. Depends on
	// nothing but its arguments and the reference frame — no shared mutable state —
	// so it is safe to run on several threads at once.
	const auto run_one_batch = [&](size_t start, size_t end, int32_t batch_index, int omp_threads) {
		const std::vector<std::string> req = build_request(start, end);
		const DistanceBlock blk = get_block(req);
		ValidateBlockMatchesRequest(blk, req, "batch");
		const uint32_t m = static_cast<uint32_t>(blk.ids.size());
		const std::vector<double> batch_coords = PcoaBlock(blk, d, seed, omp_threads, nullptr, nullptr);

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

		BatchOutput out;

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
			// The batch's POSITION, never a completion counter: it is what joins a
			// coordinate to the diagnostic describing how that coordinate was placed.
			out.diag.batch = batch_index;
			out.diag.n_samples = static_cast<uint32_t>(end - start);
			out.diag.anchor_m2 = Disparity(ref_std.data(), anchor_fit.data(), a, d);
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
		out.coords.reserve(static_cast<size_t>(nb) * d);
		for (uint32_t r = 0; r < nb; ++r) {
			for (uint32_t axis = 0; axis < d; ++axis) {
				out.coords.push_back({nb_ids[r], static_cast<int32_t>(axis), nb_fit[r * d + axis], batch_index});
			}
		}
		return out;
	};

	// Merge one batch's output. Called in batch order in every path, so the emitted
	// row order does not depend on which batch finished first.
	const auto append_batch = [&result](BatchOutput &&out) {
		result.batches.push_back(out.diag);
		result.coords.insert(result.coords.end(), std::make_move_iterator(out.coords.begin()),
		                     std::make_move_iterator(out.coords.end()));
	};

	std::vector<std::pair<size_t, size_t>> batch_ranges;
	for (size_t start = 0; start < total; start += batch_size) {
		batch_ranges.emplace_back(start, std::min(total, start + static_cast<size_t>(batch_size)));
	}
	const size_t wave_width = std::max<size_t>(1, wave_batches);
#ifdef __EMSCRIPTEN__
	// WASM builds have no threads to fan out to (libssu/libskbb are compiled
	// single-threaded there); batches run one at a time whatever the caller asked.
	(void)batch_workers;
	const size_t workers = 1;
#else
	const size_t workers = std::max<size_t>(1, batch_workers);
#endif

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
		const size_t wave_count = wave_end - wave_start;
		const size_t wave_workers = std::min(workers, wave_count);
		if (wave_workers <= 1) {
			for (size_t b = wave_start; b < wave_end; ++b) {
				// One batch at a time, so it may use the caller's whole thread budget.
				append_batch(
				    run_one_batch(batch_ranges[b].first, batch_ranges[b].second, static_cast<int32_t>(b), n_threads));
			}
			continue;
		}

		// Parallel wave. No lock of any kind: each block's ordination takes its own
		// ComputeCallScope (see PcoaBlock), which is safe concurrently because skbb's
		// fsvd is seeded per call. Parallelism therefore comes from W concurrent
		// single-threaded ordinations rather than one wide one — hence the 1 below —
		// and an unrelated query may ordinate at the same time without waiting.
		std::vector<BatchOutput> outputs(wave_count);
		std::vector<std::exception_ptr> errors(wave_count);
		std::atomic<size_t> next {0};
		const auto worker = [&]() {
			for (size_t i = next.fetch_add(1); i < wave_count; i = next.fetch_add(1)) {
				const size_t b = wave_start + i;
				try {
					outputs[i] = run_one_batch(batch_ranges[b].first, batch_ranges[b].second, static_cast<int32_t>(b),
					                           /*omp_threads=*/1);
				} catch (...) {
					errors[i] = std::current_exception();
				}
			}
		};
		// The driver takes one of the worker slots rather than blocking in join(),
		// so W workers means W threads doing work, not W plus an idle one.
		std::vector<std::thread> pool;
		pool.reserve(wave_workers - 1);
		for (size_t w = 1; w < wave_workers; ++w) {
			pool.emplace_back(worker);
		}
		worker();
		for (auto &t : pool) {
			t.join();
		}
		// Report the lowest-numbered failure — the one a serial run would have hit —
		// so a bad input always produces the same message. The wave's other batches
		// still ran; that is wasted work on an error path, and the price of not
		// letting the reported error depend on which worker lost the race.
		for (size_t i = 0; i < wave_count; ++i) {
			if (errors[i]) {
				std::rethrow_exception(errors[i]);
			}
		}
		for (auto &out : outputs) {
			append_batch(std::move(out));
		}
	}

	return result;
}

} // namespace miint::progressive
