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

#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <Eigen/LU> // determinant(), for the parity check in Rotation()

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

// How many of a wave's batches may actually run at once.
size_t ResolveWorkers(uint32_t batch_workers) {
#ifdef __EMSCRIPTEN__
	// WASM builds have no threads to fan out to (libssu/libskbb are compiled
	// single-threaded there); batches run one at a time whatever the caller asked.
	(void)batch_workers;
	return 1;
#else
	return std::max<size_t>(1, batch_workers);
#endif
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

uint32_t ChooseWaveWidthByOutput(size_t n_batches, uint32_t n_workers, uint64_t bytes_per_batch,
                                 uint64_t budget_bytes) {
	const uint64_t cap = std::max<uint64_t>(n_batches, 1);
	// The floor: never narrower than the pool, never wider than the run.
	uint64_t width = std::min<uint64_t>(cap, std::max<uint32_t>(n_workers, 1));
	if (bytes_per_batch == 0) {
		return static_cast<uint32_t>(cap); // nothing to charge, so nothing to limit
	}
	width = std::max(width, std::min<uint64_t>(cap, budget_bytes / bytes_per_batch));
	return static_cast<uint32_t>(std::min<uint64_t>(width, 4294967295ull));
}

ProgressivePcoaRun::ProgressivePcoaRun(const std::vector<std::string> &anchor_ids,
                                       const std::vector<std::string> &remaining_ids, uint32_t n_dims,
                                       uint32_t batch_size, int seed, int n_threads, BlockProvider get_block,
                                       WavePrefetch prefetch, uint32_t wave_batches, uint32_t batch_workers,
                                       InterruptCheck interrupt)
    : anchor_ids_(anchor_ids), remaining_ids_(remaining_ids), d_(n_dims), a_(static_cast<uint32_t>(anchor_ids.size())),
      batch_size_(batch_size), seed_(seed), n_threads_(n_threads), get_block_(std::move(get_block)),
      prefetch_(std::move(prefetch)), wave_width_(std::max<size_t>(1, wave_batches)),
      workers_(ResolveWorkers(batch_workers)), interrupt_(std::move(interrupt)) {
	if (n_dims < 1) {
		throw std::invalid_argument("progressive_pcoa: n_dims must be >= 1");
	}
	if (batch_size == 0) {
		throw std::invalid_argument("progressive_pcoa: batch_size must be >= 1");
	}
	if (n_threads < 1) {
		throw std::invalid_argument("progressive_pcoa: n_threads must be >= 1");
	}
	// Need a >= d + 1 for both the reference PCoA (loses one dim to centering) and
	// the procrustes fit (d + 1 points determine a d-dimensional transform). Phrase
	// the guard as `a <= d` so it stays correct even if d is near UINT32_MAX (d + 1
	// would wrap): a > d is exactly a >= d + 1 for integers.
	if (a_ <= d_) {
		throw std::invalid_argument("progressive_pcoa: need at least n_dims + 1 anchors to define the reference frame");
	}

	// Anchors must be distinct and disjoint from the streamed samples: the fit
	// pairs anchor rows 1:1 across the reference and each batch, so a duplicate or
	// an anchor that also appears as a batch sample would silently corrupt the
	// alignment (and double-emit that sample).
	anchor_set_.reserve(a_);
	for (const auto &id : anchor_ids_) {
		if (!anchor_set_.insert(id).second) {
			throw std::invalid_argument("progressive_pcoa: duplicate anchor id '" + id + "'");
		}
	}
	for (const auto &rid : remaining_ids_) {
		if (anchor_set_.count(rid)) {
			throw std::invalid_argument("progressive_pcoa: sample '" + rid + "' is both an anchor and a batch sample");
		}
	}

	for (size_t start = 0; start < remaining_ids_.size(); start += batch_size_) {
		batch_ranges_.emplace_back(start, std::min(remaining_ids_.size(), start + static_cast<size_t>(batch_size_)));
	}
}

// Block request = this batch's non-anchor samples followed by all anchors. Built
// in one place so the request announced by prefetch is byte-identical to the one
// get_block is then called with — a provider is entitled to key its cache on it.
std::vector<std::string> ProgressivePcoaRun::BuildRequest(size_t start, size_t end) const {
	std::vector<std::string> req;
	req.reserve((end - start) + a_);
	req.insert(req.end(), remaining_ids_.begin() + start, remaining_ids_.begin() + end);
	req.insert(req.end(), anchor_ids_.begin(), anchor_ids_.end());
	return req;
}

// ── Principal axes of the assembled configuration ────────────────────────────────
PrincipalAxisAccumulator::PrincipalAxisAccumulator(uint32_t d)
    : d_(d), sum_(d, 0.0), sum_xx_(static_cast<size_t>(d) * d, 0.0) {
	if (d == 0) {
		throw std::invalid_argument("progressive_pcoa: principal-axis accumulator needs d >= 1");
	}
}

void PrincipalAxisAccumulator::Add(const double *coords, size_t n_rows) {
	for (size_t r = 0; r < n_rows; ++r) {
		const double *y = coords + r * d_;
		for (uint32_t a = 0; a < d_; ++a) {
			sum_[a] += y[a];
			// Upper triangle only; symmetrized when the covariance is formed.
			for (uint32_t b = a; b < d_; ++b) {
				sum_xx_[static_cast<size_t>(a) * d_ + b] += y[a] * y[b];
			}
		}
	}
	n_ += n_rows;
}

std::vector<double> PrincipalAxisAccumulator::Mean() const {
	if (n_ == 0) {
		return {};
	}
	std::vector<double> mean(d_);
	for (uint32_t a = 0; a < d_; ++a) {
		mean[a] = sum_[a] / static_cast<double>(n_);
	}
	return mean;
}

std::vector<double> PrincipalAxisAccumulator::Rotation() const {
	if (n_ < 2) {
		throw std::invalid_argument("progressive_pcoa: principal axes need at least 2 samples, got " +
		                            std::to_string(n_));
	}
	const auto mean = Mean();
	const double n = static_cast<double>(n_);
	// C = E[yy^T] - mean mean^T. Start() standardizes the anchor frame to unit
	// Frobenius norm and every batch is mapped INTO that frame, so a coordinate is
	// O(1/sqrt(a*d)) and near-centred no matter how many samples the run has. Both
	// terms are therefore the same small size and the subtraction cannot cancel the
	// way it would on raw, off-origin data.
	Eigen::MatrixXd C(d_, d_);
	for (uint32_t a = 0; a < d_; ++a) {
		for (uint32_t b = a; b < d_; ++b) {
			const double v = sum_xx_[static_cast<size_t>(a) * d_ + b] / n - mean[a] * mean[b];
			C(a, b) = v;
			C(b, a) = v;
		}
	}
	// Symmetric by construction, so a self-adjoint solver — which also guarantees
	// real eigenvalues in ASCENDING order, hence the reversal below.
	Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(C);
	if (es.info() != Eigen::Success) {
		throw std::invalid_argument("progressive_pcoa: principal-axis eigendecomposition did not converge");
	}
	const Eigen::MatrixXd &V = es.eigenvectors();

	// R's rows are orthonormal by construction, so det(R) = +/-1. Both are orthogonal,
	// but only +1 is a ROTATION: reversing the eigenvector columns for descending
	// variance is an odd permutation whenever d_ is even, and the per-axis sign pinning
	// below is free to flip an odd number of them, so nothing so far rules out -1 — a
	// MIRRORED configuration. Procrustes M2 would not notice (it admits reflection) but
	// anything chirality-sensitive would, and "rotation" would be the wrong word for
	// what this returns. Parity is therefore corrected at the end, on the LAST axis:
	// its sign is the least consequential to flip, and doing it there leaves the
	// leading axes' sign convention exactly as pinned.
	std::vector<double> R(static_cast<size_t>(d_) * d_);
	for (uint32_t a = 0; a < d_; ++a) {
		const uint32_t src = d_ - 1 - a; // descending variance
		// Pin the sign: largest-magnitude component positive, lowest index on a tie.
		uint32_t pivot = 0;
		double best = std::abs(V(0, src));
		for (uint32_t k = 1; k < d_; ++k) {
			if (std::abs(V(k, src)) > best) {
				best = std::abs(V(k, src));
				pivot = k;
			}
		}
		const double sign = (V(pivot, src) < 0.0) ? -1.0 : 1.0;
		// Row a of R is principal direction a, so y' = R (y - mean) is a plain matvec.
		for (uint32_t k = 0; k < d_; ++k) {
			R[static_cast<size_t>(a) * d_ + k] = sign * V(k, src);
		}
	}
	// Make it a proper rotation (see above). Deterministic: the parity is a property of
	// the data, not of evaluation order.
	if (Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(R.data(), d_, d_)
	        .determinant() < 0.0) {
		const size_t last = static_cast<size_t>(d_ - 1) * d_;
		for (uint32_t k = 0; k < d_; ++k) {
			R[last + k] = -R[last + k];
		}
	}
	return R;
}

// ── Phase 0: reference PCoA on the anchors — defines the common frame ────────────
std::vector<ProgressiveCoord> ProgressivePcoaRun::Start() {
	if (started_) {
		throw std::logic_error("progressive_pcoa: Start() may only be called once per run");
	}
	// Polled before the anchor block, not after: that block plus its ordination is
	// the run's serial phase (seconds to minutes), so an already-cancelled query
	// should not pay for it.
	if (interrupt_) {
		interrupt_();
	}
	ref_block_ = get_block_(anchor_ids_);
	ValidateBlockMatchesRequest(ref_block_, anchor_ids_, "anchor");

	// Reference coordinates in ref_block_.ids order. This is also where the reported
	// eigenvalues / proportions come from (the anchor ordination — a documented
	// caveat: they describe the anchors, not the full data).
	ref_coords_ = PcoaBlock(ref_block_, d_, seed_, n_threads_, &eigvals_, &proportion_explained_);

	// Emit the reference anchor coordinates once, standardized into the shared
	// frame. ApplyToReference uses only the fit's reference translate/norm, which
	// depend solely on ref_coords_, so a self-fit produces the identical
	// standardization every batch fit would — and works even with zero batches.
	// Every batch below aligns its samples into this same standardized frame.
	const ProcrustesFit ref_fit = FitProcrustes(ref_coords_.data(), ref_coords_.data(), a_, d_);
	std::vector<double> ref_std(static_cast<size_t>(a_) * d_);
	ApplyToReference(ref_fit, ref_coords_.data(), a_, ref_std.data());
	std::vector<ProgressiveCoord> out;
	out.reserve(static_cast<size_t>(a_) * d_);
	for (uint32_t r = 0; r < a_; ++r) {
		for (uint32_t axis = 0; axis < d_; ++axis) {
			out.push_back({ref_block_.ids[r], static_cast<int32_t>(axis), ref_std[r * d_ + axis]});
		}
	}
	// Last, so a Start() that threw leaves the run unstarted rather than half-framed:
	// NextWave then fails loud instead of fitting batches onto nothing.
	started_ = true;
	return out;
}

// One batch: fetch its block, ordinate it, fit it onto the reference frame through
// the anchor overlap, and return its non-anchor coordinates. Depends on nothing but
// its arguments and the reference frame (fixed by Start(), read-only thereafter) —
// no shared mutable state — so it is safe to run on several threads at once.
BatchOutput ProgressivePcoaRun::RunOneBatch(size_t start, size_t end, int32_t batch_index, int omp_threads) const {
	// Polled per batch rather than per wave: a wave can be dozens of batches wide,
	// and a cancellation that waits for the whole wave is not a cancellation. Called
	// from the worker that will run this batch, so the hook must be thread-safe.
	if (interrupt_) {
		interrupt_();
	}
	const std::vector<std::string> req = BuildRequest(start, end);
	const DistanceBlock blk = get_block_(req);
	ValidateBlockMatchesRequest(blk, req, "batch");
	const uint32_t m = static_cast<uint32_t>(blk.ids.size());
	const std::vector<double> batch_coords = PcoaBlock(blk, d_, seed_, omp_threads, nullptr, nullptr);

	// id → batch row index.
	std::unordered_map<std::string, uint32_t> batch_row;
	batch_row.reserve(m);
	for (uint32_t r = 0; r < m; ++r) {
		batch_row.emplace(blk.ids[r], r);
	}

	// Anchor rows of this batch (raw), ordered to pair 1:1 with ref_coords_ by
	// ref_block_.ids order — so batch_anchor row r corresponds to reference row r.
	std::vector<double> batch_anchor(static_cast<size_t>(a_) * d_);
	for (uint32_t r = 0; r < a_; ++r) {
		auto it = batch_row.find(ref_block_.ids[r]);
		if (it == batch_row.end()) {
			throw std::invalid_argument("progressive_pcoa: anchor '" + ref_block_.ids[r] +
			                            "' missing from batch block");
		}
		const double *src = batch_coords.data() + static_cast<size_t>(it->second) * d_;
		std::copy(src, src + d_, batch_anchor.begin() + static_cast<size_t>(r) * d_);
	}

	// Fit the batch onto the reference through the anchor overlap.
	const ProcrustesFit fit = FitProcrustes(ref_coords_.data(), batch_anchor.data(), a_, d_);

	BatchOutput out;

	// Score that fit and report it. This is the run's own accuracy evidence
	// (see BatchDiagnostic): the anchors are the only samples this batch and
	// the reference frame have in common, so the disparity over them is what
	// says whether the batch's placement can be trusted. Costs O(a·d) against
	// a block PCoA, i.e. nothing.
	{
		std::vector<double> ref_std(static_cast<size_t>(a_) * d_);
		std::vector<double> anchor_fit(static_cast<size_t>(a_) * d_);
		ApplyToReference(fit, ref_coords_.data(), a_, ref_std.data());
		ApplyToOther(fit, batch_anchor.data(), a_, anchor_fit.data());
		// The batch's POSITION, never a completion counter: it is what joins a
		// coordinate to the diagnostic describing how that coordinate was placed.
		out.diag.batch = batch_index;
		out.diag.n_samples = static_cast<uint32_t>(end - start);
		out.diag.anchor_m2 = Disparity(ref_std.data(), anchor_fit.data(), a_, d_);
	}

	// This batch's non-anchor rows (raw), then mapped into the reference frame.
	std::vector<std::string> nb_ids;
	std::vector<double> nb_raw;
	nb_ids.reserve(end - start);
	nb_raw.reserve((end - start) * static_cast<size_t>(d_));
	for (uint32_t r = 0; r < m; ++r) {
		if (anchor_set_.count(blk.ids[r])) {
			continue; // anchors already emitted via the reference frame
		}
		nb_ids.push_back(blk.ids[r]);
		const double *src = batch_coords.data() + static_cast<size_t>(r) * d_;
		nb_raw.insert(nb_raw.end(), src, src + d_);
	}
	const uint32_t nb = static_cast<uint32_t>(nb_ids.size());
	std::vector<double> nb_fit(static_cast<size_t>(nb) * d_);
	ApplyToOther(fit, nb_raw.data(), nb, nb_fit.data());
	out.coords.reserve(static_cast<size_t>(nb) * d_);
	for (uint32_t r = 0; r < nb; ++r) {
		for (uint32_t axis = 0; axis < d_; ++axis) {
			out.coords.push_back({nb_ids[r], static_cast<int32_t>(axis), nb_fit[r * d_ + axis], batch_index});
		}
	}
	return out;
}

// ── Phase 1: one wave of the remaining samples ───────────────────────────────────
// Batches are grouped into waves: before fetching a wave's blocks we announce all
// of its requests, so a provider can serve them from one pass over its source
// instead of one pass per block (see WavePrefetch). The batch bodies are unchanged
// by this — a wave is purely an I/O grouping, and W has no effect on the
// coordinates produced.
std::vector<BatchOutput> ProgressivePcoaRun::NextWave() {
	if (!started_) {
		throw std::logic_error("progressive_pcoa: Start() must run before NextWave() — every batch is fitted onto the "
		                       "reference frame it builds");
	}
	if (Done()) {
		return {};
	}
	if (interrupt_) {
		interrupt_();
	}
	const size_t wave_start = next_wave_start_;
	const size_t wave_end = std::min(batch_ranges_.size(), wave_start + wave_width_);
	// Announced for EVERY wave, including a trailing one-batch wave: a uniform
	// contract lets a provider serve all batch blocks from its wave cache and
	// keep a direct path only for the one-off anchor (reference) block. The
	// alternative — skipping the announcement when there is nothing to amortize
	// — would force every provider to implement both paths for batch blocks.
	if (prefetch_) {
		std::vector<std::vector<std::string>> wave_requests;
		wave_requests.reserve(wave_end - wave_start);
		for (size_t b = wave_start; b < wave_end; ++b) {
			wave_requests.push_back(BuildRequest(batch_ranges_[b].first, batch_ranges_[b].second));
		}
		prefetch_(wave_requests);
	}
	const size_t wave_count = wave_end - wave_start;
	const size_t wave_workers = std::min(workers_, wave_count);
	// Returned rather than appended to a shared buffer so a wave's batches can be
	// produced in any order and handed back in batch order, which is what makes the
	// parallel path bit-identical to the serial one.
	std::vector<BatchOutput> outputs(wave_count);
	if (wave_workers <= 1) {
		for (size_t i = 0; i < wave_count; ++i) {
			const size_t b = wave_start + i;
			// One batch at a time, so it may use the caller's whole thread budget.
			outputs[i] =
			    RunOneBatch(batch_ranges_[b].first, batch_ranges_[b].second, static_cast<int32_t>(b), n_threads_);
		}
	} else {
		// Parallel wave. No lock of any kind: each block's ordination takes its own
		// ComputeCallScope (see PcoaBlock), which is safe concurrently because skbb's
		// fsvd is seeded per call. Parallelism therefore comes from W concurrent
		// single-threaded ordinations rather than one wide one — hence the 1 below —
		// and an unrelated query may ordinate at the same time without waiting.
		std::vector<std::exception_ptr> errors(wave_count);
		std::atomic<size_t> next {0};
		const auto worker = [&]() {
			for (size_t i = next.fetch_add(1); i < wave_count; i = next.fetch_add(1)) {
				const size_t b = wave_start + i;
				try {
					outputs[i] = RunOneBatch(batch_ranges_[b].first, batch_ranges_[b].second, static_cast<int32_t>(b),
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
	}
	// Advanced only on success, so a wave that threw is not silently skipped.
	next_wave_start_ = wave_end;
	return outputs;
}

ProgressivePcoaResult RunProgressivePcoa(const std::vector<std::string> &anchor_ids,
                                         const std::vector<std::string> &remaining_ids, uint32_t n_dims,
                                         uint32_t batch_size, int seed, int n_threads, const BlockProvider &get_block,
                                         const WavePrefetch &prefetch, uint32_t wave_batches, uint32_t batch_workers) {
	ProgressivePcoaRun run(anchor_ids, remaining_ids, n_dims, batch_size, seed, n_threads, get_block, prefetch,
	                       wave_batches, batch_workers);
	ProgressivePcoaResult result;
	result.d = n_dims;
	result.coords = run.Start();
	result.eigvals = run.eigvals();
	result.proportion_explained = run.proportion_explained();
	result.coords.reserve(static_cast<size_t>(anchor_ids.size() + remaining_ids.size()) * n_dims);
	result.batches.reserve(run.n_batches());
	while (!run.Done()) {
		std::vector<BatchOutput> wave = run.NextWave();
		// In batch order in every path, so the emitted row order does not depend on
		// which batch finished first.
		for (auto &out : wave) {
			result.batches.push_back(out.diag);
			result.coords.insert(result.coords.end(), std::make_move_iterator(out.coords.begin()),
			                     std::make_move_iterator(out.coords.end()));
		}
	}
	return result;
}

} // namespace miint::progressive
