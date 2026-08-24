#pragma once

#include <cstdint>
#include <functional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

// Progressive reference-anchored PCoA core — a scalable ordination that never
// materializes the full N×N distance matrix. A small set of shared "anchor"
// samples defines a common reference frame; the remaining samples are streamed
// in contiguous batches, each ordinated together with the anchors and then
// aligned back onto the reference by a partial procrustes fit on the anchor
// rows. Peak memory is one (batch + anchors)² block; cost is O(N · batch_size).
//
// This is the shared orchestration used by both progressive_pcoa_from_distances
// (block = a slice of an existing COO distance relation) and
// progressive_pcoa_from_unifrac (block = UniFrac computed on the fly). The ONLY
// difference between them is how a block is obtained, injected here as a
// BlockProvider so the core stays block-source-agnostic.
//
// The partial-procrustes technique (fit the similarity transform on the anchor
// overlap, apply it to a larger set) is McDonald, q2-diversity#338; the
// alignment math is reused from miint::procrustes (A1),
// itself ported from SciPy (BSD-3-Clause). PCoA is delegated to libskbb's
// randomized FSVD (skbb_pcoa_fsvd_fp32) — the same call the pcoa/unifrac_pcoa
// table functions use — so this header (and its implementation) is only compiled
// when UniFrac/skbb support is present (MIINT_HAS_UNIFRAC).
namespace miint::progressive {

// A dense symmetric distance block over a set of samples, in whatever physical
// row/col order the provider chose. `matrix` is m×m row-major fp32 — the exact
// layout skbb_pcoa_fsvd_fp32 consumes. The core maps between `ids` and matrix
// rows by id, so the provider is free to return the samples in any order as long
// as it returns exactly the requested set.
struct DistanceBlock {
	std::vector<std::string> ids; // size m — physical order of the block's rows/cols
	std::vector<float> matrix;    // m*m row-major fp32
};

// Given the requested sample ids, return the dense distance block over exactly
// those samples. Called once for the anchors alone (the reference block) and
// once per batch (that batch's non-anchor ids followed by all anchor ids).
//
// Thread safety: with batch_workers > 1 this is called CONCURRENTLY, once per
// batch of the wave, from worker threads (never for the anchor block, which is
// fetched on the calling thread before any fan-out). A provider that mutates
// shared state — a cache, a scratch buffer — or that calls a library with global
// state must either be made safe or be used with batch_workers = 1.
using BlockProvider = std::function<DistanceBlock(const std::vector<std::string> &requested)>;

// Optional announcement, made before the batches of one "wave" are fetched, of
// every request that wave will make — in the exact order get_block will ask for
// them.
//
// WHY: the natural implementation of BlockProvider for a stored relation is one
// filtered query per block, which re-reads the whole relation for every batch:
// with B batches that is B full scans, i.e. O(B · pairs) = O(N³ / batch_size).
// Announcing a wave lets a provider satisfy many blocks from a SINGLE pass — most
// rows are needed by exactly one block (both endpoints in the same batch), rows
// touching an anchor are needed by one block, anchor×anchor rows by all of them,
// and rows spanning two different batches by none. The wave's width is chosen by
// the caller from a memory budget, since W blocks are held at once.
//
// A provider that ignores this (or a caller that passes nothing) still works —
// get_block alone is sufficient, just at one scan per block.
using WavePrefetch = std::function<void(const std::vector<std::vector<std::string>> &requests)>;

// Optional cooperative-cancellation hook. Polled before the anchor block, at the
// start of each wave (i.e. before its prefetch), and before each batch's block is
// fetched. Throw from it to abort the run — the exception propagates out of the
// call that polled it, unchanged.
//
// A run at the scale this core exists for takes hours, so there has to be a way
// out of it; polling per batch rather than per wave is what bounds the wait to one
// block. It is called from the thread that will run that batch, so with
// batch_workers > 1 it is called CONCURRENTLY and must be thread-safe (reading an
// atomic flag is the intended use). An abort inside a parallel wave still lets the
// wave's in-flight batches finish — their work is discarded — because every other
// batch of that wave polls and throws too.
using InterruptCheck = std::function<void()>;

// How many batches one wave should hold, given the memory a caller is willing to
// spend on it. Lives beside the wave_batches parameter it feeds, because the
// sizing rule is a property of what a wave costs:
//
//   cached blocks   W · 5·(anchors + batch_size)²   fp32 matrix + fill bitmap
//   worker blocks   2 · n_workers · (same)          the block each worker holds
//                                                   plus its ordination workspace
//   scan rows       W · 20·(anchors·batch + batch²/2)  a materialized block-fill
//                                                   query, ~20 B/row
//
// and a wave never has more workers than batches, so the workers' share is
// charged per batch below n_workers and as a fixed cost above it. The result is
// always in [1, max(n_batches, 1)]: a budget that fits nothing still yields 1 (one
// block at a time is correct, just one scan per batch), and no wave is ever wider
// than the run has batches. Arithmetic is done in double so the
// (anchors + batch_size)² term cannot wrap.
//
// W is an I/O choice only — it never changes the coordinates a run produces — so
// this is free to be a heuristic. What it must never return is an illegal width.
uint32_t ChooseWaveWidth(size_t n_anchors, uint32_t batch_size, uint32_t n_workers, uint64_t budget_bytes,
                         size_t n_batches);

// The same decision for a provider that fetches and releases each block INSIDE its
// own batch (no wave staging — the progressive UniFrac path). Nothing about a wave
// is held there except the output of the batches already in it, so that is all this
// charges: `bytes_per_batch` is what one batch's coordinates occupy on the way out.
//
// Width matters because a wave is a barrier. One batch per worker — the obvious
// sizing, and what this path used to do — is its worst case: every barrier waits
// for the wave's slowest batch, so the run pays max(batch) per wave instead of
// mean(batch), and the ragged final wave leaves the rest of the pool idle to the
// end. Widening the wave amortizes both, and a wave wide enough to hold the whole
// run has no barrier at all. See the call site for the measured cost.
//
// `n_workers` is the floor rather than a term: a wave narrower than the pool idles
// workers by construction, and holding a pool's worth of output is what the
// un-widened path already did, so the floor cannot regress memory. Result is always
// in [1, max(n_batches, 1)].
uint32_t ChooseWaveWidthByOutput(size_t n_batches, uint32_t n_workers, uint64_t bytes_per_batch, uint64_t budget_bytes);

// One (sample_id, axis) coordinate in the shared standardized reference frame.
// `batch` is the 0-based index of the batch that placed this sample, or -1 for the
// anchors — they define the frame rather than being fitted into it. It joins a
// coordinate to its BatchDiagnostic, so a caller can attribute the quality of any
// sample's placement without re-deriving the batching.
struct ProgressiveCoord {
	std::string sample_id;
	int32_t axis = 0;
	double coordinate = 0.0;
	int32_t batch = -1;
};

// Per-batch quality evidence. Each batch is placed into the shared frame by a
// procrustes fit on the anchor overlap; `anchor_m2` is that fit's disparity, i.e.
// how well this batch's own view of the anchors agreed with the reference view.
// It is the only accuracy signal available without computing a full PCoA — which
// is by definition impossible at the scale this function exists for — so it is
// reported rather than discarded. ~0 means the batch slotted into the frame
// cleanly; a large value means this batch's samples are poorly determined by the
// anchor set (typically too few anchors, or anchors that don't span the region
// this batch lives in).
struct BatchDiagnostic {
	int32_t batch = 0;      // 0-based batch index, in emission order
	uint32_t n_samples = 0; // non-anchor samples placed by this batch
	double anchor_m2 = 0.0; // procrustes disparity of the anchor-overlap fit
};

// One batch's whole contribution: its coordinates and the diagnostic describing
// how they were placed. Kept together because a coordinate is only interpretable
// alongside its batch's disparity, and because a caller that consumes the run
// incrementally (see ProgressivePcoaRun) never has the full BatchDiagnostic list
// to join against.
struct BatchOutput {
	BatchDiagnostic diag;
	std::vector<ProgressiveCoord> coords;
};

struct ProgressivePcoaResult {
	std::vector<ProgressiveCoord> coords;     // anchors (once) + every batch's non-anchor coords
	std::vector<double> eigvals;              // size d — the anchor reference PCoA's eigenvalues
	std::vector<double> proportion_explained; // size d — the anchor reference PCoA's proportions
	uint32_t d = 0;                           // number of ordination axes actually emitted
	std::vector<BatchDiagnostic> batches;     // one per batch, in order (empty when every sample is an anchor)
};

// Second moments of the assembled configuration, accumulated one batch at a time,
// and the transform that takes it onto its own principal axes.
//
// WHY THIS EXISTS: every batch is fitted into the frame of the ANCHOR ordination, so
// the emitted axes are the anchor block's principal axes — not the assembled
// configuration's. "PC1" is therefore not the leading axis of the output, and the
// reported eigenvalues describe the anchors rather than the emitted coordinates.
// Rotating the finished configuration onto its own principal axes fixes the first of
// those.
//
// WHAT IT DOES AND DOES NOT BUY: a procrustes disparity is rotation-invariant, so
// this cannot move M² against any reference — claiming an accuracy gain from it
// would be wrong. What it moves is AXIS-BY-AXIS agreement, which is what a reader
// interpreting "PC1 vs PC2" actually relies on, because a full PCoA's axes are the
// principal axes of the full configuration.
//
// Only the moments are accumulated (d² + d doubles, exact, streaming); the rotation
// itself needs nothing else. What the caller must still arrange is somewhere to hold
// the coordinates, because applying the transform means revisiting rows that have
// already been emitted.
//
// The transform is y' = R (y − mean): centred as well as rotated, so the result is a
// centred principal-axis configuration like a real PCoA's, rather than a rotated one
// with an off-origin centroid. Translation does not affect a procrustes comparison
// either.
class PrincipalAxisAccumulator {
public:
	explicit PrincipalAxisAccumulator(uint32_t d);

	// `coords` is n_rows × d row-major. Call once per batch, in any order — addition
	// is commutative and the result depends only on the multiset of rows.
	void Add(const double *coords, size_t n_rows);

	uint64_t count() const {
		return n_;
	}

	// Column means (size d). Empty until at least one row has been added.
	std::vector<double> Mean() const;

	// d×d row-major rotation, applied as y' = R (y − mean). Columns are the
	// configuration's principal directions ordered by DESCENDING variance.
	//
	// Eigenvector sign is inherently arbitrary, so it is pinned: each axis is
	// oriented so that its largest-magnitude component is positive (lowest index wins
	// a tie). That is for reproducibility and readability only — it matches no
	// external convention, and a procrustes comparison is indifferent to it.
	//
	// Throws std::invalid_argument when fewer than 2 rows have been added (a
	// covariance needs at least two points).
	std::vector<double> Rotation() const;

private:
	const uint32_t d_;
	uint64_t n_ = 0;
	std::vector<double> sum_;    // d
	std::vector<double> sum_xx_; // d*d row-major
};

// Run progressive reference-anchored PCoA.
//
//   anchor_ids     distinct anchor samples shared across all batches. Must number
//                  at least n_dims + 1 (a d-dimensional procrustes fit and a
//                  d-axis PCoA each need d + 1 points).
//   remaining_ids  the non-anchor samples in streaming (batch) order; must be
//                  disjoint from anchor_ids.
//   n_dims         number of ordination axes to compute (d).
//   batch_size     non-anchor samples per batch (> 0).
//   seed           passed to every skbb PCoA call; use >= 0 for reproducibility.
//   n_threads      OpenMP fan-out for one block ordination (>= 1). Spent on a
//                  single wide ordination when batches run serially; a parallel
//                  wave instead gives each worker one OpenMP thread (see
//                  batch_workers), so this bounds total fan-out either way.
//   get_block      supplies each dense distance block (see BlockProvider).
//   prefetch       optional; when set, called once per wave with that wave's
//                  requests before any of them is fetched (see WavePrefetch).
//   wave_batches   how many batches form one wave (0 or 1 = no batching, i.e. the
//                  provider is asked for one block at a time). The caller sizes
//                  this from its memory budget: a wave holds W blocks at once.
//                  It affects only I/O, never the result — the same anchors, seed
//                  and batch_size produce identical coordinates for any W.
//   batch_workers  how many of a wave's batches to process at once (0 or 1 =
//                  serial). Batches are independent — the frame is fixed before
//                  the loop — so this only reorders execution: results are
//                  bit-identical to a serial run for any worker count, batch
//                  indices stay positional, and errors are reported for the
//                  lowest-numbered failing batch as a serial run would.
//                  Workers are drawn from one wave, so parallelism requires
//                  wave_batches > 1; the effective count is
//                  min(batch_workers, batches in the wave).
//                  Each worker's PCoA is pinned to ONE OpenMP thread and takes no
//                  process-wide lock (see ComputeCallScope), so parallelism comes
//                  from concurrent block ordinations rather than a wider fsvd, and
//                  an unrelated query may ordinate at the same time. Requires a
//                  thread-safe get_block (see BlockProvider).
//
// Flow: PCoA on the anchor block defines the reference frame; its eigenvalues and
// proportions are reported (a documented caveat — they describe the anchors, not
// the full ordination). Each batch's block is ordinated, aligned to the reference
// by FitProcrustes on the anchor rows, and its non-anchor coordinates are emitted
// via ApplyToOther. The anchor coordinates are emitted once, standardized into the
// same frame via ApplyToReference.
//
// Throws std::invalid_argument on n_dims < 1, batch_size == 0, too few anchors, or
// an anchor/remaining overlap; propagates exceptions from get_block and from the
// procrustes core (e.g. non-finite coordinates, degenerate anchor block).
ProgressivePcoaResult RunProgressivePcoa(const std::vector<std::string> &anchor_ids,
                                         const std::vector<std::string> &remaining_ids, uint32_t n_dims,
                                         uint32_t batch_size, int seed, int n_threads, const BlockProvider &get_block,
                                         const WavePrefetch &prefetch = nullptr, uint32_t wave_batches = 0,
                                         uint32_t batch_workers = 1);

// The same run, driven one wave at a time so a caller can consume each batch's
// rows as soon as that batch is placed instead of waiting for the whole run.
// RunProgressivePcoa above is a thin loop over this class and stays the API for
// callers that want the whole result at once.
//
// WHY this exists: a caller that must hold every coordinate before emitting any
// pays for the result twice — once in ProgressivePcoaResult::coords and once in
// its own row representation — and at 10M samples × d = 10 that is tens of GB, two
// copies live at the same time. Consuming per wave bounds the caller's buffer to
// one wave's rows, gives it somewhere to poll for cancellation between batches,
// and lets rows reach the consumer while later batches are still computing.
//
// Contract:
//   - `anchor_ids` and `remaining_ids` are held BY REFERENCE — at 10M samples
//     copying them is the very cost this class exists to avoid. They must outlive
//     the run. Everything else is copied.
//   - Start() exactly once, before any NextWave(); it runs the anchor block and
//     the reference PCoA (the frame every batch is fitted into) and returns the
//     standardized anchor coordinates (batch = -1). eigvals()/proportion_explained()
//     are only meaningful afterwards.
//   - NextWave() runs the next wave and returns its batches in batch order —
//     positionally identical to a serial run, for any worker count. It returns an
//     empty vector once every batch has been run; Done() reports the same thing
//     without running anything.
//   - Constructor argument validation and semantics are RunProgressivePcoa's,
//     unchanged; see its comment above for every parameter.
class ProgressivePcoaRun {
public:
	ProgressivePcoaRun(const std::vector<std::string> &anchor_ids, const std::vector<std::string> &remaining_ids,
	                   uint32_t n_dims, uint32_t batch_size, int seed, int n_threads, BlockProvider get_block,
	                   WavePrefetch prefetch = nullptr, uint32_t wave_batches = 0, uint32_t batch_workers = 1,
	                   InterruptCheck interrupt = nullptr);

	std::vector<ProgressiveCoord> Start();
	std::vector<BatchOutput> NextWave();

	bool Done() const {
		return next_wave_start_ >= batch_ranges_.size();
	}
	uint32_t d() const {
		return d_;
	}
	size_t n_batches() const {
		return batch_ranges_.size();
	}
	// The anchor reference PCoA's eigenvalues / proportions (size d), available
	// after Start(). They describe the anchors, not the full data — the documented
	// caveat of the method.
	const std::vector<double> &eigvals() const {
		return eigvals_;
	}
	const std::vector<double> &proportion_explained() const {
		return proportion_explained_;
	}

private:
	std::vector<std::string> BuildRequest(size_t start, size_t end) const;
	BatchOutput RunOneBatch(size_t start, size_t end, int32_t batch_index, int omp_threads) const;

	const std::vector<std::string> &anchor_ids_;
	const std::vector<std::string> &remaining_ids_;
	const uint32_t d_;
	const uint32_t a_;
	const uint32_t batch_size_;
	const int seed_;
	const int n_threads_;
	const BlockProvider get_block_;
	const WavePrefetch prefetch_;
	const size_t wave_width_;
	const size_t workers_;
	const InterruptCheck interrupt_;
	std::unordered_set<std::string> anchor_set_;
	std::vector<std::pair<size_t, size_t>> batch_ranges_;
	// Fixed by Start(), then read-only — which is what lets a wave's batches run
	// concurrently with no synchronization at all.
	DistanceBlock ref_block_;
	std::vector<double> ref_coords_;
	std::vector<double> eigvals_;
	std::vector<double> proportion_explained_;
	bool started_ = false;
	size_t next_wave_start_ = 0;
};

} // namespace miint::progressive
