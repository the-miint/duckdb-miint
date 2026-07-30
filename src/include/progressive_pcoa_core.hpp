#pragma once

#include <cstdint>
#include <functional>
#include <string>
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
// overlap, apply it to a larger set) is the author's own work (Daniel McDonald,
// q2-diversity#338); the alignment math is reused from miint::procrustes (A1),
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

struct ProgressivePcoaResult {
	std::vector<ProgressiveCoord> coords;     // anchors (once) + every batch's non-anchor coords
	std::vector<double> eigvals;              // size d — the anchor reference PCoA's eigenvalues
	std::vector<double> proportion_explained; // size d — the anchor reference PCoA's proportions
	uint32_t d = 0;                           // number of ordination axes actually emitted
	std::vector<BatchDiagnostic> batches;     // one per batch, in order (empty when every sample is an anchor)
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
//   n_threads      OpenMP fan-out for each skbb PCoA (>= 1); also drives the
//                  process-wide OmpThreadScope that serializes libssu/skbb calls
//                  against other concurrent DuckDB queries.
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
//                  Each worker's PCoA is pinned to ONE OpenMP thread (the wave
//                  holds a single OmpThreadScope), so parallelism comes from
//                  concurrent block ordinations rather than a wider fsvd — see
//                  OmpThreadPin for why libssu-backed providers must not use it.
//                  Requires a thread-safe get_block (see BlockProvider).
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

} // namespace miint::progressive
