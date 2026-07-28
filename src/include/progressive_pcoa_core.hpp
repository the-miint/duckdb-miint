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
using BlockProvider = std::function<DistanceBlock(const std::vector<std::string> &requested)>;

// One (sample_id, axis) coordinate in the shared standardized reference frame.
struct ProgressiveCoord {
	std::string sample_id;
	int32_t axis = 0;
	double coordinate = 0.0;
};

struct ProgressivePcoaResult {
	std::vector<ProgressiveCoord> coords;     // anchors (once) + every batch's non-anchor coords
	std::vector<double> eigvals;              // size d — the anchor reference PCoA's eigenvalues
	std::vector<double> proportion_explained; // size d — the anchor reference PCoA's proportions
	uint32_t d = 0;                           // number of ordination axes actually emitted
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
                                         uint32_t batch_size, int seed, int n_threads, const BlockProvider &get_block);

} // namespace miint::progressive
