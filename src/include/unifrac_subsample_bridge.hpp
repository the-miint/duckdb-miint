#pragma once

#include <cstdint>

#include "unifrac_support_biom.hpp"

namespace miint::unifrac {

// Subsample a UnifracSupportBiomView and return a fresh
// UnifracSupportBiomView over the subsampled table.
//
// Why a bridge: libssu's subsample_table_inmem_seeded returns an
// opaque_biom_inmem_t* whose accessors expose
// (subsampled_n_samples, subsampled_n_obs, subsampled_get_sample_id,
//  subsampled_get_obs_id, subsampled_get_obs_data) but no direct
// support_biom_t. faith_pd_inmem (and any future support_biom_t-only
// API) needs that flat layout, so we materialize the opaque table's
// contents through COO rows into a new UnifracSupportBiomView,
// reusing the Phase-1 adapter's canonical-ordering and dedup logic.
//
// Behaviour:
// - seed >= 0: deterministic; same seed → byte-identical reconstruction.
// - seed < 0: libssu falls back to the global skbb RNG. Callers that
//   need reproducibility must always pass a non-negative seed.
// - subsample_table_inmem_seeded drops samples whose total count is
//   below `depth`. The returned view's n_samples may be smaller than
//   the input view's n_samples.
// - OTUs that end up with zero counts in every surviving sample are
//   dropped (preserves the sparse-storage invariant).
//
// Cost: O(n_obs × n_samples) — we walk every cell of the opaque table.
// Acceptable for the v1 single-table Faith PD path; revisit if profiling
// shows it dominates.
//
// Throws std::runtime_error on libssu compute errors. Throws
// std::invalid_argument if the subsampled table is empty.
UnifracSupportBiomView BridgeSubsample(const UnifracSupportBiomView &biom_view, uint32_t subsample_depth,
                                       bool with_replacement, int seed);

} // namespace miint::unifrac
