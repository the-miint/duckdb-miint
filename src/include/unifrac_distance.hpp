#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "unifrac_bptree.hpp"
#include "unifrac_libssu.hpp" // mat_full_fp32_t
#include "unifrac_support_biom.hpp"

namespace miint::unifrac {

// RAII wrapper around a libssu fp32 distance matrix.
//
// Phase-4 PCoA and Phase-5 PERMANOVA both run the same compute path:
// acquire the libssu RNG mutex, set the seed, call
// one_off_matrix_inmem_fp32_v3, then read mat->matrix / mat->n_samples /
// mat->sample_ids. This wrapper owns the libssu pointer (frees via
// destroy_mat_full_fp32) and copies the sample_ids char** into an owned
// std::vector<std::string> so callers can keep using them after the
// matrix is destroyed by a downstream destructor.
//
// Subsampling note: when subsample_depth > 0, libssu drops samples
// whose total counts fall below the depth. n_samples() and
// sample_ids() reflect the post-subsample view, NOT the input table's
// size. Always trust these values, never the pre-subsample
// `ordered_sample_ids` from the input biom view.
//
// Thread safety: Compute() takes the process-wide libssu/OpenMP mutex
// (via OmpThreadScope) internally — callers don't need to. The mutex is
// held only across ssu_set_random_seed + one_off_matrix_inmem_fp32_v3,
// not for the caller's downstream work on the returned matrix.
//
// `n_threads` pins the libssu OpenMP fan-out (must be >= 1). Callers
// should resolve from DuckDB's TaskScheduler::NumberOfThreads() via
// ResolveThreadsParameter; unit tests pass 1 to keep test wall time
// independent of host core count.
class UnifracDistanceMatrix {
public:
	static UnifracDistanceMatrix Compute(const UnifracSupportBiomView &biom_view, const UnifracBptreeView &bptree_view,
	                                     const std::string &variant_fp32, bool variance_adjust, double alpha,
	                                     bool bypass_tips, bool normalize_sample_counts, uint32_t subsample_depth,
	                                     bool subsample_with_replacement, int seed, int n_threads);

	~UnifracDistanceMatrix();
	UnifracDistanceMatrix(UnifracDistanceMatrix &&other) noexcept;
	UnifracDistanceMatrix &operator=(UnifracDistanceMatrix &&other) noexcept;
	UnifracDistanceMatrix(const UnifracDistanceMatrix &) = delete;
	UnifracDistanceMatrix &operator=(const UnifracDistanceMatrix &) = delete;

	const float *matrix() const {
		return mat_->matrix;
	}
	uint32_t n_samples() const {
		return mat_->n_samples;
	}
	const std::vector<std::string> &sample_ids() const {
		return sample_ids_;
	}

private:
	UnifracDistanceMatrix() = default;
	mat_full_fp32_t *mat_ = nullptr;
	std::vector<std::string> sample_ids_;
};

} // namespace miint::unifrac
