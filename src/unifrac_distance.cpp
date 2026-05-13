#include "unifrac_distance.hpp"

#include <mutex>
#include <stdexcept>
#include <string>

namespace miint::unifrac {

namespace {

// Process-wide serialization for libssu calls that touch the global skbb
// RNG. one_off_matrix_inmem_fp32_v3 with subsample_depth > 0 reads the
// RNG state set by ssu_set_random_seed; concurrent invocations from
// multiple DuckDB connections would otherwise interleave seed updates.
// We always hold the mutex around the libssu call even when
// subsample_depth == 0 — uniform call site, negligible cost relative to
// the distance computation itself.
std::mutex &GlobalRngMutex() {
	static std::mutex m;
	return m;
}

const char *ComputeStatusName(ComputeStatus s) {
	switch (s) {
	case okay:
		return "okay";
	case tree_missing:
		return "tree_missing";
	case table_missing:
		return "table_missing";
	case table_empty:
		return "table_empty";
	case unknown_method:
		return "unknown_method";
	case table_and_tree_do_not_overlap:
		return "table_and_tree_do_not_overlap";
	case output_error:
		return "output_error";
	case invalid_method:
		return "invalid_method";
	case grouping_missing:
		return "grouping_missing";
	}
	return "unknown";
}

} // namespace

UnifracDistanceMatrix UnifracDistanceMatrix::Compute(const UnifracSupportBiomView &biom_view,
                                                     const UnifracBptreeView &bptree_view,
                                                     const std::string &variant_fp32, bool variance_adjust,
                                                     double alpha, bool bypass_tips, bool normalize_sample_counts,
                                                     uint32_t subsample_depth, bool subsample_with_replacement,
                                                     int seed) {
	mat_full_fp32_t *mat = nullptr;
	{
		std::lock_guard<std::mutex> lock(GlobalRngMutex());
		if (seed >= 0) {
			ssu_set_random_seed(static_cast<unsigned int>(seed));
		}
		ComputeStatus s = one_off_matrix_inmem_fp32_v3(
		    biom_view.support_biom(), bptree_view.support_bptree(), variant_fp32.c_str(), variance_adjust, alpha,
		    bypass_tips, normalize_sample_counts, /*n_substeps*/ 1, subsample_depth, subsample_with_replacement,
		    /*mmap_dir*/ nullptr, &mat);
		if (s != okay) {
			throw std::runtime_error(std::string("libssu one_off_matrix_inmem_fp32_v3 returned status ") +
			                         ComputeStatusName(s));
		}
	}

	UnifracDistanceMatrix out;
	out.mat_ = mat;
	out.sample_ids_.reserve(mat->n_samples);
	for (uint32_t i = 0; i < mat->n_samples; ++i) {
		out.sample_ids_.emplace_back(mat->sample_ids[i]);
	}
	return out;
}

UnifracDistanceMatrix::~UnifracDistanceMatrix() {
	if (mat_ != nullptr) {
		destroy_mat_full_fp32(&mat_);
		mat_ = nullptr;
	}
}

UnifracDistanceMatrix::UnifracDistanceMatrix(UnifracDistanceMatrix &&other) noexcept
    : mat_(other.mat_), sample_ids_(std::move(other.sample_ids_)) {
	other.mat_ = nullptr;
}

UnifracDistanceMatrix &UnifracDistanceMatrix::operator=(UnifracDistanceMatrix &&other) noexcept {
	if (this != &other) {
		if (mat_ != nullptr) {
			destroy_mat_full_fp32(&mat_);
		}
		mat_ = other.mat_;
		sample_ids_ = std::move(other.sample_ids_);
		other.mat_ = nullptr;
	}
	return *this;
}

} // namespace miint::unifrac
