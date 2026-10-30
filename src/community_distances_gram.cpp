#include "community_distances_gram.hpp"

#include <Eigen/Core>

namespace miint {

void CommunityDistancesGramLower(const double *operand, uint32_t n, uint32_t f, uint32_t n_cached_tail, double *gram) {
	using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	const Eigen::Map<const RowMajor> X(operand, n, f);
	Eigen::Map<Eigen::MatrixXd> S(gram, n, n);

	// With a cached tail the product splits into the head triangle and the tail x head
	// rectangle, leaving the tail's own triangle uncomputed. `triangularView<Lower>() =
	// product` dispatches to the same kernel with the same blocking as rankUpdate, so the
	// values are bit-identical.
	const uint32_t head = n - n_cached_tail;
	if (n_cached_tail == 0) {
		S.triangularView<Eigen::Lower>() = X * X.transpose();
	} else if (head > 0) {
		// Reached only when n_cached_tail != 0 (the branch above takes 0), so the tail
		// rectangle is mandatory here, not optional. `head > 0` is defensive: the sparse
		// entry point returns early when the whole block is cached, so a zero head does
		// not arrive from there.
		S.topLeftCorner(head, head).triangularView<Eigen::Lower>() = X.topRows(head) * X.topRows(head).transpose();
		S.bottomLeftCorner(n_cached_tail, head).noalias() = X.bottomRows(n_cached_tail) * X.topRows(head).transpose();
	}
}

} // namespace miint
