#include "community_distances_gram.hpp"

#include <Eigen/Core>

namespace miint {

void CommunityDistancesGramLower(const double *operand, uint32_t n, uint32_t f, uint32_t n_cached_tail, double *gram) {
	using RowMajor = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
	const Eigen::Map<const RowMajor> X(operand, n, f);
	Eigen::Map<Eigen::MatrixXd> S(gram, n, n);

	// With a cached tail only the head triangle and the tail x head rectangle are
	// formed; a block that is all tail writes nothing.
	const uint32_t head = n - n_cached_tail;
	if (n_cached_tail == 0) {
		S.triangularView<Eigen::Lower>() = X * X.transpose();
	} else if (head > 0) {
		S.topLeftCorner(head, head).triangularView<Eigen::Lower>() = X.topRows(head) * X.topRows(head).transpose();
		S.bottomLeftCorner(n_cached_tail, head).noalias() = X.bottomRows(n_cached_tail) * X.topRows(head).transpose();
	}
}

} // namespace miint
