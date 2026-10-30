#pragma once
//
// The Eigen product behind community_distances.cpp's dense-Gram path. It has its own
// file because community_distances.cpp is built with -ffp-contract=off and so must not
// instantiate Eigen's templates (CMakeLists.txt).

#include <cstdint>

namespace miint {

//! Writes the lower triangle of X * X^T into `gram`, an n x n column-major buffer,
//! where X is the row-major n x f `operand`. With `n_cached_tail` > 0 the tail x tail
//! triangle is left unwritten, as is the strict upper triangle in every case.
void CommunityDistancesGramLower(const double *operand, uint32_t n, uint32_t f, uint32_t n_cached_tail, double *gram);

} // namespace miint
