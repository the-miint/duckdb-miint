#pragma once

#include "duckdb/common/arrow/arrow.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

//! Build Arrow C Data Interface arrays from C++ vectors, and release them again.
//!
//! The producer half of the interface, with nothing domain-specific in it: an
//! exported array owns its buffers through `private_data` and frees them when the
//! consumer invokes `release`. Only the three layouts this codebase hands across
//! a boundary are here -- Int64, Float64 and Utf8 -- each dense and non-null,
//! which is the layout a borrowing consumer can read without copying.

//! Export `values` as an Arrow `Int64` array (format "l"): [validity, data].
void ExportInt64(ArrowArray &array, ArrowSchema &schema, std::vector<int64_t> values);

//! Export `values` as an Arrow `Float64` array (format "g"): [validity, data].
void ExportFloat64(ArrowArray &array, ArrowSchema &schema, std::vector<double> values);

//! Export `values` as an Arrow `Utf8` array (format "u"): [validity, offsets, data].
//!
//! Offsets hold `n + 1` entries and count BYTES, not characters; row `i` is
//! `chars[offsets[i] .. offsets[i + 1]]`.
void ExportUtf8(ArrowArray &array, ArrowSchema &schema, const std::vector<std::string> &values);

//! Release an array/schema pair if it was ever filled. A released array has a
//! NULL callback, so this is a no-op on an out-slot nothing wrote to.
void ReleaseIfLive(ArrowArray &array, ArrowSchema &schema);

} // namespace miint
