#pragma once

// Guarded wrapper around unifrac-binaries' api.hpp.
//
// Upstream api.hpp ships without #include guards or #pragma once, so any
// translation unit that pulls it in twice (e.g. via two of our adapter
// headers) hits redefinition errors on every typedef. This header is the
// single allowed include site for api.hpp inside duckdb-miint.

#include "api.hpp"
