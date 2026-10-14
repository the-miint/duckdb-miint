#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_version() -> VARCHAR: the version of the linked sc random-forest library.
//!
//! Deliberately the thinnest possible binding to sc: no Arrow, no context, no
//! model. It exists so a build can prove the sc_* symbols are present and
//! callable before any of the COO/Arrow marshaling lands on top of them.
class ScVersionFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
