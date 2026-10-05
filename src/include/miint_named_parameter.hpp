#pragma once

#include "duckdb/function/table_function.hpp"

namespace duckdb {

//! Declares a named option on a table function.
//! DuckDB v2.0 replaced `TableFunction::named_parameters` with a typed "**kwargs" parameter on the signature.
//! A signature may carry only one, and ExtendTypedKwargs throws until it exists, so the first option creates it
//! and later ones (including those added by shared helpers) extend it.
inline void AddNamedParameter(TableFunction &fn, const char *name, LogicalType type) {
	auto add = [&](TypedKwargs &options) {
		options.Add(name, std::move(type));
	};
	auto &signature = fn.GetSignature();
	if (signature.GetTypedKwargs()) {
		signature.ExtendTypedKwargs(add);
	} else {
		signature.WithTypedKwargs("options", add);
	}
}

} // namespace duckdb
