#pragma once

// What miint's installers (`install_gpl_boundary()`, `install_skb_arrow()`)
// return, and the `force` argument they take. Header-only, so the WASM and
// Windows stubs, which compile neither installer, return the same STRUCT.

#include "duckdb.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/string_vector.hpp"
#include "duckdb/common/vector/struct_vector.hpp"

namespace duckdb {

/// STRUCT(installed BOOLEAN, path VARCHAR, version VARCHAR, message VARCHAR).
/// The field order is part of the API: callers project by position.
inline LogicalType InstallReportType() {
	return LogicalType::STRUCT({{"installed", LogicalType::BOOLEAN},
	                            {"path", LogicalType::VARCHAR},
	                            {"version", LogicalType::VARCHAR},
	                            {"message", LogicalType::VARCHAR}});
}

/// The `force` argument: false for the 0-arg overload, and for NULL (the
/// convention of `mask_dust(seq, NULL)`). Read from the first row only; the
/// installer runs once per call.
inline bool InstallForceArg(DataChunk &args) {
	if (args.ColumnCount() == 0) {
		return false;
	}
	UnifiedVectorFormat fmt;
	args.data[0].ToUnifiedFormat(fmt);
	const auto idx = fmt.sel->get_index(0);
	if (!fmt.validity.RowIsValid(idx)) {
		return false;
	}
	return UnifiedVectorFormat::GetData<bool>(fmt)[idx];
}

/// An installer's overloads: 0-arg (find, else install) and 1-arg BOOLEAN
/// `force`. The stubs build theirs here too, so SQL that passes `force`
/// type-checks on every build. Volatile: an installer has side effects, so
/// DuckDB must never fold one into a plan, where EXPLAIN or PREPARE would run it.
inline ScalarFunctionSet InstallFunctionSet(const char *name, const scalar_function_t &execute) {
	ScalarFunctionSet set(name);
	for (const auto &arguments : {vector<LogicalType> {}, vector<LogicalType> {LogicalType::BOOLEAN}}) {
		ScalarFunction fn(arguments, InstallReportType(), execute);
		fn.SetVolatile();
		set.AddFunction(fn);
	}
	return set;
}

/// Make `result` a constant vector of one report, so `SELECT installer()
/// FROM range(n)` installs once, not n times.
inline void SetInstallReport(Vector &result, bool installed, const string &path, const string &version,
                             const string &message) {
	auto &entries = StructVector::GetEntries(result);
	FlatVector::GetDataMutable<bool>(entries[0])[0] = installed;
	FlatVector::GetDataMutable<string_t>(entries[1])[0] = StringVector::AddString(entries[1], path);
	FlatVector::GetDataMutable<string_t>(entries[2])[0] = StringVector::AddString(entries[2], version);
	FlatVector::GetDataMutable<string_t>(entries[3])[0] = StringVector::AddString(entries[3], message);
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
}

} // namespace duckdb
