#include "sc_version_function.hpp"

#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include "sc.h"

namespace duckdb {

static void ScVersionScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	// sc_version() returns a pointer into a compile-time byte string that is
	// valid for the whole program lifetime and must NOT be freed (see sc.h).
	// AddString copies into the result vector's own string heap, so the
	// returned string_t does not alias sc's static storage.
	const char *version = sc_version();

	// Zero-argument and value-stable for the life of the process, so one
	// constant vector answers every row in the chunk.
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::GetData<string_t>(result)[0] = StringVector::AddString(result, version);
}

void ScVersionFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(ScalarFunction("sc_version", {}, LogicalType::VARCHAR, ScVersionScalarFunction));
}

} // namespace duckdb
