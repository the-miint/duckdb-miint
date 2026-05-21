#include "qc_functions.hpp"

#include "duckdb/common/vector_operations/vector_operations.hpp"
#include "duckdb/function/scalar_function.hpp"

namespace duckdb {

static constexpr const char *QC_VERSION_STRING = "qc 0.1.0 (port of fastp algorithms)";

static void QcVersionFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto out = ConstantVector::GetData<string_t>(result);
	out[0] = StringVector::AddString(result, QC_VERSION_STRING);
}

void QcFunctions::Register(ExtensionLoader &loader) {
	ScalarFunction qc_version("qc_version", {}, LogicalType::VARCHAR, QcVersionFunction);
	loader.RegisterFunction(qc_version);
}

} // namespace duckdb
