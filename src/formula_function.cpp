#include "formula_function.hpp"

#include "formula_parser.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>

namespace duckdb {

static void FormulaScalarFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	UnaryExecutor::Execute<string_t, double>(args.data[0], result, args.size(), [](string_t input) {
		try {
			return miint::ParseFormula(input.GetString());
		} catch (const std::runtime_error &e) {
			throw InvalidInputException(e.what());
		}
	});
}

void FormulaFunction::Register(ExtensionLoader &loader) {
	ScalarFunction formula_func("formula", {LogicalType::VARCHAR}, LogicalType::DOUBLE, FormulaScalarFunction);
	loader.RegisterFunction(formula_func);
}

} // namespace duckdb
