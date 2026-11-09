#include "lineage_rank.hpp"

#include "taxonomy_lineage.hpp"

#include "duckdb/common/optional.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"

#include <string>

namespace duckdb {

// The mapping is FindRankColumn's, the one read_ncbi_lineage collapses with and the
// taxonomy_lineage macro calls this function for, so the three cannot disagree.
static void LineageRankExecute(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(),
	                                           [](string_t ncbi_rank) -> optional<string_t> {
		                                           const auto *column = miint::FindRankColumn(ncbi_rank.GetString());
		                                           if (!column) {
			                                           return std::nullopt;
		                                           }
		                                           // A string literal: it outlives the result.
		                                           return string_t(column->column);
	                                           });
}

void LineageRankFunction::Register(ExtensionLoader &loader) {
	ScalarFunction fn("lineage_rank", {FunctionParameter(Identifier("ncbi_rank"), LogicalType::VARCHAR)},
	                  LogicalType::VARCHAR, LineageRankExecute);
	loader.RegisterFunction(fn);
}

} // namespace duckdb
