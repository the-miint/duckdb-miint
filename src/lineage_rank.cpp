#include "lineage_rank.hpp"

#include "taxonomy_lineage.hpp"

#include "duckdb/common/optional.hpp"
#include "duckdb/common/vector_operations/unary_executor.hpp"
#include "duckdb/function/scalar_function.hpp"

#include <string_view>

namespace duckdb {

// The mapping is FindRankColumn's, the one read_ncbi_lineage collapses with; the
// taxonomy_lineage macro calls this function, so none of them maps a rank differently.
static optional<string_t> RankColumn(string_t ncbi_rank) {
	const auto *column = miint::FindRankColumn(std::string_view(ncbi_rank.GetData(), ncbi_rank.GetSize()));
	if (!column) {
		return std::nullopt;
	}
	// A string literal: it outlives the result.
	return string_t(column->column);
}

static void LineageRankExecute(DataChunk &args, ExpressionState &, Vector &result) {
	UnaryExecutor::Execute<string_t, string_t>(args.data[0], result, args.size(), RankColumn);
}

void LineageRankFunction::Register(ExtensionLoader &loader) {
	ScalarFunction fn("lineage_rank", {FunctionParameter(Identifier("ncbi_rank"), LogicalType::VARCHAR)},
	                  LogicalType::VARCHAR, LineageRankExecute);
	loader.RegisterFunction(fn);
}

} // namespace duckdb
