#include "unifrac_table_functions.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/function/table_function.hpp"

namespace duckdb {

static unique_ptr<FunctionData> UnifracPermanovaBind(ClientContext &, TableFunctionBindInput &, vector<LogicalType> &,
                                                     vector<string> &) {
	throw NotImplementedException("unifrac_permanova is not yet implemented (Phase 5)");
}

static void UnifracPermanovaExecute(ClientContext &, TableFunctionInput &, DataChunk &) {
}

void RegisterUnifracPermanova(ExtensionLoader &loader) {
	TableFunction fn("unifrac_permanova", {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR},
	                 UnifracPermanovaExecute, UnifracPermanovaBind);
	loader.RegisterFunction(fn);
}

} // namespace duckdb
