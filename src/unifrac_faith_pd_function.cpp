#include "unifrac_table_functions.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/function/table_function.hpp"

namespace duckdb {

static unique_ptr<FunctionData> UnifracFaithPDBind(ClientContext &, TableFunctionBindInput &, vector<LogicalType> &,
                                                   vector<string> &) {
	throw NotImplementedException("unifrac_faith_pd is not yet implemented (Phase 6)");
}

static void UnifracFaithPDExecute(ClientContext &, TableFunctionInput &, DataChunk &) {
}

void RegisterUnifracFaithPD(ExtensionLoader &loader) {
	TableFunction fn("unifrac_faith_pd", {LogicalType::VARCHAR, LogicalType::VARCHAR}, UnifracFaithPDExecute,
	                 UnifracFaithPDBind);
	loader.RegisterFunction(fn);
}

} // namespace duckdb
