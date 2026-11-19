#include "read_hmm.hpp"

namespace duckdb {

unique_ptr<FunctionData> ReadHmmTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<Identifier> &names) {
	throw InvalidInputException("read_hmm: not implemented");
}

unique_ptr<LocalTableFunctionState> ReadHmmTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LineScanLocalState>();
}

void ReadHmmTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	output.SetChildCardinality(0);
}

void ReadHmmTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeLineScanFunction("read_hmm", Execute, Bind, InitLocal));
}

} // namespace duckdb
