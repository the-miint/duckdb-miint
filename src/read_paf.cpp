#include "read_paf.hpp"

namespace duckdb {

namespace {

enum Column : idx_t {
	READ_ID,
	QUERY_LENGTH,
	QUERY_START,
	QUERY_STOP,
	STRAND,
	REFERENCE,
	REFERENCE_LENGTH,
	POSITION,
	STOP_POSITION,
	RESIDUE_MATCHES,
	BLOCK_LENGTH,
	MAPQ,
	TAGS
};

void SetString(Vector &vec, idx_t row, const std::string &value) {
	FlatVector::GetData<string_t>(vec)[row] = StringVector::AddString(vec, value);
}

void SetInt64(Vector &vec, idx_t row, int64_t value) {
	FlatVector::GetData<int64_t>(vec)[row] = value;
}

} // namespace

unique_ptr<FunctionData> ReadPafTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<std::string> &names) {
	names = {"read_id",  "query_length",  "query_start",     "query_stop",   "strand", "reference", "reference_length",
	         "position", "stop_position", "residue_matches", "block_length", "mapq",   "tags"};
	return_types = {LogicalType::VARCHAR,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::VARCHAR,
	                LogicalType::VARCHAR,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::BIGINT,
	                LogicalType::UTINYINT,
	                LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)};
	return BindLineScan(context, input, names, return_types);
}

unique_ptr<LocalTableFunctionState> ReadPafTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadPafTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &local = data_p.local_state->Cast<LocalState>();
	const auto &r = local.record;
	ScanLines(
	    context, data_p, output,
	    [&](std::string_view line) {
		    local.parser.Parse(line, local.record);
		    return true;
	    },
	    [&](idx_t row) {
		    SetString(output.data[READ_ID], row, r.read_id);
		    SetInt64(output.data[QUERY_LENGTH], row, r.query_length);
		    SetInt64(output.data[QUERY_START], row, r.query_start);
		    SetInt64(output.data[QUERY_STOP], row, r.query_stop);
		    SetString(output.data[STRAND], row, r.strand);
		    SetString(output.data[REFERENCE], row, r.reference);
		    SetInt64(output.data[REFERENCE_LENGTH], row, r.reference_length);
		    SetInt64(output.data[POSITION], row, r.position);
		    SetInt64(output.data[STOP_POSITION], row, r.stop_position);
		    SetInt64(output.data[RESIDUE_MATCHES], row, r.residue_matches);
		    SetInt64(output.data[BLOCK_LENGTH], row, r.block_length);
		    FlatVector::GetData<uint8_t>(output.data[MAPQ])[row] = r.mapq;
		    SetTagMap(output.data[TAGS], row, r.tags);
	    });
}

void ReadPafTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeLineScanFunction("read_paf", Execute, Bind, InitLocal));
}

} // namespace duckdb
