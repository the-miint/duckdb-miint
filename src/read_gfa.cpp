#include "read_gfa.hpp"

namespace duckdb {

namespace {

enum Column : idx_t { LINE_NUMBER, RECORD_TYPE, FIELDS, TAGS };

void SetFields(Vector &list_vec, idx_t row, const std::vector<std::string> &fields) {
	auto offset = ListVector::GetListSize(list_vec);
	ListVector::Reserve(list_vec, offset + fields.size());
	auto &child = ListVector::GetEntry(list_vec);
	for (size_t i = 0; i < fields.size(); i++) {
		FlatVector::GetData<string_t>(child)[offset + i] = StringVector::AddString(child, fields[i]);
	}
	auto &entry = FlatVector::GetData<list_entry_t>(list_vec)[row];
	entry.offset = offset;
	entry.length = fields.size();
	ListVector::SetListSize(list_vec, offset + fields.size());
}

} // namespace

unique_ptr<FunctionData> ReadGfaTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<std::string> &names) {
	names = {"line_number", "record_type", "fields", "tags"};
	return_types = {LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)};
	return BindLineScan(context, input, names, return_types);
}

unique_ptr<LocalTableFunctionState> ReadGfaTableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadGfaTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &local = data_p.local_state->Cast<LocalState>();
	const auto &r = local.record;
	ScanLines(
	    context, data_p, output, [&](std::string_view line) { return local.parser.Parse(line, local.record); },
	    [&](idx_t row) {
		    FlatVector::GetData<int64_t>(output.data[LINE_NUMBER])[row] =
		        static_cast<int64_t>(local.file->lines.LineNumber());
		    FlatVector::GetData<string_t>(output.data[RECORD_TYPE])[row] =
		        StringVector::AddString(output.data[RECORD_TYPE], r.record_type);
		    SetFields(output.data[FIELDS], row, r.fields);
		    if (r.has_tags) {
			    SetTagMap(output.data[TAGS], row, r.tags);
		    } else {
			    FlatVector::SetNull(output.data[TAGS], row, true);
		    }
	    });
}

void ReadGfaTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeLineScanFunction("read_gfa", Execute, Bind, InitLocal));
}

} // namespace duckdb
