#include "read_gtdb.hpp"

#include "GtdbFiles.hpp"
#include "catalog_utils.hpp"
#include "miint_named_parameter.hpp"
#include "miint_streaming_query.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <optional>
#include <stdexcept>

namespace duckdb {

namespace {

enum class GtdbTable { METADATA, QC_FAILED };

const char *FunctionName(GtdbTable table) {
	return table == GtdbTable::METADATA ? "read_gtdb_metadata" : "read_gtdb_qc_failed";
}

// One file the reader reads: its path, whether its bytes are gzip, and for metadata, the
// domain its name implies.
struct GtdbInput {
	string path;
	string domain;
	bool gzip = false;
};

struct ReadGtdbData : public TableFunctionData {
	GtdbTable table = GtdbTable::METADATA;
	vector<GtdbInput> inputs;
	// Every input's header line, which binding found equal, and its column names.
	string header_line;
	vector<string> columns;
};

struct ReadGtdbGlobalState : public GlobalTableFunctionState {
	unique_ptr<Connection> conn;
	unique_ptr<StreamingQuery> stream;
	unique_ptr<DataChunk> chunk;

	idx_t MaxThreads() const override {
		return 1;
	}
};

struct ReadGtdbLocalState : public LocalTableFunctionState {};

// Runs one of GtdbFiles's rules, which report a file that does not fit as
// std::runtime_error, as a DuckDB error naming the function and `where`. A DuckDB
// exception, itself a std::runtime_error, passes through unchanged.
template <class FN>
auto ApplyRule(const char *fn_name, const string &where, FN &&fn) -> decltype(fn()) {
	try {
		return fn();
	} catch (const Exception &) {
		throw;
	} catch (const std::runtime_error &e) {
		throw InvalidInputException("%s: %s: %s", fn_name, where, e.what());
	}
}

std::optional<string> OptionalString(TableFunctionBindInput &input, const char *name) {
	auto it = input.named_parameters.find(Identifier(name));
	if (it == input.named_parameters.end() || it->second.IsNull()) {
		return std::nullopt;
	}
	return it->second.ToString();
}

bool IsHttpUrl(const string &path) {
	return StringUtil::StartsWith(path, "http://") || StringUtil::StartsWith(path, "https://");
}

// The first line of a file and whether it is gzip, reading no further than the line.
struct Header {
	string line;
	bool gzip;
};

Header ReadHeader(ClientContext &context, const char *fn_name, const string &path) {
	if (IsHttpUrl(path)) {
		throw NotImplementedException("%s: reading %s", fn_name, path);
	}
	auto &fs = FileSystem::GetFileSystem(context);
	auto handle = fs.OpenFile(path, FileFlags::FILE_FLAGS_READ);
	miint::gtdb::HeaderLineReader reader;
	vector<char> buffer(65536);
	return ApplyRule(fn_name, path, [&] {
		while (true) {
			const auto n = handle->Read(buffer.data(), buffer.size());
			if (n <= 0 || reader.Feed(buffer.data(), static_cast<size_t>(n))) {
				break;
			}
		}
		auto line = reader.Finish();
		return Header {std::move(line), reader.IsGzip()};
	});
}

// The files a source names: the file itself, or the one release's files in a directory.
vector<GtdbInput> ResolveSource(ClientContext &context, GtdbTable table, const string &source) {
	const auto fn_name = FunctionName(table);
	auto &fs = FileSystem::GetFileSystem(context);
	vector<string> paths;
	if (IsHttpUrl(source)) {
		paths.push_back(source);
	} else if (fs.DirectoryExists(source)) {
		vector<string> names;
		fs.ListFiles(source, [&](const string &name, bool is_dir) {
			if (!is_dir) {
				names.push_back(name);
			}
		});
		const auto matched = ApplyRule(fn_name, source, [&] {
			return table == GtdbTable::METADATA ? miint::gtdb::MatchMetadataFiles(names)
			                                    : vector<string> {miint::gtdb::MatchQcFailedFile(names)};
		});
		for (const auto &name : matched) {
			paths.push_back(fs.JoinPath(source, name));
		}
	} else if (fs.FileExists(source)) {
		paths.push_back(source);
	} else {
		throw IOException("%s: source not found: %s", fn_name, source);
	}

	vector<GtdbInput> inputs;
	for (auto &path : paths) {
		GtdbInput input;
		if (table == GtdbTable::METADATA) {
			const auto name = string(miint::gtdb::FileName(path));
			auto domain = miint::gtdb::MetadataDomain(name);
			if (!domain) {
				throw InvalidInputException("%s: cannot tell the domain of '%s': a metadata file's name starts with "
				                            "bac120_ or ar53_ (ar122_ before release 207)",
				                            fn_name, name);
			}
			input.domain = std::move(*domain);
		}
		input.path = std::move(path);
		inputs.push_back(std::move(input));
	}
	return inputs;
}

unique_ptr<FunctionData> ReadGtdbBind(ClientContext &context, TableFunctionBindInput &input,
                                      vector<LogicalType> &return_types, vector<Identifier> &names, GtdbTable table) {
	const auto fn_name = FunctionName(table);
	const auto source = OptionalString(input, "source");
	const auto release = OptionalString(input, "release");
	if (source.has_value() == release.has_value()) {
		throw InvalidInputException("%s: pass exactly one of source and release", fn_name);
	}

	auto data = make_uniq<ReadGtdbData>();
	data->table = table;
	if (source) {
		if (source->empty()) {
			throw InvalidInputException("%s: source must not be empty", fn_name);
		}
		data->inputs = ResolveSource(context, table, *source);
	} else {
		try {
			miint::gtdb::ParseRelease(*release);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("%s: %s", fn_name, e.what());
		}
		throw NotImplementedException("%s: release", fn_name);
	}

	// The header names the columns. Read together, the two domains' files must agree on
	// them, as one release's do.
	for (auto &in : data->inputs) {
		auto header = ReadHeader(context, fn_name, in.path);
		in.gzip = header.gzip;
		auto columns = ApplyRule(fn_name, in.path, [&] { return miint::gtdb::SplitHeader(header.line); });
		if (&in == &data->inputs.front()) {
			data->header_line = std::move(header.line);
			data->columns = std::move(columns);
			continue;
		}
		if (header.line == data->header_line) {
			continue;
		}
		const auto first = string(miint::gtdb::FileName(data->inputs.front().path));
		const auto other = string(miint::gtdb::FileName(in.path));
		if (columns.size() != data->columns.size()) {
			throw InvalidInputException("%s: %s and %s have different columns: %d and %d", fn_name, first, other,
			                            data->columns.size(), columns.size());
		}
		for (idx_t i = 0; i < columns.size(); ++i) {
			if (columns[i] != data->columns[i]) {
				throw InvalidInputException("%s: %s and %s have different columns: column %d is '%s' and '%s'", fn_name,
				                            first, other, i + 1, data->columns[i], columns[i]);
			}
		}
	}

	if (table == GtdbTable::METADATA) {
		for (const auto &column : data->columns) {
			if (StringUtil::CIEquals(column, "domain")) {
				throw InvalidInputException("%s: %s already has a column named '%s', which the reader adds", fn_name,
				                            string(miint::gtdb::FileName(data->inputs.front().path)), column);
			}
		}
		names.emplace_back("domain");
		return_types.push_back(LogicalType::VARCHAR);
	}
	for (const auto &column : data->columns) {
		names.emplace_back(column);
		return_types.push_back(LogicalType::VARCHAR);
	}
	return std::move(data);
}

unique_ptr<FunctionData> ReadGtdbMetadataBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<Identifier> &names) {
	return ReadGtdbBind(context, input, return_types, names, GtdbTable::METADATA);
}

unique_ptr<FunctionData> ReadGtdbQcFailedBind(ClientContext &context, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<Identifier> &names) {
	return ReadGtdbBind(context, input, return_types, names, GtdbTable::QC_FAILED);
}

string Quoted(const string &text) {
	return KeywordHelper::WriteQuoted(text, '\'');
}

// The query that reads the rows: DuckDB's CSV reader on each file, in GTDB's dialect (tab
// separated, a value holding '"' quoted with the quote doubled), with the bound columns,
// all VARCHAR. Metadata adds the domain, and reads GTDB's 'none' as NULL.
string RowsQuery(const ReadGtdbData &data, const vector<GtdbInput> &inputs) {
	string columns;
	for (const auto &column : data.columns) {
		columns += (columns.empty() ? "" : ", ") + Quoted(column) + ": 'VARCHAR'";
	}
	string query;
	for (const auto &in : inputs) {
		query += query.empty() ? "SELECT " : " UNION ALL SELECT ";
		if (data.table == GtdbTable::METADATA) {
			query += Quoted(in.domain) + " AS domain, ";
		}
		query += "* FROM read_csv(" + Quoted(in.path) +
		         ", delim = '\\t', quote = '\"', escape = '\"', header = true, auto_detect = false, columns = {" +
		         columns + "}, compression = " + (in.gzip ? "'gzip'" : "'none'");
		if (data.table == GtdbTable::METADATA) {
			query += ", nullstr = ['none', '']";
		}
		query += ")";
	}
	return query;
}

unique_ptr<GlobalTableFunctionState> ReadGtdbInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<ReadGtdbData>();
	const auto fn_name = FunctionName(data.table);

	// The rows are read under the bound column names, so a file whose header changed since
	// (rewritten under a prepared statement, say) would be read under the wrong ones.
	auto inputs = data.inputs;
	for (auto &in : inputs) {
		auto header = ReadHeader(context, fn_name, in.path);
		if (header.line != data.header_line) {
			throw InvalidInputException(
			    "%s: the header of '%s' changed after the query was bound; run or prepare it again", fn_name, in.path);
		}
		in.gzip = header.gzip;
	}

	auto gstate = make_uniq<ReadGtdbGlobalState>();
	gstate->conn = make_uniq<Connection>(MakeReadOnlyHelperConnection(context));
	gstate->stream = SubmitStream(*gstate->conn, RowsQuery(data, inputs));
	if (gstate->stream->HasError()) {
		gstate->stream->GetErrorObject().Throw(string(fn_name) + ": ");
	}
	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState> ReadGtdbInitLocal(ExecutionContext &, TableFunctionInitInput &,
                                                      GlobalTableFunctionState *) {
	return make_uniq<ReadGtdbLocalState>();
}

void ReadGtdbExecute(ClientContext &, TableFunctionInput &input, DataChunk &output) {
	auto &data = input.bind_data->Cast<ReadGtdbData>();
	auto &gstate = input.global_state->Cast<ReadGtdbGlobalState>();
	do {
		gstate.chunk = gstate.stream->Fetch();
	} while (gstate.chunk && gstate.chunk->size() == 0);
	if (!gstate.chunk) {
		if (gstate.stream->HasError()) {
			gstate.stream->GetErrorObject().Throw(string(FunctionName(data.table)) + ": ");
		}
		output.SetChildCardinality(0);
		return;
	}
	output.Reference(*gstate.chunk);
}

TableFunction MakeFunction(const char *name, table_function_bind_t bind) {
	TableFunction fn(name, {}, ReadGtdbExecute, bind, ReadGtdbInitGlobal, ReadGtdbInitLocal);
	AddNamedParameter(fn, "source", LogicalType::VARCHAR);
	AddNamedParameter(fn, "release", LogicalType::VARCHAR);
	AddNamedParameter(fn, "refresh", LogicalType::BOOLEAN);
	return fn;
}

} // namespace

void ReadGtdbFunctions::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(MakeFunction("read_gtdb_metadata", ReadGtdbMetadataBind));
	loader.RegisterFunction(MakeFunction("read_gtdb_qc_failed", ReadGtdbQcFailedBind));
}

} // namespace duckdb
