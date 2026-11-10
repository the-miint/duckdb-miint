#include "read_gtdb.hpp"

#include "GtdbFiles.hpp"
#include "catalog_utils.hpp"
#include "http_client_common.hpp"
#include "miint_cache.hpp"
#include "miint_named_parameter.hpp"
#include "miint_streaming_query.hpp"
#include "duckdb/common/error_data.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/common/sql_identifier.hpp"

#include <optional>
#include <stdexcept>

namespace duckdb {

namespace {

enum class GtdbTable { METADATA, QC_FAILED };

const char *FunctionName(GtdbTable table) {
	return table == GtdbTable::METADATA ? "read_gtdb_metadata" : "read_gtdb_qc_failed";
}

// One file the reader reads: where the rows are read from (a local file, a URL, or the
// cache's copy of a release's file), for a release's file the URL it is downloaded from,
// for metadata the domain its name implies, and whether its bytes are gzip, which is
// read from them when the query runs.
struct GtdbInput {
	string path;
	string url;
	string domain;
	bool gzip = false;
};

struct ReadGtdbData : public TableFunctionData {
	GtdbTable table = GtdbTable::METADATA;
	vector<GtdbInput> inputs;
	// For a release: its name ("232.0"), and whether to download its files though cached.
	string release;
	bool refresh = false;
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

// The first line of a file and whether it is gzip, reading no further than the line.
struct Header {
	string line;
	bool gzip;
};

// A URL's header is read by a GET stopped once the line is in, so binding does not fetch
// the file (GTDB's bacterial metadata is hundreds of MB).
Header ReadUrlHeader(ClientContext &context, const char *fn_name, const string &url) {
	miint::gtdb::HeaderLineReader reader;
	return ApplyRule(fn_name, url, [&] {
		miint::HttpGetStreaming(
		    context, url, fn_name, [&](const HTTPResponse &) { reader = miint::gtdb::HeaderLineReader(); },
		    [&](const char *data, size_t size) { return !reader.Feed(data, size); });
		auto line = reader.Finish();
		return Header {std::move(line), reader.IsGzip()};
	});
}

Header ReadHeader(ClientContext &context, const char *fn_name, const string &path) {
	if (miint::IsHttpUrl(path)) {
		return ReadUrlHeader(context, fn_name, path);
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
vector<GtdbInput> ResolveSource(ClientContext &context, GtdbTable table, const string &given) {
	const auto fn_name = FunctionName(table);
	auto &fs = FileSystem::GetFileSystem(context);
	// '~' is expanded here, in the caller's context: the rows are read on a helper
	// connection, which would expand it by its own settings.
	const auto source = miint::IsHttpUrl(given) ? given : fs.ExpandPath(given);
	vector<string> paths;
	if (miint::IsHttpUrl(source)) {
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

// Where downloaded releases are kept: the miint_gtdb_cache_dir setting, or by default
// ${XDG_CACHE_HOME:-~/.cache}/miint/gtdb.
string GtdbCacheDir(ClientContext &context, const char *fn_name) {
	Value value;
	string dir;
	if (context.TryGetCurrentSetting("miint_gtdb_cache_dir", value) && !value.IsNull()) {
		dir = value.ToString();
	}
	if (dir.empty()) {
		dir = miint::MiintCacheDir("gtdb");
	} else {
		dir = FileSystem::GetFileSystem(context).ExpandPath(dir);
	}
	if (dir.empty()) {
		throw InvalidInputException("%s: cannot tell where to cache GTDB releases: set miint_gtdb_cache_dir, or HOME "
		                            "or XDG_CACHE_HOME",
		                            fn_name);
	}
	return dir;
}

// The files a release needs, each read from its cached copy, which is downloaded from
// GTDB (the miint_gtdb_url setting) when absent.
vector<GtdbInput> ResolveRelease(ClientContext &context, GtdbTable table, const miint::gtdb::Release &release) {
	using miint::gtdb::ReleaseFile;
	const auto base_url = miint::UrlSetting(context, "miint_gtdb_url");
	const auto cache_dir = GtdbCacheDir(context, FunctionName(table));
	const auto files = table == GtdbTable::METADATA
	                       ? vector<ReleaseFile> {ReleaseFile::BAC120_METADATA, ReleaseFile::AR53_METADATA}
	                       : vector<ReleaseFile> {ReleaseFile::QC_FAILED};
	vector<GtdbInput> inputs;
	for (const auto file : files) {
		GtdbInput input;
		input.path = miint::gtdb::ReleaseCachePath(cache_dir, file, release);
		input.url = miint::gtdb::ReleaseFileUrl(base_url, file, release);
		if (table == GtdbTable::METADATA) {
			input.domain = miint::gtdb::MetadataDomain(miint::gtdb::ReleaseFileName(file, release)).value();
		}
		inputs.push_back(std::move(input));
	}
	return inputs;
}

// Runs `fn`, which reads a file from `url`; for a release's file, a 404 says why GTDB
// may not have it.
template <class FN>
auto WithReleaseHint(const char *fn_name, const string &release, const string &url, FN &&fn) -> decltype(fn()) {
	try {
		return fn();
	} catch (const HTTPException &e) {
		ErrorData error(e);
		auto status = error.ExtraInfo().find("status_code");
		if (release.empty() || status == error.ExtraInfo().end() || status->second != "404") {
			throw;
		}
		throw IOException("%s: release %s has no file at %s (HTTP 404): releases before 214.1 ship .tar.gz archives; "
		                  "extract the .tsv and pass it as source",
		                  fn_name, release, url);
	}
}

// Whether a release's file must be fetched: it is not cached, or refresh asked for it.
bool NeedsDownload(ClientContext &context, const ReadGtdbData &data, const GtdbInput &input) {
	return !input.url.empty() && (data.refresh || !FileSystem::GetFileSystem(context).FileExists(input.path));
}

// Downloads `url` to `path`. The body streams to a file of its own beside `path`, which is
// renamed to `path` only once it is complete, so a file at its final name is always
// whole, and processes downloading at once each write their own. A failed or interrupted
// download removes its file.
void DownloadToCache(ClientContext &context, const char *fn_name, const string &url, const string &path) {
	auto &fs = FileSystem::GetFileSystem(context);
	fs.CreateDirectoriesRecursive(path.substr(0, path.find_last_of('/')));
	const auto part = path + "." + UUID::ToString(UUID::GenerateRandomUUID()) + ".part";
	unique_ptr<FileHandle> out;
	idx_t written = 0;
	try {
		miint::HttpGetStreaming(
		    context, url, fn_name,
		    // A retried body starts over, so the file does too.
		    [&](const HTTPResponse &) {
			    out.reset();
			    out = fs.OpenFile(part, FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW);
			    written = 0;
		    },
		    [&](const char *data, size_t size) {
			    out->Write(const_cast<char *>(data), static_cast<int64_t>(size));
			    written += size;
			    return true;
		    });
		// No GTDB file is empty, and an empty one cached would be trusted until refresh.
		if (written == 0) {
			throw IOException("%s: %s answered with an empty body", fn_name, url);
		}
		out->Sync();
		out->Close();
		out.reset();
		fs.MoveFile(part, path);
	} catch (...) {
		out.reset();
		fs.TryRemoveFile(part);
		throw;
	}
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
	data->refresh = GetBoolParam(input, "refresh", false);
	if (source) {
		if (source->empty()) {
			throw InvalidInputException("%s: source must not be empty", fn_name);
		}
		if (data->refresh) {
			throw InvalidInputException("%s: refresh applies to release; a source is read as it is, never cached",
			                            fn_name);
		}
		data->inputs = ResolveSource(context, table, *source);
	} else {
		miint::gtdb::Release parsed;
		try {
			parsed = miint::gtdb::ParseRelease(*release);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("%s: %s", fn_name, e.what());
		}
		data->release = parsed.ToString();
		data->inputs = ResolveRelease(context, table, parsed);
	}

	// The header names the columns. Read together, the two domains' files must agree on
	// them, as one release's do. A release's file not yet cached (or to be fetched again)
	// has its header read from GTDB.
	for (auto &in : data->inputs) {
		const auto header_from_url = NeedsDownload(context, *data, in);
		const auto &from = header_from_url ? in.url : in.path;
		auto header =
		    WithReleaseHint(fn_name, data->release, in.url, [&] { return ReadHeader(context, fn_name, from); });
		auto columns = ApplyRule(fn_name, from, [&] { return miint::gtdb::SplitHeader(header.line); });
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
	return SQLString::ToString(text);
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

	// A release's files are downloaded now, when the query runs, never while binding.
	// Then the rows are read under the bound column names, so a file whose header changed
	// since (rewritten under a prepared statement, say) would be read under the wrong ones.
	auto inputs = data.inputs;
	for (auto &in : inputs) {
		if (NeedsDownload(context, data, in)) {
			WithReleaseHint(fn_name, data.release, in.url, [&] { DownloadToCache(context, fn_name, in.url, in.path); });
		}
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
