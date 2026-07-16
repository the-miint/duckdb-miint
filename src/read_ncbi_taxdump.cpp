#include "read_ncbi_taxdump.hpp"

#include "ensure_httpfs.hpp"
#include "miint_cache.hpp"
#include "remote_file_helper.hpp"
#include "taxdump_archive.hpp"
#include "duckdb/common/exception/http_exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/http_util.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/function_set.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <thread>

namespace duckdb {

namespace {

constexpr size_t READ_BUFFER_SIZE = 65536;

// Read a whole file (local or remote) into a string via the DuckDB FileSystem.
std::string ReadWholeFile(FileSystem &fs, const std::string &path) {
	auto handle = fs.OpenFile(path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));
	std::string content;
	char buf[READ_BUFFER_SIZE];
	while (true) {
		auto n = handle->Read(buf, sizeof(buf));
		if (n <= 0) {
			break;
		}
		content.append(buf, static_cast<size_t>(n));
	}
	return content;
}

// NCBI's canonical taxonomy dump. Used when read_ncbi_taxdump() is called with no
// explicit source — nobody remembers this URL offhand, so it is the default.
constexpr const char *DEFAULT_TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz";

bool IsHttpUrl(const std::string &path) {
	return path.rfind("http://", 0) == 0 || path.rfind("https://", 0) == 0;
}

// Download `url` into memory via DuckDB's core HTTP client (HTTPUtil), the same
// mechanism read_ncbi uses. This deliberately avoids the httpfs filesystem
// extension, which the extension does not auto-load, so read_ncbi_taxdump() works
// out of the box. Retries transient failures with exponential backoff.
std::string HttpGetToString(ClientContext &context, const std::string &url) {
	miint::EnsureHttpfsLoaded(context);
	DatabaseInstance &db = *context.db;
	HTTPHeaders headers(db);
	auto &http_util = HTTPUtil::Get(db);
	auto params = http_util.InitializeParameters(db, url);
	GetRequestInfo get_request(url, headers, *params, nullptr, nullptr);
	get_request.try_request = true;

	constexpr int MAX_RETRIES = 3;
	int retry_delay_ms = 1000;
	for (int attempt = 0; attempt <= MAX_RETRIES; ++attempt) {
		auto response = http_util.Request(get_request);
		if (response->Success()) {
			return response->body;
		}
		bool retryable = !response->HasRequestError() &&
		                 (int(response->status) >= 500 || int(response->status) == 429 || int(response->status) == 408);
		if (attempt < MAX_RETRIES && retryable) {
			std::this_thread::sleep_for(std::chrono::milliseconds(retry_delay_ms));
			retry_delay_ms *= 2;
			continue;
		}
		if (response->HasRequestError()) {
			throw IOException("read_ncbi_taxdump: download failed: %s (URL: %s)", response->GetRequestError(), url);
		}
		throw HTTPException(*response, "read_ncbi_taxdump: download failed with HTTP %d (URL: %s)",
		                    int(response->status), url);
	}
	throw IOException("read_ncbi_taxdump: download failed after retries (URL: %s)", url);
}

// Read a taxdump archive's raw bytes from `source`: HTTP(S) URLs go through the
// core HTTP client; anything else is read through the DuckDB FileSystem (local
// paths, or other schemes if the user has the relevant extension loaded).
std::string ReadArchiveBytes(ClientContext &context, FileSystem &fs, const std::string &source) {
	if (IsHttpUrl(source)) {
		return HttpGetToString(context, source);
	}
	return ReadWholeFile(fs, source);
}

bool ParseRefreshParam(const named_parameter_map_t &params) {
	auto it = params.find("refresh");
	if (it != params.end() && !it->second.IsNull()) {
		return it->second.GetValue<bool>();
	}
	return false;
}

// Where the auto-downloaded taxdump is cached (extracted). An explicit override
// wins (used by tests to avoid touching the user's real cache); otherwise the
// shared miint cache location.
std::string TaxonomyCacheDir() {
	if (const char *o = ::getenv("MIINT_TAXONOMY_CACHE_DIR"); o && o[0] != '\0') {
		return std::string(o);
	}
	return miint::MiintCacheDir("taxonomy");
}

// Write <dir>/<name> atomically: write a temp file then rename it into place, so a
// concurrent or later reader never observes a half-written cache file.
void WriteCacheFile(FileSystem &fs, const std::string &dir, const char *name, const std::string &content) {
	std::string final_path = fs.JoinPath(dir, name);
	std::string tmp_path = final_path + ".tmp";
	{
		auto handle = fs.OpenFile(
		    tmp_path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_WRITE | FileOpenFlags::FILE_FLAGS_FILE_CREATE_NEW));
		if (!content.empty()) {
			handle->Write(const_cast<char *>(content.data()), static_cast<int64_t>(content.size()));
		}
	}
	fs.MoveFile(tmp_path, final_path);
}

// Ensure the default taxdump is present in the cache and return the cache directory
// (which read_ncbi_taxdump then reads like any other extracted-.dmp directory).
std::string EnsureTaxdumpCache(ClientContext &context, bool refresh) {
	FileSystem &fs = FileSystem::GetFileSystem(context);
	std::string dir = TaxonomyCacheDir();
	if (dir.empty()) {
		throw IOException("read_ncbi_taxdump: cannot determine a cache directory; set HOME, XDG_CACHE_HOME, or "
		                  "MIINT_TAXONOMY_CACHE_DIR, or pass an explicit taxdump path/URL");
	}
	std::string nodes = fs.JoinPath(dir, "nodes.dmp");
	std::string names = fs.JoinPath(dir, "names.dmp");
	if (!refresh && fs.FileExists(nodes) && fs.FileExists(names)) {
		return dir; // cache hit
	}

	fs.CreateDirectoriesRecursive(dir);
	std::string archive_bytes = HttpGetToString(context, DEFAULT_TAXDUMP_URL);
	miint::TaxdumpFiles files;
	try {
		files = miint::TaxdumpArchive::ExtractTaxdump(archive_bytes);
	} catch (const std::exception &e) {
		throw IOException("read_ncbi_taxdump: failed to extract taxdump from " + std::string(DEFAULT_TAXDUMP_URL) +
		                  ": " + e.what());
	}
	WriteCacheFile(fs, dir, "nodes.dmp", files.nodes);
	WriteCacheFile(fs, dir, "names.dmp", files.names);
	WriteCacheFile(fs, dir, "merged.dmp", files.merged);
	WriteCacheFile(fs, dir, "delnodes.dmp", files.delnodes);
	return dir;
}

} // namespace

miint::TaxdumpFiles LoadTaxdumpFiles(ClientContext &context, const std::string &source) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	// A directory of already-extracted .dmp files.
	if (fs.DirectoryExists(source)) {
		auto read_member = [&](const char *name, bool required) -> std::string {
			std::string path = fs.JoinPath(source, name);
			if (!fs.FileExists(path)) {
				if (required) {
					throw IOException("read_ncbi_taxdump: required file not found: " + path);
				}
				return {};
			}
			return ReadWholeFile(fs, path);
		};
		miint::TaxdumpFiles files;
		files.nodes = read_member("nodes.dmp", true);
		files.names = read_member("names.dmp", true);
		files.merged = read_member("merged.dmp", false);
		files.delnodes = read_member("delnodes.dmp", false);
		return files;
	}

	// Otherwise treat `source` as a .tar.gz archive (local file or http(s) URL).
	if (!miint::RemoteFileHelper::IsRemotePath(source) && !fs.FileExists(source)) {
		throw IOException("read_ncbi_taxdump: source not found: " + source);
	}
	std::string archive_bytes = ReadArchiveBytes(context, fs, source);
	try {
		return miint::TaxdumpArchive::ExtractTaxdump(archive_bytes);
	} catch (const std::exception &e) {
		throw IOException("read_ncbi_taxdump: failed to read taxdump archive '" + source + "': " + e.what());
	}
}

// ── read_ncbi_taxdump ───────────────────────────────────────────────────────

unique_ptr<FunctionData> ReadNCBITaxdumpTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                            vector<LogicalType> &return_types,
                                                            vector<std::string> &names) {
	auto data = make_uniq<Data>();
	// No positional arg, or an explicit NULL, => auto-download the default taxdump
	// into the cache. A NULL source is the taxonomy_lineage(source := NULL) default,
	// which must route to the cached-download path (not error) so the macro can
	// delegate here uniformly.
	if (!input.inputs.empty() && !input.inputs[0].IsNull()) {
		data->source = input.inputs[0].ToString();
		if (data->source.empty()) {
			throw InvalidInputException("read_ncbi_taxdump: source must not be empty");
		}
	}
	data->refresh = ParseRefreshParam(input.named_parameters);

	names = {"node_index", "parent_index", "name", "rank", "is_tip"};
	return_types = {LogicalType::BIGINT, LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR,
	                LogicalType::BOOLEAN};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadNCBITaxdumpTableFunction::InitGlobal(ClientContext &context,
                                                                              TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	std::string source = data.source.empty() ? EnsureTaxdumpCache(context, data.refresh) : data.source;
	auto files = LoadTaxdumpFiles(context, source);
	gstate->nodes = miint::TaxdumpParser::ParseNodes(files.nodes, files.names);
	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState>
ReadNCBITaxdumpTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void ReadNCBITaxdumpTableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, gstate.nodes.size() - gstate.cursor);

	auto node_index = FlatVector::GetData<int64_t>(output.data[0]);
	auto &parent_vec = output.data[1];
	auto parent_index = FlatVector::GetData<int64_t>(parent_vec);
	auto &name_vec = output.data[2];
	auto &rank_vec = output.data[3];
	auto is_tip = FlatVector::GetData<bool>(output.data[4]);

	for (idx_t i = 0; i < count; ++i) {
		const auto &node = gstate.nodes[gstate.cursor + i];
		node_index[i] = node.taxid;
		if (node.parent_taxid.has_value()) {
			parent_index[i] = node.parent_taxid.value();
		} else {
			FlatVector::Validity(parent_vec).SetInvalid(i);
		}
		FlatVector::GetData<string_t>(name_vec)[i] = StringVector::AddString(name_vec, node.name);
		FlatVector::GetData<string_t>(rank_vec)[i] = StringVector::AddString(rank_vec, node.rank);
		is_tip[i] = node.is_tip;
	}

	gstate.cursor += count;
	output.SetCardinality(count);
}

void ReadNCBITaxdumpTableFunction::Register(ExtensionLoader &loader) {
	TableFunctionSet set("read_ncbi_taxdump");
	// read_ncbi_taxdump() -> default remote taxdump; read_ncbi_taxdump('path'|'url') -> explicit source.
	for (const auto &args : {vector<LogicalType> {}, vector<LogicalType> {LogicalType::VARCHAR}}) {
		TableFunction tf("read_ncbi_taxdump", args, Execute, Bind, InitGlobal, InitLocal);
		tf.named_parameters["refresh"] = LogicalType::BOOLEAN;
		set.AddFunction(tf);
	}
	loader.RegisterFunction(set);
}

// ── read_ncbi_taxdump_merged ─────────────────────────────────────────────────

unique_ptr<FunctionData> ReadNCBITaxdumpMergedTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                                  vector<LogicalType> &return_types,
                                                                  vector<std::string> &names) {
	auto data = make_uniq<Data>();
	// No positional arg, or an explicit NULL, => auto-download the default taxdump
	// (mirrors read_ncbi_taxdump so both accept the NULL-source default).
	if (!input.inputs.empty() && !input.inputs[0].IsNull()) {
		data->source = input.inputs[0].ToString();
		if (data->source.empty()) {
			throw InvalidInputException("read_ncbi_taxdump_merged: source must not be empty");
		}
	}
	data->refresh = ParseRefreshParam(input.named_parameters);

	names = {"old_taxid", "new_taxid"};
	return_types = {LogicalType::BIGINT, LogicalType::BIGINT};
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadNCBITaxdumpMergedTableFunction::InitGlobal(ClientContext &context,
                                                                                    TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();
	std::string source = data.source.empty() ? EnsureTaxdumpCache(context, data.refresh) : data.source;
	auto files = LoadTaxdumpFiles(context, source);
	gstate->merged = miint::TaxdumpParser::ParseMerged(files.merged);
	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState> ReadNCBITaxdumpMergedTableFunction::InitLocal(ExecutionContext &,
                                                                                  TableFunctionInitInput &,
                                                                                  GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void ReadNCBITaxdumpMergedTableFunction::Execute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	idx_t count = std::min<idx_t>(STANDARD_VECTOR_SIZE, gstate.merged.size() - gstate.cursor);

	auto old_taxid = FlatVector::GetData<int64_t>(output.data[0]);
	auto new_taxid = FlatVector::GetData<int64_t>(output.data[1]);
	for (idx_t i = 0; i < count; ++i) {
		const auto &m = gstate.merged[gstate.cursor + i];
		old_taxid[i] = m.old_taxid;
		new_taxid[i] = m.new_taxid;
	}

	gstate.cursor += count;
	output.SetCardinality(count);
}

void ReadNCBITaxdumpMergedTableFunction::Register(ExtensionLoader &loader) {
	TableFunctionSet set("read_ncbi_taxdump_merged");
	for (const auto &args : {vector<LogicalType> {}, vector<LogicalType> {LogicalType::VARCHAR}}) {
		TableFunction tf("read_ncbi_taxdump_merged", args, Execute, Bind, InitGlobal, InitLocal);
		tf.named_parameters["refresh"] = LogicalType::BOOLEAN;
		set.AddFunction(tf);
	}
	loader.RegisterFunction(set);
}

} // namespace duckdb
