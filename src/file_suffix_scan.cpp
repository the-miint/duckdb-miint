#include "file_suffix_scan.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/file_system.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/function/replacement_scan.hpp"
#include "duckdb/main/config.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/parser/expression/constant_expression.hpp"
#include "duckdb/parser/expression/function_expression.hpp"
#include "duckdb/parser/tableref/table_function_ref.hpp"

namespace duckdb {

namespace {

// A suffix claims only the compressed forms its reader opens, so a file the reader cannot read
// is reported by DuckDB as one no reader handles, not as a parse error from the reader.
// read_jplace and read_gff read zstd because they are macros over read_json and read_csv.
struct SuffixReader {
	const char *suffix;
	const char *reader;
	bool reads_gzip;
	bool reads_zstd;
};

constexpr SuffixReader SUFFIX_READERS[] = {
    {".fa", "read_fastx", true, false},       {".fasta", "read_fastx", true, false},
    {".fna", "read_fastx", true, false},      {".faa", "read_fastx", true, false},
    {".ffn", "read_fastx", true, false},      {".frn", "read_fastx", true, false},
    {".fq", "read_fastx", true, false},       {".fastq", "read_fastx", true, false},
    {".sam", "read_alignments", true, false}, {".bam", "read_alignments", false, false},
    {".nwk", "read_newick", true, false},     {".newick", "read_newick", true, false},
    {".tre", "read_newick", true, false},     {".tree", "read_newick", true, false},
    {".paf", "read_paf", true, false},        {".gfa", "read_gfa", true, false},
    {".jplace", "read_jplace", true, true},   {".gff", "read_gff", true, true},
    {".gff3", "read_gff", true, true},        {".mzml", "read_mzml", false, false},
    {".mzxml", "read_mzxml", false, false},   {".sff", "read_sequences_sff", false, false},
#ifdef MIINT_HAS_HDF5
    {".biom", "read_biom", false, false},
#endif
};

// The reader for path's suffix, or nullptr. The format's suffix is matched in any case
// (data.mzML is read_mzml's), the compression's as written: read_newick, read_paf and
// read_gfa take a file for gzip only by a lower-case .gz.
const char *ReaderFor(const string &path) {
	const bool gzip = IsGzipped(path);
	const bool zstd = StringUtil::EndsWith(path, ".zst");
	auto name = StringUtil::Lower(path);
	if (gzip) {
		name.resize(name.size() - 3);
	} else if (zstd) {
		name.resize(name.size() - 4);
	}
	for (const auto &entry : SUFFIX_READERS) {
		if (StringUtil::EndsWith(name, entry.suffix)) {
			if ((gzip && !entry.reads_gzip) || (zstd && !entry.reads_zstd)) {
				return nullptr;
			}
			return entry.reader;
		}
	}
	return nullptr;
}

// As DuckDB's parquet and CSV scans: the whole qualified name is the path, so FROM x.fq.gz
// (catalog x, schema fq, table gz) reads x.fq.gz, and the table is named after the file.
unique_ptr<TableRef> FileSuffixScan(ClientContext &context, ReplacementScanInput &input,
                                    optional_ptr<ReplacementScanData> data) {
	auto path = ReplacementScan::GetFullPath(input);
	auto reader = ReaderFor(path);
	if (!reader) {
		return nullptr;
	}
	auto table_function = make_uniq<TableFunctionRef>();
	vector<unique_ptr<ParsedExpression>> children;
	children.push_back(ConstantExpression::String(path));
	table_function->function = make_uniq<FunctionExpression>(Identifier(reader), std::move(children));

	if (!FileSystem::HasGlob(path)) {
		auto &fs = FileSystem::GetFileSystem(context);
		table_function->alias = Identifier(fs.ExtractBaseName(path));
	}
	return std::move(table_function);
}

} // namespace

void RegisterFileSuffixScan(ExtensionLoader &loader) {
	DBConfig::GetConfig(loader.GetDatabaseInstance()).replacement_scans.emplace_back(FileSuffixScan);
}

} // namespace duckdb
