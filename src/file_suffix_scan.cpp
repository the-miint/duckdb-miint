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

// A reader claims only the compressed forms it opens, so any other falls through to DuckDB's
// own error instead of reaching a reader that would parse its compressed bytes. read_jplace and
// read_gff read zstd because they are macros over read_json and read_csv.
struct FormatReader {
	const char *reader;
	bool reads_gzip;
	bool reads_zstd;
	const char *suffixes[8];
};

constexpr FormatReader FORMAT_READERS[] = {
    {"read_fastx", true, false, {".fa", ".fasta", ".fna", ".faa", ".ffn", ".frn", ".fq", ".fastq"}},
    {"read_alignments", true, false, {".sam"}},
    {"read_alignments", false, false, {".bam"}},
    {"read_newick", true, false, {".nwk", ".newick", ".tre", ".tree"}},
    {"read_paf", true, false, {".paf"}},
    {"read_gfa", true, false, {".gfa"}},
    {"read_jplace", true, true, {".jplace"}},
    {"read_gff", true, true, {".gff", ".gff3"}},
    {"read_mzml", false, false, {".mzml"}},
    {"read_mzxml", false, false, {".mzxml"}},
    {"read_sequences_sff", false, false, {".sff"}},
#ifdef MIINT_HAS_HDF5
    {"read_biom", false, false, {".biom"}},
#endif
};

// The reader for path, or nullptr. The format's suffix is matched in any case (data.mzML is
// read_mzml's), the compression's as written: a reader that decides gzip by name does so with
// IsGzipped, which takes only a lower-case .gz.
const char *ReaderFor(const string &path) {
	const bool gzip = IsGzipped(path);
	const bool zstd = StringUtil::EndsWith(path, ".zst");
	auto name = StringUtil::Lower(path);
	if (gzip || zstd) {
		name.resize(name.rfind('.'));
	}
	for (const auto &format : FORMAT_READERS) {
		for (const char *suffix : format.suffixes) {
			if (suffix && StringUtil::EndsWith(name, suffix)) {
				const bool opens = (!gzip || format.reads_gzip) && (!zstd || format.reads_zstd);
				return opens ? format.reader : nullptr;
			}
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
