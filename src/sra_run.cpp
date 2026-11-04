#include "sra_run.hpp"

#include <algorithm>

namespace miint {
namespace sra {

namespace {

// The columns the reader takes reads, qualities, read boundaries and spot names from. Only their
// encodings must be supported; a run may store other columns (454 signal, say) any way it likes.
const std::set<std::string> &ColumnsRead() {
	static const std::set<std::string> COLUMNS = {"READ",     "ALTREAD",   "QUALITY", "ORIGINAL_QUALITY", "READ_START",
	                                              "READ_LEN", "READ_TYPE", "NAME",    "NAME_FMT",         "SPOT_NAME",
	                                              "RAW_NAME", "X",         "Y"};
	return COLUMNS;
}

// Sequence-table types whose decoding has been checked against fasterq-dump; the version after
// '#' may vary.
bool SupportedTableType(const std::string &type) {
	static const std::set<std::string> TYPES = {"NCBI:SRA:Illumina:tbl:phred:v2", "NCBI:SRA:GenericFastq:sequence",
	                                            "NCBI:SRA:GenericFastq:sequence_nanopore", "NCBI:align:tbl:seq"};
	return TYPES.count(type.substr(0, type.find('#'))) > 0;
}

// Raw, zip_encoding#1 (with optional deflate parameters) and izip_encoding#1: the only encodings
// the surveyed runs used on the columns that are read. The element type in angle brackets and
// any spacing are not part of the encoding.
bool SupportedEncoding(const std::string &expr) {
	std::string e;
	for (char c : expr) {
		if (c != ' ') {
			e += c;
		}
	}
	if (!e.empty() && e[0] == '<') {
		auto close = e.find('>');
		if (close == std::string::npos) {
			return false;
		}
		e = e.substr(close + 1);
	}
	if (e.empty() || e == "izip_encoding#1") {
		return true;
	}
	const std::string zip = "zip_encoding#1";
	if (e.compare(0, zip.size(), zip) != 0) {
		return false;
	}
	std::string params = e.substr(zip.size());
	return params.empty() || (params.front() == '<' && params.back() == '>' && params.find('>') == params.size() - 1);
}

bool Has(const std::vector<std::string> &v, const std::string &x) {
	return std::find(v.begin(), v.end(), x) != v.end();
}

} // namespace

RunFacts GatherRunFacts(const SraTable &table) {
	RunFacts f;
	f.in_database = table.InDatabase();
	f.database_tables = table.DatabaseTables();
	f.table_type = table.TableType();
	f.delite = table.Metadata().Find("SOFTWARE/delite") != nullptr ||
	           (table.DatabaseMetadata() && table.DatabaseMetadata()->Find("SOFTWARE/delite") != nullptr);
	for (const auto &name : table.PhysicalColumns()) {
		f.physical_columns.insert(name);
		std::string md_path = table.Prefix() + "col/" + name + "/md/cur";
		std::string expr = "(no schema)";
		if (table.Archive().Find(md_path) != nullptr) {
			MetadataNode md = ParseMetadata(table.Archive().ReadFile(md_path));
			if (const MetadataNode *schema = md.Find("schema")) {
				auto it = schema->attrs.find("expr");
				expr = it == schema->attrs.end() ? std::string() : it->second;
			}
		}
		f.encodings[name] = expr;
	}
	for (const auto &kv : table.StaticColumns()) {
		f.static_columns.insert(kv.first);
	}
	if (const StaticColumn *pai = table.Static("PRIMARY_ALIGNMENT_ID")) {
		f.static_primary_alignment_id = pai->row;
	}
	return f;
}

void ClassifyRun(const RunFacts &f) {
	// Aligned data (cSRA): reads stored against a reference must be restored from alignments.
	for (const char *t : {"PRIMARY_ALIGNMENT", "SECONDARY_ALIGNMENT", "REFERENCE"}) {
		if (Has(f.database_tables, t)) {
			throw SraUnsupportedError(std::string("cSRA (aligned reads, with a ") + t + " table) is not supported");
		}
	}
	if (f.physical_columns.count("CMP_READ") || f.physical_columns.count("PRIMARY_ALIGNMENT_ID") ||
	    f.static_primary_alignment_id.find_first_not_of('\0') != std::string::npos) {
		throw SraUnsupportedError("cSRA (reads stored against alignments) is not supported");
	}
	// Native PacBio and Nanopore consensus databases keep their reads elsewhere than SEQUENCE.
	for (const char *t : {"CONSENSUS", "ZMW_METRICS", "PASSES"}) {
		if (Has(f.database_tables, t)) {
			throw SraUnsupportedError(std::string("a PacBio or Nanopore database with a ") + t +
			                          " table (consensus or polymerase reads) is not supported");
		}
	}
	if (f.in_database && !Has(f.database_tables, "SEQUENCE")) {
		throw SraUnsupportedError("the database has no SEQUENCE table");
	}
	// SRA Lite replaces original qualities with synthesized ones.
	auto has_column = [&](const char *c) {
		return f.physical_columns.count(c) || f.static_columns.count(c);
	};
	if (f.delite || (!has_column("QUALITY") && !has_column("ORIGINAL_QUALITY"))) {
		throw SraUnsupportedError("SRA Lite (no original qualities) is not supported");
	}
	if (!SupportedTableType(f.table_type)) {
		throw SraUnsupportedError("table type '" + f.table_type + "' is not supported");
	}
	for (const auto &kv : f.encodings) {
		if (ColumnsRead().count(kv.first) && !SupportedEncoding(kv.second)) {
			throw SraUnsupportedError("column '" + kv.first + "' uses encoding '" + kv.second +
			                          "', which is not supported");
		}
	}
}

} // namespace sra
} // namespace miint
