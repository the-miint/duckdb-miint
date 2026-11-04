#pragma once

#include "sra_archive.hpp"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace miint {
namespace sra {

// What classification looks at, gathered from an opened table. Kept separate from SraTable so a
// test can describe a run no committed fixture is (a native PacBio database, say) without
// building an archive.
struct RunFacts {
	bool in_database = false;
	std::vector<std::string> database_tables;
	std::string table_type;
	// SOFTWARE/delite in the table's or the database's metadata: the run was converted to SRA
	// Lite and its qualities replaced.
	bool delite = false;
	std::set<std::string> physical_columns;
	std::set<std::string> static_columns;
	// Encoding expression of each physical column (empty when unencoded).
	std::map<std::string, std::string> encodings;
	// Bytes of the static PRIMARY_ALIGNMENT_ID cell, if that column is static.
	std::string static_primary_alignment_id;
};

RunFacts GatherRunFacts(const SraTable &table);

// Throws SraUnsupportedError, naming what the run is, unless it is a kind of run miint reads:
// unaligned reads in a SEQUENCE table of a supported type, with original qualities, stored with
// supported encodings. The decision follows the tables present, as fasterq-dump's does, not the
// platform: PacBio and Nanopore runs loaded as ordinary sequence tables are read.
void ClassifyRun(const RunFacts &facts);

} // namespace sra
} // namespace miint
