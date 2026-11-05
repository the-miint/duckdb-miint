#pragma once

#include "sra_source.hpp"

#include <cstdint>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace miint {
namespace sra {

// One entry of a KAR archive's table of contents. `path` is the full path inside the archive
// ("tbl/SEQUENCE/col/READ/data"). For a file, `offset` is absolute in the .sra file.
struct KarEntry {
	enum class Type { Directory, File, EmptyFile };
	std::string path;
	Type type;
	uint64_t offset = 0;
	uint64_t size = 0;
};

// The KAR container that wraps every .sra: a 24-byte header, a persisted tree of entries, then
// the file bodies. Only the entry types seen in public runs are read (directory, file, empty
// file); chunked files and links raise SraUnsupportedError, as does a big-endian archive.
class KarArchive {
public:
	explicit KarArchive(std::shared_ptr<const SraByteSource> source);

	const SraByteSource &Source() const {
		return *source_;
	}
	// Absolute offset where the file bodies start; the TOC fills [24, DataOffset()).
	uint64_t DataOffset() const {
		return data_offset_;
	}
	const std::vector<KarEntry> &Entries() const {
		return entries_;
	}
	// nullptr if there is no entry at that path.
	const KarEntry *Find(const std::string &path) const;
	bool IsDirectory(const std::string &path) const;
	// Names (not paths) of a directory's direct children, sorted. Empty if not a directory.
	std::vector<std::string> Children(const std::string &dir) const;
	// The whole file, or n bytes at off within it. Throws SraFormatError if the path is not a
	// file or the range runs past it.
	std::vector<uint8_t> ReadFile(const std::string &path) const;
	std::vector<uint8_t> ReadFile(const std::string &path, uint64_t off, uint64_t n) const;

private:
	std::shared_ptr<const SraByteSource> source_;
	uint64_t data_offset_ = 0;
	std::vector<KarEntry> entries_;
	std::map<std::string, size_t> by_path_;
};

// A node of a KDB metadata tree (an `md/cur` file): a name, string attributes, children, and a
// value of raw bytes.
struct MetadataNode {
	std::string name;
	std::map<std::string, std::string> attrs;
	std::vector<MetadataNode> children;
	std::string value;

	// Descends by '/'-separated child names; nullptr if any step is missing.
	const MetadataNode *Find(const std::string &path) const;
	// The value as an unsigned integer stored in 1, 2, 4 or 8 bytes, as KDB writes counts.
	uint64_t ValueAsU64() const;
	int64_t ValueAsI64() const;
};

// Parses an `md/cur` file. The returned root has no name; its children are the top-level nodes.
MetadataNode ParseMetadata(const std::vector<uint8_t> &bytes);

// A column whose value is the same in every row. VDB stores it in the table's metadata under
// col/<NAME> instead of as a column directory.
struct StaticColumn {
	std::string type; // the column's stored type, e.g. "INSDC:coord:len"
	std::string row;  // the cell's bytes, the same for every row
	int64_t start_id = 1;
	uint64_t row_count = 0;
	// The cell's length in bits, when its bytes do not say it: a 2na cell of 574 bases fills
	// 144 bytes but holds 1148 bits.
	std::optional<uint64_t> size_bits;
};

// The table holding the reads: tbl/SEQUENCE inside a database, or the archive itself when the
// run is a flat table (no database).
class SraTable {
public:
	explicit SraTable(std::shared_ptr<const KarArchive> archive);

	const KarArchive &Archive() const {
		return *archive_;
	}
	bool InDatabase() const {
		return in_database_;
	}
	// Tables of the database (empty for a flat table).
	const std::vector<std::string> &DatabaseTables() const {
		return database_tables_;
	}
	// Where the table's files live in the archive: "tbl/SEQUENCE/" or "".
	const std::string &Prefix() const {
		return prefix_;
	}
	const MetadataNode &Metadata() const {
		return metadata_;
	}
	// The database's metadata; nullptr for a flat table.
	const MetadataNode *DatabaseMetadata() const {
		return db_metadata_ ? &*db_metadata_ : nullptr;
	}
	// The schema type the table was written with, e.g. "NCBI:SRA:Illumina:tbl:phred:v2#2".
	std::string TableType() const;
	// Columns stored as column directories, sorted.
	const std::vector<std::string> &PhysicalColumns() const {
		return physical_columns_;
	}
	bool HasPhysical(const std::string &name) const;
	// nullptr if the column is not static.
	const StaticColumn *Static(const std::string &name) const;
	const std::map<std::string, StaticColumn> &StaticColumns() const {
		return static_columns_;
	}

private:
	std::shared_ptr<const KarArchive> archive_;
	bool in_database_ = false;
	std::vector<std::string> database_tables_;
	std::string prefix_;
	MetadataNode metadata_;
	std::optional<MetadataNode> db_metadata_;
	std::vector<std::string> physical_columns_;
	std::map<std::string, StaticColumn> static_columns_;
};

} // namespace sra
} // namespace miint
