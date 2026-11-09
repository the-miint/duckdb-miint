#pragma once

#include "sra_archive.hpp"

#include <cstdint>
#include <map>
#include <memory>
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

// One spot (one row of the SEQUENCE table) as fasterq-dump reads it: the bases of all its reads
// end to end, one quality per base, and where each read starts.
struct Spot {
	int64_t row_id = 0;
	std::string bases;            // IUPAC letters, e.g. "ACGTN"
	std::vector<uint8_t> quality; // phred values, one per base
	std::vector<int32_t> read_start;
	std::vector<uint32_t> read_len;
	std::vector<uint8_t> read_type; // empty when the table has no READ_TYPE column
	bool has_read_type = false;
	// The spot's name ($sn); has_name is false when the table has no way to name spots, and
	// fasterq-dump then leaves the name out of its deflines.
	bool has_name = false;
	std::string name;
};

// Reads spots of a classified run. The logical columns are taken the way the run's schema
// produces them: READ from 2na bases ORed with ALTREAD's ambiguity codes, QUALITY from
// ORIGINAL_QUALITY or QUALITY, and each column from its directory or, when its value never
// changes, from the table's metadata. Not thread-safe; use one reader per thread.
class SpotReader {
public:
	// A reader of the same run as `share_from` takes its name index rather than reading another
	// copy: each decoding thread has a reader, and the index of a run that stores whole names
	// holds tens of bytes per spot.
	explicit SpotReader(std::shared_ptr<const SraTable> table, const SpotReader *share_from = nullptr);
	~SpotReader();

	// The table's row ids are [FirstRow(), FirstRow() + RowCount()).
	int64_t FirstRow() const;
	uint64_t RowCount() const;
	// Reads one spot and applies fasterq-dump's per-spot checks, throwing SraFormatError where
	// it would stop: READ_START and READ_LEN counts differ, a read runs past the spot, or the
	// spot's bases and qualities differ in length.
	void Read(int64_t row_id, Spot &out);
	// Hints to the byte source that rows [first, end) will be read soon: on S3 the blocks of a
	// work unit are fetched before a thread reaches them.
	void PrefetchRows(int64_t first, int64_t end);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

// fasterq-dump's split modes that the reader offers.
enum class SplitMode {
	// --split-3: mates of a spot with two biological reads go to files 1 and 2; a spot with
	// fewer goes to the unpaired file. Technical reads are always skipped.
	Split3,
	// --split-spot: every read on its own.
	Spot,
};

// A read as fasterq-dump writes it.
struct SraRead {
	int64_t spot = 0;             // $si
	uint32_t read_number = 0;     // $ri: the read's position in its spot, 1-based
	uint32_t destination = 0;     // 0: the unpaired (or only) file; n: file _n
	bool has_name = false;        // the spot's
	std::string name;             // $sn
	std::string bases;            // $rl is its length
	std::vector<uint8_t> quality; // phred, capped at 93 as fasterq-dump prints it ('~')
};

// The reads fasterq-dump writes for one spot, in its order (tbl_join.c). include_technical only
// applies to SplitMode::Spot: split-3 always skips technical reads. Appends to out.
void RouteSpot(const Spot &spot, SplitMode mode, bool include_technical, std::vector<SraRead> &out);

// A row of read_sra's output: one read, or under split-3 a spot's two mates.
struct SraRow {
	int64_t spot = 0;
	uint32_t read_number = 0; // of the first read
	SraRead read1;
	bool paired = false;
	SraRead read2;
};

// Turns one spot's reads (RouteSpot's output, in its order) into rows. Under split-3 a spot is one
// row: mates 1 and 2 together, or its single unpaired read. A spot with a third biological read has
// no place in that shape and raises SraUnsupportedError. Under split-spot each read is a row.
void RowsForSpot(std::vector<SraRead> &reads, SplitMode mode, std::vector<SraRow> &out);

} // namespace sra
} // namespace miint
