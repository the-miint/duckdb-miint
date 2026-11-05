#include "sra_run.hpp"

#include "sra_codecs.hpp"
#include "sra_column.hpp"

#include <algorithm>
#include <cstring>

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

// ---------------------------------------------------------------------------------------------
// Spots

namespace {

// One logical column of the table, from its directory or from a static cell, with the decoded
// blob of the last row read kept for the next.
class ColumnSource {
public:
	ColumnSource(const SraTable &table, const std::string &name, const std::string &expected_type) : name_(name) {
		std::string type;
		if (table.HasPhysical(name)) {
			physical_ = std::make_unique<PhysicalColumn>(table, name);
			type = physical_->Type();
		} else {
			static_ = table.Static(name);
			type = static_->type;
		}
		if (type != expected_type) {
			throw SraUnsupportedError("column '" + name + "' stores '" + type + "', not '" + expected_type + "'");
		}
		type_ = ElementTypeOf(type);
		if (static_) {
			uint64_t bits = static_->size_bits ? *static_->size_bits : uint64_t(static_->row.size()) * 8;
			if (bits % type_.bits != 0 || bits > uint64_t(static_->row.size()) * 8) {
				throw SraFormatError("static column '" + name + "' has a cell of " + std::to_string(bits) + " bits");
			}
			static_count_ = bits / type_.bits;
		}
	}

	static bool Exists(const SraTable &table, const std::string &name) {
		return table.HasPhysical(name) || table.Static(name) != nullptr;
	}

	// The row ids the column covers.
	void Range(int64_t &first, uint64_t &count) const {
		if (static_) {
			first = static_->start_id;
			count = static_->row_count;
			return;
		}
		const auto &blobs = physical_->Index().blobs;
		if (blobs.empty()) {
			first = 1;
			count = 0;
			return;
		}
		first = blobs.front().start_id;
		count = uint64_t(blobs.back().start_id + int64_t(blobs.back().id_count) - first);
	}

	// The row's elements: `count` elements of the column's width, starting `offset` elements
	// into `data` (bit-packed below 8 bits).
	void Cell(int64_t row, const uint8_t *&data, uint64_t &offset, uint64_t &count) {
		if (static_) {
			if (row < static_->start_id || uint64_t(row - static_->start_id) >= static_->row_count) {
				throw SraFormatError("static column '" + name_ + "' has no row " + std::to_string(row));
			}
			data = reinterpret_cast<const uint8_t *>(static_->row.data());
			offset = 0;
			count = static_count_;
			return;
		}
		if (current_ == nullptr || row < current_->start_id ||
		    row >= current_->start_id + int64_t(current_->id_count)) {
			current_ = &physical_->Index().Locate(row);
			blob_ = DecodeBlob(physical_->ReadVdbBlob(*current_), physical_->Encoding(), type_);
		}
		size_t i = size_t(row - current_->start_id);
		data = blob_.data.data();
		offset = blob_.page_map.data_offset[i];
		count = blob_.page_map.row_length[i];
	}

	template <typename T>
	void Values(int64_t row, std::vector<T> &out) {
		const uint8_t *data;
		uint64_t offset, count;
		Cell(row, data, offset, count);
		out.resize(count);
		if (count) {
			std::memcpy(out.data(), data + offset * sizeof(T), count * sizeof(T));
		}
	}

private:
	std::string name_;
	ElementType type_;
	const StaticColumn *static_ = nullptr;
	uint64_t static_count_ = 0;
	std::unique_ptr<PhysicalColumn> physical_;
	const BlobLoc *current_ = nullptr;
	DecodedBlob blob_;
};

// INSDC 4na: each bit a base (A=1, C=2, G=4, T=8), so an ambiguity code is the OR of its bases.
const char FOUR_NA_TEXT[17] = ".ACMGRSVTWYHKDBN";

} // namespace

struct SpotReader::Impl {
	std::shared_ptr<const SraTable> table;
	std::unique_ptr<ColumnSource> read, altread, quality, read_start, read_len, read_type;
	int64_t first = 1;
	uint64_t count = 0;
	std::vector<uint8_t> bases4na;
};

SpotReader::SpotReader(std::shared_ptr<const SraTable> table) : impl_(std::make_unique<Impl>()) {
	Impl &m = *impl_;
	m.table = std::move(table);
	const SraTable &t = *m.table;
	auto need = [&](const char *name, const char *type) {
		if (!ColumnSource::Exists(t, name)) {
			throw SraUnsupportedError(std::string("the table has no ") + name + " column");
		}
		return std::make_unique<ColumnSource>(t, name, type);
	};
	m.read = need("READ", "INSDC:2na:packed");
	if (ColumnSource::Exists(t, "ALTREAD")) {
		m.altread = std::make_unique<ColumnSource>(t, "ALTREAD", "INSDC:4na:bin");
	}
	// The schemas read original qualities first; older tables have only QUALITY.
	const char *qual = ColumnSource::Exists(t, "ORIGINAL_QUALITY") ? "ORIGINAL_QUALITY" : "QUALITY";
	m.quality = need(qual, "INSDC:quality:phred");
	m.read_start = need("READ_START", "INSDC:coord:zero");
	m.read_len = need("READ_LEN", "INSDC:coord:len");
	if (ColumnSource::Exists(t, "READ_TYPE")) {
		m.read_type = std::make_unique<ColumnSource>(t, "READ_TYPE", "INSDC:SRA:xread_type");
	}
	if (const MetadataNode *rna = t.Metadata().Find("RNA_FLAG")) {
		if (rna->value == "1") {
			throw SraUnsupportedError("RNA_FLAG is set (bases would print T as U), which no public run surveyed did");
		}
	}
	// fasterq-dump iterates over the rows of READ.
	m.read->Range(m.first, m.count);
}

SpotReader::~SpotReader() = default;

int64_t SpotReader::FirstRow() const {
	return impl_->first;
}

uint64_t SpotReader::RowCount() const {
	return impl_->count;
}

void SpotReader::Read(int64_t row_id, Spot &out) {
	Impl &m = *impl_;
	out.row_id = row_id;
	const std::string where = "spot " + std::to_string(row_id);

	// READ = 2na bases, two bits each, most significant first, widened to 4na and ORed with the
	// right-aligned ALTREAD (ambiguity codes): bit_or<ALIGN_RIGHT> in the schema.
	const uint8_t *data;
	uint64_t offset, n;
	m.read->Cell(row_id, data, offset, n);
	m.bases4na.resize(n);
	for (uint64_t i = 0; i < n; i++) {
		uint64_t bit = (offset + i) * 2;
		uint8_t two_na = (data[bit >> 3] >> (6 - (bit & 7))) & 3;
		m.bases4na[i] = uint8_t(1u << two_na);
	}
	if (m.altread) {
		const uint8_t *alt;
		uint64_t alt_off, alt_n;
		m.altread->Cell(row_id, alt, alt_off, alt_n);
		if (alt_n > n) {
			throw SraFormatError(where + ": ALTREAD is longer than READ");
		}
		for (uint64_t j = 0; j < alt_n; j++) {
			m.bases4na[n - alt_n + j] |= alt[alt_off + j];
		}
	}
	out.bases.resize(n);
	for (uint64_t i = 0; i < n; i++) {
		out.bases[i] = FOUR_NA_TEXT[m.bases4na[i] & 15];
	}

	m.quality->Values(row_id, out.quality);
	m.read_start->Values(row_id, out.read_start);
	m.read_len->Values(row_id, out.read_len);
	out.has_read_type = m.read_type != nullptr;
	if (m.read_type) {
		m.read_type->Values(row_id, out.read_type);
	} else {
		out.read_type.clear();
	}

	// fasterq-dump's per-spot checks (fq_seq_ua_iter.c); each stops the run.
	if (out.read_start.size() != out.read_len.size()) {
		throw SraFormatError(where + ": READ_START has " + std::to_string(out.read_start.size()) +
		                     " reads but READ_LEN has " + std::to_string(out.read_len.size()));
	}
	for (size_t i = 0; i < out.read_len.size(); i++) {
		if (out.read_start[i] < 0 || uint64_t(out.read_start[i]) + out.read_len[i] > n) {
			throw SraFormatError(where + ": read " + std::to_string(i + 1) + " runs past the spot's " +
			                     std::to_string(n) + " bases");
		}
	}
	if (out.quality.size() != n) {
		throw SraFormatError(where + ": " + std::to_string(n) + " bases but " + std::to_string(out.quality.size()) +
		                     " qualities");
	}
}

namespace {
constexpr uint8_t READ_TYPE_BIOLOGICAL = 1;
constexpr uint8_t MAX_PRINTED_QUALITY = 93; // fasterq-dump prints every q >= 93 as '~'

void Emit(const Spot &spot, uint32_t read_number, uint32_t destination, size_t start, size_t len,
          std::vector<SraRead> &out) {
	SraRead r;
	r.spot = spot.row_id;
	r.read_number = read_number;
	r.destination = destination;
	r.bases = spot.bases.substr(start, len);
	r.quality.resize(len);
	for (size_t i = 0; i < len; i++) {
		r.quality[i] = std::min(spot.quality[start + i], MAX_PRINTED_QUALITY);
	}
	out.push_back(std::move(r));
}

bool Biological(const Spot &spot, size_t i) {
	return (spot.read_type[i] & READ_TYPE_BIOLOGICAL) == READ_TYPE_BIOLOGICAL;
}
} // namespace

void RouteSpot(const Spot &spot, SplitMode mode, bool include_technical, std::vector<SraRead> &out) {
	// split-3 always skips technical reads (perform_fastq_split_3_join).
	bool skip_tech = mode == SplitMode::Split3 || !include_technical;
	const std::string where = "spot " + std::to_string(spot.row_id);
	bool use_types = skip_tech && spot.has_read_type;
	if (use_types && spot.read_type.size() != spot.read_len.size()) {
		throw SraFormatError(where + ": READ_TYPE has " + std::to_string(spot.read_type.size()) +
		                     " reads but READ_LEN has " + std::to_string(spot.read_len.size()));
	}
	const size_t reads = spot.read_len.size();

	// A one-read spot is written whole, tested on its whole length (filter1).
	if (reads == 1) {
		if (spot.bases.empty()) {
			return;
		}
		if (use_types && !spot.read_type.empty() && !Biological(spot, 0)) {
			return;
		}
		Emit(spot, 1, 0, 0, spot.bases.size(), out);
		return;
	}

	if (mode == SplitMode::Spot) {
		// print_fastq_n_reads_split: every non-empty read that passes the filter, in place.
		for (size_t i = 0; i < reads; i++) {
			if (spot.read_len[i] == 0 || (use_types && !Biological(spot, i))) {
				continue;
			}
			Emit(spot, uint32_t(i + 1), 0, size_t(spot.read_start[i]), spot.read_len[i], out);
		}
		return;
	}

	// print_fastq_n_reads_split_3: mates are numbered by the reads written, and a spot with
	// fewer than two non-empty biological reads goes to the unpaired file.
	if (!spot.has_read_type) {
		throw SraUnsupportedError("split-3 needs READ_TYPE, which this table lacks");
	}
	uint32_t valid_reads = 0, valid_bio_reads = 0;
	for (size_t i = 0; i < reads; i++) {
		if (spot.read_len[i] > 0) {
			valid_reads++;
			if (Biological(spot, i)) {
				valid_bio_reads++;
			}
		}
	}
	if (valid_reads == 0) {
		return;
	}
	uint32_t write_id = 1;
	for (size_t i = 0; i < reads; i++) {
		if (spot.read_len[i] == 0 || !Biological(spot, i)) {
			continue;
		}
		if (valid_bio_reads < 2) {
			write_id = 0;
		}
		Emit(spot, uint32_t(i + 1), write_id, size_t(spot.read_start[i]), spot.read_len[i], out);
		if (write_id > 0) {
			write_id++;
		}
	}
}

void RowsForSpot(std::vector<SraRead> &reads, SplitMode mode, std::vector<SraRow> &out) {
	if (reads.empty()) {
		return;
	}
	if (mode == SplitMode::Spot) {
		for (auto &r : reads) {
			SraRow row;
			row.spot = r.spot;
			row.read_number = r.read_number;
			row.read1 = std::move(r);
			out.push_back(std::move(row));
		}
		return;
	}
	// split-3 writes one spot's reads either all to the unpaired file (at most one read: the only
	// non-empty biological one) or to files 1, 2, ...
	SraRow row;
	row.spot = reads[0].spot;
	row.read_number = reads[0].read_number;
	for (auto &r : reads) {
		if (r.destination > 2) {
			throw SraUnsupportedError("spot " + std::to_string(r.spot) +
			                          " has a third biological read, which a split_3 row cannot hold; "
			                          "use split => 'spot'");
		}
	}
	if (reads[0].destination == 0) {
		if (reads.size() != 1) {
			throw SraFormatError("spot " + std::to_string(row.spot) + " routes " + std::to_string(reads.size()) +
			                     " reads to the unpaired file");
		}
		row.read1 = std::move(reads[0]);
	} else {
		if (reads.size() != 2 || reads[0].destination != 1 || reads[1].destination != 2) {
			throw SraFormatError("spot " + std::to_string(row.spot) + " does not route to mates 1 and 2");
		}
		row.read1 = std::move(reads[0]);
		row.read2 = std::move(reads[1]);
		row.paired = true;
	}
	out.push_back(std::move(row));
}

} // namespace sra
} // namespace miint
