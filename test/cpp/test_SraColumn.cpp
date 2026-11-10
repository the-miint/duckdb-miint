#include <catch2/catch_test_macros.hpp>

#include "sra_column.hpp"
#include "sra_test_helpers.hpp"

#include <zlib.h>

#include <algorithm>
#include <cstring>
#include <filesystem>
#include <string>
#include <vector>

using namespace miint::sra;
using miint_test::OpenSraFixture;
using miint_test::RowLengths;
using miint_test::RowOffsets;
using miint_test::SlurpBytes;
using miint_test::SraArchiveFromBytes;
using miint_test::ThrowsTypeMentioning;

// Expected values were read from the fixtures with an independent parser written while the
// format was worked out (data/sra/README.md describes the fixtures).

namespace {
void PutU32(std::vector<uint8_t> &out, uint32_t v) {
	for (int i = 0; i < 4; i++) {
		out.push_back(uint8_t(v >> (8 * i)));
	}
}
void PutU64(std::vector<uint8_t> &out, uint64_t v) {
	for (int i = 0; i < 8; i++) {
		out.push_back(uint8_t(v >> (8 * i)));
	}
}
void PutKdbHeader(std::vector<uint8_t> &out, uint32_t version) {
	PutU32(out, 0x05031988);
	PutU32(out, version);
}

// A version 3 column header (the `idx` file).
std::vector<uint8_t> ColumnHeader(uint64_t data_eof, uint64_t idx2_eof, uint32_t idx0_count, uint32_t num_blocks,
                                  uint8_t checksum, uint32_t version = 3) {
	std::vector<uint8_t> idx;
	PutKdbHeader(idx, version);
	PutU64(idx, data_eof);
	PutU64(idx, idx2_eof);
	PutU32(idx, idx0_count);
	PutU32(idx, num_blocks);
	PutU32(idx, 1); // page size
	idx.push_back(checksum);
	idx.insert(idx.end(), 3, 0);
	return idx;
}

// A 24-byte locator, as idx0 (blob) and idx1 (block) store them.
void PutLoc(std::vector<uint8_t> &out, uint64_t pg, uint32_t gen, uint32_t id_range, int64_t start_id) {
	PutU64(out, pg);
	PutU32(out, gen);
	PutU32(out, id_range);
	PutU64(out, uint64_t(start_id));
}

std::vector<uint8_t> Idx1Empty() {
	std::vector<uint8_t> idx1;
	PutKdbHeader(idx1, 3);
	return idx1;
}

void PutVlenU(std::vector<uint8_t> &out, uint64_t v) {
	uint8_t tmp[10];
	int n = 0;
	do {
		tmp[n++] = uint8_t(v & 0x7f);
		v >>= 7;
	} while (v != 0);
	for (int i = n - 1; i >= 0; i--) {
		out.push_back(uint8_t(tmp[i] | (i > 0 ? 0x80 : 0)));
	}
}

std::vector<uint8_t> RawDeflate(const std::vector<uint8_t> &in) {
	z_stream zs {};
	REQUIRE(deflateInit2(&zs, Z_BEST_SPEED, Z_DEFLATED, -15, 9, Z_DEFAULT_STRATEGY) == Z_OK);
	std::vector<uint8_t> out(in.size() + 64);
	zs.next_in = const_cast<Bytef *>(in.data());
	zs.avail_in = uInt(in.size());
	zs.next_out = out.data();
	zs.avail_out = uInt(out.size());
	REQUIRE(deflate(&zs, Z_FINISH) == Z_STREAM_END);
	out.resize(zs.total_out);
	deflateEnd(&zs);
	return out;
}

// A v2 blob with no encoding headers: header byte, hdr_size, map_size, page map, payload.
std::vector<uint8_t> V2Blob(const std::vector<uint8_t> &page_map, const std::vector<uint8_t> &payload) {
	std::vector<uint8_t> b;
	b.push_back(0x80); // v2, little-endian, variant 0, whole final byte
	b.push_back(0);
	b.push_back(uint8_t(page_map.size()));
	b.insert(b.end(), page_map.begin(), page_map.end());
	b.insert(b.end(), payload.begin(), payload.end());
	return b;
}

// Every committed .sra (data/sra/README.md), accepted or not: the container layer must read them
// all, since classification itself reads the container.
std::vector<std::string> AllFixtures() {
	std::vector<std::string> files;
	for (const auto &e : std::filesystem::directory_iterator("data/sra")) {
		auto ext = e.path().extension().string();
		if (ext == ".sra" || ext == ".sralite") {
			files.push_back(e.path().filename().string());
		}
	}
	std::sort(files.begin(), files.end());
	return files;
}
} // namespace

TEST_CASE("Crc32Ncbi is VDB's CRC, not zlib's", "[sra]") {
	// Check value from the reference implementation; zlib's crc32 of the same input is
	// 0xcbf43926, so a reader using it would reject every checksummed blob.
	const char *s = "123456789";
	CHECK(Crc32Ncbi(reinterpret_cast<const uint8_t *>(s), 9) == 0x89a1897f);
	CHECK(Crc32Ncbi(nullptr, 0) == 0);
}

TEST_CASE("PhysicalColumn reads a column's type, encoding and index", "[sra]") {
	SraTable t(OpenSraFixture("SRR4424282.sra"));
	PhysicalColumn read(t, "READ");
	CHECK(read.Type() == "INSDC:2na:packed");
	CHECK(read.Encoding().empty());
	CHECK(read.Index().data_eof == 382);
	CHECK(read.Index().page_size == 1);
	CHECK(read.Index().checksum == ChecksumType::Crc32);
	REQUIRE(read.Index().blobs.size() == 1);
	const BlobLoc &b = read.Index().blobs[0];
	CHECK(b.start_id == 1);
	CHECK(b.id_count == 3);
	CHECK(b.pg == 0);
	CHECK(b.size == 382);

	PhysicalColumn qual(t, "QUALITY");
	CHECK(qual.Type() == "INSDC:quality:phred");
	CHECK(qual.Encoding() == "<INSDC:quality:phred>zip_encoding#1");
}

TEST_CASE("A v1 blob gives every row the fixed length its header records", "[sra]") {
	SraTable t(OpenSraFixture("SRR4424282.sra"));
	PhysicalColumn read(t, "READ");
	REQUIRE_FALSE(read.Index().blobs.empty());
	VdbBlob blob = read.ReadVdbBlob(read.Index().blobs[0]);
	CHECK(blob.frames.empty());
	CHECK(blob.big_endian);
	CHECK(RowLengths(blob.page_map) == std::vector<uint32_t> {500, 500, 500});
	CHECK(RowOffsets(blob.page_map) == std::vector<uint64_t> {0, 500, 1000});
	// 382 stored bytes = 3 header bytes + 375 bytes of 2-bit bases + a 4-byte CRC.
	CHECK(blob.PayloadSize() == 375);
	CHECK(blob.payload_bits == 3000);
}

TEST_CASE("A v2 blob's page map shares one record among repeated rows", "[sra]") {
	// ERR932865's READ: 148 rows of 40 bases, but identical consecutive reads are stored once,
	// so only 132 records follow. Rows must map onto records, not onto 148 slots.
	SraTable t(OpenSraFixture("ERR932865.sra"));
	PhysicalColumn read(t, "READ");
	REQUIRE_FALSE(read.Index().blobs.empty());
	VdbBlob blob = read.ReadVdbBlob(read.Index().blobs[0]);
	REQUIRE(blob.page_map.row_count == 148);
	auto offsets = RowOffsets(blob.page_map);
	CHECK(offsets[0] == 0);
	CHECK(offsets[1] == 40);
	CHECK(offsets[2] == 80);
	CHECK(offsets[3] == 80); // data_run [1, 1, 2, ...]: rows 2 and 3 share
	CHECK(blob.page_map.ElementsSpanned() == 132 * 40);
	CHECK(blob.payload_bits == 132 * 40 * 2);
}

TEST_CASE("A v2 blob with variable row lengths and padding bits", "[sra]") {
	SraTable t(OpenSraFixture("SRR32151205.sra"));
	PhysicalColumn read(t, "READ");
	REQUIRE(read.Index().blobs.size() == 1);
	VdbBlob blob = read.ReadVdbBlob(read.Index().blobs[0]);
	CHECK(blob.page_map.row_count == 1690);
	// 227,706 payload bytes, the last 6 bits padding: 910,821 bases.
	CHECK(blob.payload_bits == 227706u * 8 - 6);
	CHECK(blob.page_map.ElementsSpanned() * 2 == blob.payload_bits);
}

TEST_CASE("PrefetchRows hints the bytes of the blobs holding the rows", "[sra]") {
	// On S3 the scan hints a work unit's blobs when the unit is claimed, so they are fetched
	// before the thread reaches them. The hint must name the blobs' bytes in the file.
	struct HintRecorder : SraFileSource {
		using SraFileSource::SraFileSource;
		void Prefetch(uint64_t pos, uint64_t n) const override {
			hints.emplace_back(pos, n);
		}
		mutable std::vector<std::pair<uint64_t, uint64_t>> hints;
	};
	auto source = std::make_shared<HintRecorder>("data/sra/SRR4424282.sra");
	SraTable t(std::make_shared<const KarArchive>(source));
	PhysicalColumn read(t, "READ");
	read.PrefetchRows(1, 4);
	// READ/data starts at 2612; its one blob holds all three rows in 382 bytes.
	REQUIRE(source->hints.size() == 1);
	CHECK(source->hints[0] == std::make_pair<uint64_t, uint64_t>(2612, 382));
	source->hints.clear();
	read.PrefetchRows(4, 9); // past the last row: nothing
	CHECK(source->hints.empty());
}

TEST_CASE("ReadBlob rejects a blob whose CRC32 does not match", "[sra]") {
	auto bytes = SlurpBytes("data/sra/SRR4424282.sra");
	// READ/data starts at 2612; flip one base inside its only blob.
	bytes[2612 + 10] ^= 0x01;
	SraTable t(SraArchiveFromBytes(bytes));
	PhysicalColumn read(t, "READ");
	REQUIRE_FALSE(read.Index().blobs.empty());
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { read.ReadBlob(read.Index().blobs[0]); }, "checksum"));
}

TEST_CASE("Every blob of every fixture is located, verified and framed", "[sra]") {
	// The container layer's done-criterion: each blob's page map covers its rows, and the
	// payload holds exactly the elements the page map addresses. For an unencoded payload that
	// is its bit length; for an encoded one it is the size the first decoding stage declares.
	size_t columns = 0;
	size_t blobs = 0;
	std::vector<std::string> skipped;
	for (const auto &file : AllFixtures()) {
		INFO(file);
		SraTable t(OpenSraFixture(file));
		for (const auto &name : t.PhysicalColumns()) {
			INFO(name);
			PhysicalColumn col(t, name);
			uint32_t elem_bits;
			try {
				elem_bits = ElementBits(col.Type());
			} catch (const SraUnsupportedError &) {
				skipped.push_back(file + ":" + name); // a type no supported run stores
				continue;
			}
			const auto &locs = col.Index().blobs;
			REQUIRE_FALSE(locs.empty());
			for (size_t i = 1; i < locs.size(); i++) {
				CHECK(locs[i - 1].start_id + int64_t(locs[i - 1].id_count) == locs[i].start_id);
			}
			bool single_stage = col.Encoding().empty() || col.Encoding().find("zip_encoding#1") != std::string::npos;
			for (const auto &loc : locs) {
				VdbBlob blob = col.ReadVdbBlob(loc);
				REQUIRE(blob.page_map.row_count == loc.id_count);
				uint64_t elems = blob.page_map.ElementsSpanned();
				if (blob.frames.empty() || (blob.frames[0].flags & 1)) {
					CHECK(blob.payload_bits == elems * elem_bits);
				} else if (single_stage) {
					CHECK(uint64_t(blob.frames[0].osize) == (elems * elem_bits + 7) / 8);
				}
				blobs++;
			}
			columns++;
		}
	}
	// The fixtures are fixed, so these counts are too. Only the rejection fixtures store types the
	// reader has no use for (454 flow signal, alignment ids).
	CHECK(columns == 91);
	CHECK(blobs > columns);
	for (const auto &s : skipped) {
		INFO(s);
		CHECK((s.rfind("SRR882123.sra:", 0) == 0 || s.rfind("SRR26779001.sra:", 0) == 0 ||
		       s.rfind("SRR2513321.sra:", 0) == 0));
	}
}

TEST_CASE("ColumnIndex reads blobs recorded only in idx0", "[sra]") {
	// About 2.4% of surveyed runs keep blobs in the idx0 journal with nothing in idx1/idx2; no
	// committed fixture does (data/sra/README.md), so the layout is built here.
	std::vector<uint8_t> idx0;
	PutLoc(idx0, 0, 100, 10, 1);   // rows 1..10 in 100 bytes at 0
	PutLoc(idx0, 100, 80, 5, 11);  // rows 11..15 in 80 bytes at 100
	PutLoc(idx0, 180, 120, 10, 1); // rows 1..10 rewritten: a later entry replaces an earlier one
	ColumnIndex ix = ColumnIndex::Parse(ColumnHeader(300, 0, 3, 0, 0), idx0, Idx1Empty(), {});
	REQUIRE(ix.blobs.size() == 2);
	CHECK(ix.blobs[0].start_id == 1);
	CHECK(ix.blobs[0].pg == 180);
	CHECK(ix.blobs[0].size == 120);
	CHECK(ix.blobs[1].start_id == 11);
	CHECK(ix.Locate(13).start_id == 11);
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { ix.Locate(16); }, "16"));
}

TEST_CASE("ColumnIndex rejects idx0 entries that overlap without matching", "[sra]") {
	std::vector<uint8_t> idx0;
	PutLoc(idx0, 0, 100, 10, 1);
	PutLoc(idx0, 100, 100, 10, 5);
	CHECK_THROWS_AS(ColumnIndex::Parse(ColumnHeader(200, 0, 2, 0, 0), idx0, Idx1Empty(), {}), SraFormatError);
}

TEST_CASE("ColumnIndex refuses index variants absent from public runs", "[sra]") {
	SECTION("column header version 2") {
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
		    [&] { ColumnIndex::Parse(ColumnHeader(0, 0, 0, 0, 0, 2), {}, Idx1Empty(), {}); }, "version 2"));
	}
	SECTION("a Random-type idx2 block") {
		std::vector<uint8_t> idx1 = Idx1Empty();
		PutLoc(idx1, 0, 24, 10, 1); // id and pg types 0 (Random)
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
		    [&] { ColumnIndex::Parse(ColumnHeader(0, 24, 0, 1, 0), {}, idx1, std::vector<uint8_t>(24)); },
		    "idx2 block"));
	}
	SECTION("a compressed idx2 block") {
		std::vector<uint8_t> idx1 = Idx1Empty();
		PutLoc(idx1, 0, 12u | (3u << 27) | (3u << 29) | (1u << 31), 10, 1);
		CHECK(ThrowsTypeMentioning<SraUnsupportedError>(
		    [&] { ColumnIndex::Parse(ColumnHeader(0, 12, 0, 1, 0), {}, idx1, std::vector<uint8_t>(12)); },
		    "compressed"));
	}
}

TEST_CASE("A random-access page map gives each row its own offset", "[sra]") {
	// Version 2 page maps record one data offset per row, so rows can share or reorder stored
	// data freely. Variant 0: fixed row length 2; rows 0 and 2 share the record at 0.
	std::vector<uint8_t> body;
	for (uint64_t off : {0, 2, 0}) {
		PutVlenU(body, off);
	}
	std::vector<uint8_t> pm {uint8_t((2 << 2) | 0)};
	PutVlenU(pm, 2);
	auto deflated = RawDeflate(body);
	pm.insert(pm.end(), deflated.begin(), deflated.end());
	VdbBlob blob = ParseVdbBlob(V2Blob(pm, std::vector<uint8_t>(16)), 3, 32);
	CHECK(RowLengths(blob.page_map) == std::vector<uint32_t> {2, 2, 2});
	CHECK(RowOffsets(blob.page_map) == std::vector<uint64_t> {0, 2, 0});
	CHECK(blob.page_map.ElementsSpanned() == 4);
}

TEST_CASE("A version 0 page map is read without inflating", "[sra]") {
	// Variant 2 (variable lengths), stored raw: lengths [3, 5] repeated [1, 2] times.
	std::vector<uint8_t> pm {uint8_t((0 << 2) | 2)};
	PutVlenU(pm, 2);
	for (uint64_t v : {3, 5, 1, 2}) {
		PutVlenU(pm, v);
	}
	VdbBlob blob = ParseVdbBlob(V2Blob(pm, std::vector<uint8_t>(13)), 3, 8);
	CHECK(RowLengths(blob.page_map) == std::vector<uint32_t> {3, 5, 5});
	CHECK(RowOffsets(blob.page_map) == std::vector<uint64_t> {0, 3, 8});
	CHECK(blob.payload_bits == 13 * 8);
}

TEST_CASE("A page map whose runs do not cover the blob's rows is damage", "[sra]") {
	std::vector<uint8_t> pm {uint8_t((0 << 2) | 2)};
	PutVlenU(pm, 1);
	PutVlenU(pm, 3);
	PutVlenU(pm, 2); // two rows of 3, but the blob holds three rows
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { ParseVdbBlob(V2Blob(pm, std::vector<uint8_t>(9)), 3, 8); },
	                                           "page map"));
}

TEST_CASE("ElementBits knows the stored types a sequence table uses", "[sra]") {
	CHECK(ElementBits("INSDC:2na:packed") == 2);
	CHECK(ElementBits("INSDC:4na:bin") == 8);
	CHECK(ElementBits("INSDC:quality:phred") == 8);
	CHECK(ElementBits("INSDC:coord:len") == 32);
	CHECK(ElementBits("INSDC:coord:zero") == 32);
	CHECK(ElementBits("INSDC:coord:val") == 32);
	CHECK(ElementBits("INSDC:SRA:xread_type") == 8);
	CHECK(ElementBits("ascii") == 8);
	CHECK_THROWS_AS(ElementBits("NCBI:isamp1"), SraUnsupportedError);
}

// Counts in a file are only trusted as far as the bytes behind them: a damaged or crafted file
// must fail on its own bytes, quickly, and never make the reader allocate in proportion to a
// count it merely claims.

TEST_CASE("A page map of empty rows needs no memory per row", "[sra]") {
	// VDB cuts blobs by size, so a column of mostly empty cells (ALTREAD) can put millions of rows
	// in one blob, described by a handful of bytes. 2^40 rows must not mean 2^40 entries.
	std::vector<uint8_t> pm {uint8_t((0 << 2) | 0)}; // version 0, variant 0: every row the same
	PutVlenU(pm, 0);                                 // length 0
	const uint64_t rows = uint64_t(1) << 40;
	VdbBlob blob = ParseVdbBlob(V2Blob(pm, {}), rows, 8);
	CHECK(blob.page_map.row_count == rows);
	uint32_t length = 7;
	uint64_t offset = 7;
	blob.page_map.Row(rows - 1, length, offset);
	CHECK(length == 0);
}

TEST_CASE("A page map's counts are bounded by its bytes", "[sra]") {
	// 2^31 length runs claimed, with a stream that inflates to a few bytes: each stored value
	// takes at least a byte, so the claim is refused before any buffer is sized from it.
	std::vector<uint8_t> pm {uint8_t((1 << 2) | 2)}; // version 1, variant 2
	PutVlenU(pm, uint64_t(1) << 31);
	auto deflated = RawDeflate({1, 1});
	pm.insert(pm.end(), deflated.begin(), deflated.end());
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { ParseVdbBlob(V2Blob(pm, {}), 1, 8); }, "claims"));
}

TEST_CASE("A page map's counts cannot wrap past the bound", "[sra]") {
	// 2^63 length runs: doubled for lengths and run counts, the total wraps to 0 and would pass a
	// bound applied to the sum. Each count is bounded before they are added.
	std::vector<uint8_t> pm {uint8_t((0 << 2) | 2)}; // version 0, variant 2
	PutVlenU(pm, uint64_t(1) << 63);
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { ParseVdbBlob(V2Blob(pm, {}), 1, 8); }, "too many values"));
}

TEST_CASE("A page map larger than the first inflate buffer is read whole", "[sra]") {
	// Random-access offsets that repeat deflate to back-references; 5,000 of them inflate past the
	// first 4 KiB buffer, so the buffer grows mid-stream and must keep what it already holds.
	const uint64_t rows = 5000;
	std::vector<uint8_t> body;
	std::vector<uint64_t> expected;
	for (uint64_t r = 0; r < rows; r++) {
		uint64_t off = (r % 7) * 300;
		expected.push_back(off);
		PutVlenU(body, off);
	}
	REQUIRE(body.size() > 4096);
	std::vector<uint8_t> pm {uint8_t((2 << 2) | 0)}; // version 2 (random access), variant 0
	PutVlenU(pm, 2);                                 // row length 2
	auto deflated = RawDeflate(body);
	pm.insert(pm.end(), deflated.begin(), deflated.end());
	std::vector<uint8_t> blob {0x80, 0, 0, 0, 0}; // v2, variant 2 header: 4-byte map size
	blob[0] = 0x80 | (2 << 4);
	const uint32_t map_size = uint32_t(pm.size());
	blob[2] = uint8_t(map_size);
	blob[3] = uint8_t(map_size >> 8);
	blob[4] = uint8_t(map_size >> 16);
	blob.push_back(uint8_t(map_size >> 24));
	blob.insert(blob.end(), pm.begin(), pm.end());
	blob.resize(blob.size() + 6 * 300 + 2, 0);
	VdbBlob parsed = ParseVdbBlob(blob, rows, 8);
	CHECK(RowOffsets(parsed.page_map) == expected);
}

TEST_CASE("An index block cannot hold more blobs than rows", "[sra]") {
	// A predictable block records its blob count in the locator; every blob holds a row at least.
	std::vector<uint8_t> idx1 = Idx1Empty();
	PutLoc(idx1, 0, (uint32_t(1) << 20) | (3u << 27) | (3u << 29), 10, 1); // 2^20 blobs, 10 rows
	std::vector<uint8_t> idx2(12);
	CHECK(ThrowsTypeMentioning<SraFormatError>(
	    [&] { ColumnIndex::Parse(ColumnHeader(1 << 30, 12, 0, 1, 0), {}, idx1, idx2); }, "more blobs"));
}

TEST_CASE("A page map that inflates past what its values can need is damage", "[sra]") {
	// One random-access offset needs at most ten bytes, yet the stream holds a hundred: the reader
	// stops at what the map's own counts allow instead of following the stream.
	std::vector<uint8_t> pm {uint8_t((2 << 2) | 0)}; // version 2 (random access), variant 0
	PutVlenU(pm, 1);                                 // row length 1
	auto deflated = RawDeflate(std::vector<uint8_t>(100, 0));
	pm.insert(pm.end(), deflated.begin(), deflated.end());
	CHECK(ThrowsTypeMentioning<SraFormatError>([&] { ParseVdbBlob(V2Blob(pm, std::vector<uint8_t>(4)), 1, 8); },
	                                           "inflates past"));
}
