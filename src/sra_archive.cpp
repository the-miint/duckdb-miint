#include "sra_archive.hpp"

#include <algorithm>
#include <cstring>

namespace miint {
namespace sra {

namespace {

// Byte-order tags. Every layer of an SRA file (archive, metadata, column index) starts with one;
// read as little-endian, a file written on a big-endian machine shows the swapped value.
constexpr uint32_t ORDER_NATIVE = 0x05031988;
constexpr uint32_t ORDER_SWAPPED = 0x88190305;

// A persisted binary search tree's node offsets are one byte wide while the payload fits in 256
// bytes, two while it fits in 64 KiB, and four otherwise.
constexpr uint32_t PBST_MAX_U8_PAYLOAD = uint32_t(UINT8_MAX) + 1;
constexpr uint32_t PBST_MAX_U16_PAYLOAD = uint32_t(UINT16_MAX) + 1;

// The KAR header: the magic, a byte-order tag, the version and the offset of the file bodies.
constexpr char KAR_MAGIC[] = "NCBI.sra";
constexpr size_t KAR_MAGIC_SIZE = sizeof(KAR_MAGIC) - 1;
constexpr size_t KAR_BYTE_ORDER_AT = 8;
constexpr size_t KAR_VERSION_AT = 12;
constexpr size_t KAR_DATA_OFFSET_AT = 16;
constexpr size_t KAR_HEADER_SIZE = 24;
constexpr uint32_t KAR_VERSION = 1;
// Public archives' tables of contents are a few KB.
constexpr uint64_t KAR_MAX_TOC = uint64_t(64) << 20;

// A TOC node: a 2-byte name length, the name, then an 8-byte mtime, a 4-byte mode and a 1-byte
// entry type, then the type's body. A file's body is its offset and size, 8 bytes each.
constexpr size_t TOC_NAME_LEN_SIZE = 2;
constexpr size_t TOC_TYPE_AT = 12; // after the name
constexpr size_t TOC_ATTRS_SIZE = 13;
constexpr size_t TOC_FILE_BODY_SIZE = 16;
enum TocType : uint8_t {
	TOC_DIRECTORY = 1,
	TOC_FILE = 2,
	TOC_CHUNKED_FILE = 3,
	TOC_SOFT_LINK = 4,
	TOC_HARD_LINK = 5,
	TOC_EMPTY_FILE = 6
};

// A metadata file: a byte-order tag and a version, then a tree of nodes. A node's first byte holds
// its name length less one, shifted left by two, and two flags.
constexpr uint32_t METADATA_VERSION = 2;
constexpr uint8_t METADATA_HAS_ATTRS = 1;
constexpr uint8_t METADATA_HAS_CHILDREN = 2;
constexpr size_t METADATA_NODE_HEADER_SIZE = 1;
constexpr int METADATA_NAME_LEN_SHIFT = 2;
constexpr int MAX_METADATA_DEPTH = 64;

struct Span {
	const uint8_t *data;
	size_t size;
};

// A persisted binary search tree (KDB's PBSTree), the container of both the archive TOC and
// metadata trees: a node count, the payload size, one offset per node (1, 2 or 4 bytes wide by
// payload size), then the payloads. Returns the tree's length in bytes.
size_t ParsePbstree(const uint8_t *buf, size_t len, std::vector<Span> &nodes, const std::string &what) {
	if (len < PBST_COUNT_SIZE) {
		throw SraFormatError(what + ": truncated tree");
	}
	uint32_t count = LoadU32(buf);
	if (count == 0) {
		return PBST_COUNT_SIZE;
	}
	if (len < PBST_HEADER_SIZE) {
		throw SraFormatError(what + ": truncated tree");
	}
	uint32_t data_size = LoadU32(buf + PBST_COUNT_SIZE);
	size_t width = data_size <= PBST_MAX_U8_PAYLOAD    ? sizeof(uint8_t)
	               : data_size <= PBST_MAX_U16_PAYLOAD ? sizeof(uint16_t)
	                                                   : sizeof(uint32_t);
	size_t index_end = PBST_HEADER_SIZE + width * size_t(count);
	if (index_end > len || data_size > len - index_end) {
		throw SraFormatError(what + ": tree runs past its bytes");
	}
	const uint8_t *payload = buf + index_end;
	auto offset_of = [&](size_t i) -> size_t {
		const uint8_t *p = buf + PBST_HEADER_SIZE + width * i;
		return width == sizeof(uint8_t) ? *p : width == sizeof(uint16_t) ? LoadU16(p) : LoadU32(p);
	};
	for (size_t i = 0; i < count; i++) {
		size_t begin = offset_of(i);
		size_t end = i + 1 < count ? offset_of(i + 1) : data_size;
		if (begin > end || end > data_size) {
			throw SraFormatError(what + ": tree node offsets out of order");
		}
		nodes.push_back(Span {payload + begin, end - begin});
	}
	return index_end + data_size;
}

std::string ParentOf(const std::string &path) {
	auto slash = path.rfind('/');
	return slash == std::string::npos ? std::string() : path.substr(0, slash);
}

} // namespace

void CheckByteOrder(uint32_t tag, const std::string &what) {
	if (tag == ORDER_SWAPPED) {
		throw SraUnsupportedError(what + " is big-endian, which no public run surveyed was");
	}
	if (tag != ORDER_NATIVE) {
		throw SraFormatError(what + " has an invalid byte-order tag");
	}
}

KarArchive::KarArchive(std::shared_ptr<const SraByteSource> source) : source_(std::move(source)) {
	const uint64_t file_size = source_->Size();
	if (file_size < KAR_HEADER_SIZE) {
		throw SraFormatError("not an SRA archive: " + std::to_string(file_size) + " bytes is shorter than its header");
	}
	auto header = source_->Read(0, KAR_HEADER_SIZE);
	if (std::memcmp(header.data(), KAR_MAGIC, KAR_MAGIC_SIZE) != 0) {
		throw SraFormatError(std::string("not an SRA archive: the file does not start with ") + KAR_MAGIC);
	}
	CheckByteOrder(LoadU32(header.data() + KAR_BYTE_ORDER_AT), "the archive");
	uint32_t version = LoadU32(header.data() + KAR_VERSION_AT);
	if (version == 0) {
		throw SraFormatError("the archive header records version 0");
	}
	if (version > KAR_VERSION) {
		throw SraUnsupportedError("SRA archive version " + std::to_string(version) + " is not supported");
	}
	data_offset_ = LoadU64(header.data() + KAR_DATA_OFFSET_AT);
	if (data_offset_ < KAR_HEADER_SIZE + PBST_COUNT_SIZE || data_offset_ > file_size) {
		throw SraFormatError("the archive header places its data at offset " + std::to_string(data_offset_) +
		                     ", outside the file");
	}

	// A table of contents past the bound is refused rather than read whole, so a header placing it
	// far into the file cannot make the reader take the file into memory.
	if (data_offset_ - KAR_HEADER_SIZE > KAR_MAX_TOC) {
		throw SraUnsupportedError("a table of contents of " + std::to_string(data_offset_ - KAR_HEADER_SIZE) +
		                          " bytes (public archives' are a few KB)");
	}
	auto toc = source_->Read(KAR_HEADER_SIZE, size_t(data_offset_ - KAR_HEADER_SIZE));
	// A directory's body is another tree, holding its children.
	std::vector<std::pair<Span, std::string>> pending; // (tree bytes, parent path)
	pending.push_back({Span {toc.data(), toc.size()}, std::string()});
	while (!pending.empty()) {
		auto [tree, parent] = pending.back();
		pending.pop_back();
		std::vector<Span> nodes;
		ParsePbstree(tree.data, tree.size, nodes, "the archive's table of contents");
		for (const auto &node : nodes) {
			const uint8_t *p = node.data;
			size_t left = node.size;
			auto need = [&](size_t n) {
				if (left < n) {
					throw SraFormatError("the archive's table of contents is truncated");
				}
			};
			need(TOC_NAME_LEN_SIZE);
			size_t name_len = LoadU16(p);
			const size_t fixed = TOC_NAME_LEN_SIZE + name_len + TOC_ATTRS_SIZE;
			need(fixed);
			std::string name(reinterpret_cast<const char *>(p + TOC_NAME_LEN_SIZE), name_len);
			uint8_t type = p[TOC_NAME_LEN_SIZE + name_len + TOC_TYPE_AT];
			p += fixed;
			left -= fixed;
			KarEntry entry;
			entry.path = parent.empty() ? name : parent + "/" + name;
			switch (type) {
			case TOC_DIRECTORY:
				entry.type = KarEntry::Type::Directory;
				pending.push_back({Span {p, left}, entry.path});
				break;
			case TOC_FILE: {
				need(TOC_FILE_BODY_SIZE);
				entry.type = KarEntry::Type::File;
				uint64_t relative = LoadU64(p);
				entry.size = LoadU64(p + sizeof(uint64_t));
				// Compared before adding, so a huge offset cannot wrap into the file.
				if (relative > file_size - data_offset_ || entry.size > file_size - data_offset_ - relative) {
					throw SraFormatError("archive entry '" + entry.path + "' runs past the end of the file");
				}
				entry.offset = data_offset_ + relative;
				break;
			}
			case TOC_EMPTY_FILE:
				entry.type = KarEntry::Type::EmptyFile;
				break;
			case TOC_CHUNKED_FILE:
				throw SraUnsupportedError("archive entry '" + entry.path +
				                          "' is a chunked file, which no public run surveyed used");
			case TOC_SOFT_LINK:
			case TOC_HARD_LINK:
				throw SraUnsupportedError("archive entry '" + entry.path +
				                          "' is a link, which no public run surveyed used");
			default:
				throw SraFormatError("archive entry '" + entry.path + "' has unknown type " + std::to_string(type));
			}
			entries_.push_back(std::move(entry));
		}
	}
	std::sort(entries_.begin(), entries_.end(), [](const KarEntry &a, const KarEntry &b) { return a.path < b.path; });
	for (size_t i = 0; i < entries_.size(); i++) {
		if (!by_path_.emplace(entries_[i].path, i).second) {
			throw SraFormatError("the archive lists '" + entries_[i].path + "' twice");
		}
	}
}

const KarEntry *KarArchive::Find(const std::string &path) const {
	auto it = by_path_.find(path);
	return it == by_path_.end() ? nullptr : &entries_[it->second];
}

bool KarArchive::IsDirectory(const std::string &path) const {
	const KarEntry *e = Find(path);
	return e != nullptr && e->type == KarEntry::Type::Directory;
}

std::vector<std::string> KarArchive::Children(const std::string &dir) const {
	std::vector<std::string> out;
	if (!IsDirectory(dir)) {
		return out;
	}
	for (const auto &e : entries_) {
		if (ParentOf(e.path) == dir) {
			out.push_back(e.path.substr(dir.size() + 1));
		}
	}
	return out;
}

std::vector<uint8_t> KarArchive::ReadFile(const std::string &path) const {
	const KarEntry *e = Find(path);
	if (e == nullptr || e->type == KarEntry::Type::Directory) {
		throw SraFormatError("the archive has no file '" + path + "'");
	}
	return ReadFile(path, 0, e->size);
}

std::vector<uint8_t> KarArchive::ReadFile(const std::string &path, uint64_t off, uint64_t n) const {
	const KarEntry *e = Find(path);
	if (e == nullptr || e->type == KarEntry::Type::Directory) {
		throw SraFormatError("the archive has no file '" + path + "'");
	}
	if (off > e->size || n > e->size - off) {
		throw SraFormatError("read of " + std::to_string(n) + " bytes at " + std::to_string(off) + " runs past '" +
		                     path + "' (" + std::to_string(e->size) + " bytes)");
	}
	return source_->Read(e->offset + off, size_t(n));
}

// ---------------------------------------------------------------------------------------------
// Metadata

const MetadataNode *MetadataNode::Find(const std::string &path) const {
	const MetadataNode *node = this;
	size_t start = 0;
	while (start <= path.size()) {
		size_t slash = path.find('/', start);
		std::string part = path.substr(start, slash == std::string::npos ? std::string::npos : slash - start);
		const MetadataNode *next = nullptr;
		for (const auto &child : node->children) {
			if (child.name == part) {
				next = &child;
				break;
			}
		}
		if (next == nullptr) {
			return nullptr;
		}
		node = next;
		if (slash == std::string::npos) {
			break;
		}
		start = slash + 1;
	}
	return node;
}

uint64_t MetadataNode::ValueAsU64() const {
	const auto *p = reinterpret_cast<const uint8_t *>(value.data());
	switch (value.size()) {
	case sizeof(uint8_t):
		return p[0];
	case sizeof(uint16_t):
		return LoadU16(p);
	case sizeof(uint32_t):
		return LoadU32(p);
	case sizeof(uint64_t):
		return LoadU64(p);
	default:
		throw SraFormatError("metadata node '" + name + "' holds " + std::to_string(value.size()) +
		                     " bytes, not an integer");
	}
}

int64_t MetadataNode::ValueAsI64() const {
	uint64_t u = ValueAsU64();
	size_t bits = value.size() * BITS_PER_BYTE;
	if (bits < sizeof(u) * BITS_PER_BYTE && (u >> (bits - 1)) & 1) {
		u |= ~uint64_t(0) << bits; // sign-extend
	}
	return int64_t(u);
}

namespace {
// A version 2 metadata node: a byte packing the name length and two flags (has attributes, has
// children), the name, an attribute tree, a child tree, then the value.
MetadataNode ParseMetadataNode(const uint8_t *p, size_t len, int depth) {
	if (depth > MAX_METADATA_DEPTH) {
		throw SraFormatError("metadata nests too deeply");
	}
	if (len < METADATA_NODE_HEADER_SIZE) {
		throw SraFormatError("empty metadata node");
	}
	uint8_t bits = p[0];
	size_t name_len = (bits >> METADATA_NAME_LEN_SHIFT) + 1;
	if (len < METADATA_NODE_HEADER_SIZE + name_len) {
		throw SraFormatError("truncated metadata node");
	}
	MetadataNode node;
	node.name.assign(reinterpret_cast<const char *>(p + METADATA_NODE_HEADER_SIZE), name_len);
	size_t pos = METADATA_NODE_HEADER_SIZE + name_len;
	if (bits & METADATA_HAS_ATTRS) {
		std::vector<Span> attrs;
		pos += ParsePbstree(p + pos, len - pos, attrs, "metadata node '" + node.name + "'");
		for (const auto &a : attrs) {
			const auto *nul = static_cast<const uint8_t *>(std::memchr(a.data, 0, a.size));
			if (nul == nullptr) {
				throw SraFormatError("metadata attribute without a terminated name");
			}
			std::string key(reinterpret_cast<const char *>(a.data), size_t(nul - a.data));
			node.attrs[key].assign(reinterpret_cast<const char *>(nul + 1), a.size - size_t(nul - a.data) - 1);
		}
	}
	if (bits & METADATA_HAS_CHILDREN) {
		std::vector<Span> kids;
		pos += ParsePbstree(p + pos, len - pos, kids, "metadata node '" + node.name + "'");
		for (const auto &k : kids) {
			node.children.push_back(ParseMetadataNode(k.data, k.size, depth + 1));
		}
	}
	node.value.assign(reinterpret_cast<const char *>(p + pos), len - pos);
	return node;
}
} // namespace

MetadataNode ParseMetadata(const std::vector<uint8_t> &bytes) {
	if (bytes.size() < KDB_FILE_HEADER_SIZE) {
		throw SraFormatError("metadata file is shorter than its header");
	}
	CheckByteOrder(LoadU32(bytes.data()), "a metadata file");
	uint32_t version = LoadU32(bytes.data() + sizeof(uint32_t));
	if (version != METADATA_VERSION) {
		throw SraUnsupportedError("metadata version " + std::to_string(version) +
		                          " is not supported (public runs use version 2)");
	}
	std::vector<Span> top;
	ParsePbstree(bytes.data() + KDB_FILE_HEADER_SIZE, bytes.size() - KDB_FILE_HEADER_SIZE, top, "metadata");
	MetadataNode root;
	for (const auto &n : top) {
		root.children.push_back(ParseMetadataNode(n.data, n.size, 0));
	}
	return root;
}

// ---------------------------------------------------------------------------------------------
// Table

SraTable::SraTable(std::shared_ptr<const KarArchive> archive) : archive_(std::move(archive)) {
	const KarArchive &ar = *archive_;
	if (ar.IsDirectory(KDB_TABLES)) {
		in_database_ = true;
		if (ar.IsDirectory(KDB_DATABASES)) {
			throw SraUnsupportedError("the database holds nested databases");
		}
		database_tables_ = ar.Children(KDB_TABLES);
		if (std::find(database_tables_.begin(), database_tables_.end(), SEQUENCE_TABLE) == database_tables_.end()) {
			std::string names;
			for (const auto &t : database_tables_) {
				names += (names.empty() ? "" : ", ") + t;
			}
			throw SraUnsupportedError("the database has no SEQUENCE table (tables: " + names + ")");
		}
		prefix_ = std::string(KDB_TABLES) + "/" + SEQUENCE_TABLE + "/";
		if (ar.Find(KDB_METADATA) != nullptr) {
			db_metadata_ = ParseMetadata(ar.ReadFile(KDB_METADATA));
		}
	} else if (!ar.IsDirectory(KDB_COLUMNS)) {
		throw SraFormatError("the archive holds neither a database nor a table");
	}
	if (ar.Find(prefix_ + KDB_METADATA) == nullptr) {
		throw SraUnsupportedError("the table has no md/cur metadata (a pre-release table)");
	}
	metadata_ = ParseMetadata(ar.ReadFile(prefix_ + KDB_METADATA));

	const std::string col_dir = prefix_ + KDB_COLUMNS;
	for (const auto &name : ar.Children(col_dir)) {
		if (ar.IsDirectory(ColumnDir(name))) {
			physical_columns_.push_back(name);
		}
	}
	if (const MetadataNode *cols = metadata_.Find(KDB_COLUMNS)) {
		for (const auto &node : cols->children) {
			if (HasPhysical(node.name)) {
				continue; // VDB reads a static column only when no column directory exists
			}
			StaticColumn sc;
			auto type = node.attrs.find("type");
			if (type != node.attrs.end()) {
				sc.type = type->second;
			}
			const MetadataNode *row = node.Find("row");
			const MetadataNode *row_count = node.Find("row_count");
			if (row_count == nullptr) {
				throw SraUnsupportedError("static column '" + node.name + "' has no row count (a legacy layout)");
			}
			if (row != nullptr) {
				sc.row = row->value;
			}
			const MetadataNode *start_id = node.Find("start_id");
			sc.start_id = start_id ? start_id->ValueAsI64() : VDB_FIRST_ROW_ID;
			sc.row_count = row_count->ValueAsU64();
			if (const MetadataNode *size = node.Find("size")) {
				sc.size_bits = size->ValueAsU64();
			}
			static_columns_.emplace(node.name, std::move(sc));
		}
	}
}

std::string SraTable::TableType() const {
	const MetadataNode *schema = metadata_.Find("schema");
	if (schema == nullptr) {
		return std::string();
	}
	auto it = schema->attrs.find("name");
	return it == schema->attrs.end() ? std::string() : it->second;
}

bool SraTable::HasPhysical(const std::string &name) const {
	return std::binary_search(physical_columns_.begin(), physical_columns_.end(), name);
}

const StaticColumn *SraTable::Static(const std::string &name) const {
	auto it = static_columns_.find(name);
	return it == static_columns_.end() ? nullptr : &it->second;
}

} // namespace sra
} // namespace miint
