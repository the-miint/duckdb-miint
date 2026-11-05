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

void CheckByteOrder(uint32_t tag, const std::string &what) {
	if (tag == ORDER_SWAPPED) {
		throw SraUnsupportedError(what + " is big-endian, which no public run surveyed was");
	}
	if (tag != ORDER_NATIVE) {
		throw SraFormatError(what + " has an invalid byte-order tag");
	}
}

struct Span {
	const uint8_t *data;
	size_t size;
};

// A persisted binary search tree (KDB's PBSTree), the container of both the archive TOC and
// metadata trees: a node count, the payload size, one offset per node (1, 2 or 4 bytes wide by
// payload size), then the payloads. Returns the tree's length in bytes.
size_t ParsePbstree(const uint8_t *buf, size_t len, std::vector<Span> &nodes, const std::string &what) {
	if (len < 4) {
		throw SraFormatError(what + ": truncated tree");
	}
	uint32_t count = LoadU32(buf);
	if (count == 0) {
		return 4;
	}
	if (len < 8) {
		throw SraFormatError(what + ": truncated tree");
	}
	uint32_t data_size = LoadU32(buf + 4);
	size_t width = data_size <= 256 ? 1 : data_size <= 65536 ? 2 : 4;
	size_t index_end = 8 + width * size_t(count);
	if (index_end > len || data_size > len - index_end) {
		throw SraFormatError(what + ": tree runs past its bytes");
	}
	const uint8_t *payload = buf + index_end;
	auto offset_of = [&](size_t i) -> size_t {
		const uint8_t *p = buf + 8 + width * i;
		return width == 1 ? *p : width == 2 ? LoadU16(p) : LoadU32(p);
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

KarArchive::KarArchive(std::shared_ptr<const SraByteSource> source) : source_(std::move(source)) {
	constexpr size_t HEADER = 24;
	const uint64_t file_size = source_->Size();
	if (file_size < HEADER) {
		throw SraFormatError("not an SRA archive: " + std::to_string(file_size) + " bytes is shorter than its header");
	}
	auto header = source_->Read(0, HEADER);
	if (std::memcmp(header.data(), "NCBI.sra", 8) != 0) {
		throw SraFormatError("not an SRA archive: the file does not start with NCBI.sra");
	}
	CheckByteOrder(LoadU32(header.data() + 8), "the archive");
	uint32_t version = LoadU32(header.data() + 12);
	if (version == 0) {
		throw SraFormatError("the archive header records version 0");
	}
	if (version > 1) {
		throw SraUnsupportedError("SRA archive version " + std::to_string(version) + " is not supported");
	}
	data_offset_ = LoadU64(header.data() + 16);
	if (data_offset_ < HEADER + 4 || data_offset_ > file_size) {
		throw SraFormatError("the archive header places its data at offset " + std::to_string(data_offset_) +
		                     ", outside the file");
	}

	auto toc = source_->Read(HEADER, size_t(data_offset_ - HEADER));
	// Each TOC node: name length, name, mtime, mode, type, then a type-specific body. A directory's
	// body is another tree, holding its children.
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
			need(2);
			size_t name_len = LoadU16(p);
			need(2 + name_len + 13);
			std::string name(reinterpret_cast<const char *>(p + 2), name_len);
			uint8_t type = p[2 + name_len + 12];
			p += 2 + name_len + 13;
			left -= 2 + name_len + 13;
			KarEntry entry;
			entry.path = parent.empty() ? name : parent + "/" + name;
			switch (type) {
			case 1: // directory
				entry.type = KarEntry::Type::Directory;
				pending.push_back({Span {p, left}, entry.path});
				break;
			case 2: { // file
				need(16);
				entry.type = KarEntry::Type::File;
				entry.offset = data_offset_ + LoadU64(p);
				entry.size = LoadU64(p + 8);
				if (entry.offset > file_size || entry.size > file_size - entry.offset) {
					throw SraFormatError("archive entry '" + entry.path + "' runs past the end of the file");
				}
				break;
			}
			case 6: // empty file
				entry.type = KarEntry::Type::EmptyFile;
				break;
			case 3:
				throw SraUnsupportedError("archive entry '" + entry.path +
				                          "' is a chunked file, which no public run surveyed used");
			case 4:
			case 5:
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
	case 1:
		return p[0];
	case 2:
		return LoadU16(p);
	case 4:
		return LoadU32(p);
	case 8:
		return LoadU64(p);
	default:
		throw SraFormatError("metadata node '" + name + "' holds " + std::to_string(value.size()) +
		                     " bytes, not an integer");
	}
}

int64_t MetadataNode::ValueAsI64() const {
	uint64_t u = ValueAsU64();
	size_t bits = value.size() * 8;
	if (bits < 64 && (u >> (bits - 1)) & 1) {
		u |= ~uint64_t(0) << bits; // sign-extend
	}
	return int64_t(u);
}

namespace {
// A version 2 metadata node: a byte packing the name length and two flags (has attributes, has
// children), the name, an attribute tree, a child tree, then the value.
MetadataNode ParseMetadataNode(const uint8_t *p, size_t len, int depth) {
	if (depth > 64) {
		throw SraFormatError("metadata nests too deeply");
	}
	if (len < 1) {
		throw SraFormatError("empty metadata node");
	}
	uint8_t bits = p[0];
	size_t name_len = (bits >> 2) + 1;
	if (len < 1 + name_len) {
		throw SraFormatError("truncated metadata node");
	}
	MetadataNode node;
	node.name.assign(reinterpret_cast<const char *>(p + 1), name_len);
	size_t pos = 1 + name_len;
	if (bits & 1) {
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
	if (bits & 2) {
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
	if (bytes.size() < 8) {
		throw SraFormatError("metadata file is shorter than its header");
	}
	CheckByteOrder(LoadU32(bytes.data()), "a metadata file");
	uint32_t version = LoadU32(bytes.data() + 4);
	if (version != 2) {
		throw SraUnsupportedError("metadata version " + std::to_string(version) +
		                          " is not supported (public runs use version 2)");
	}
	std::vector<Span> top;
	ParsePbstree(bytes.data() + 8, bytes.size() - 8, top, "metadata");
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
	if (ar.IsDirectory("tbl")) {
		in_database_ = true;
		if (ar.IsDirectory("db")) {
			throw SraUnsupportedError("the database holds nested databases");
		}
		database_tables_ = ar.Children("tbl");
		if (std::find(database_tables_.begin(), database_tables_.end(), "SEQUENCE") == database_tables_.end()) {
			std::string names;
			for (const auto &t : database_tables_) {
				names += (names.empty() ? "" : ", ") + t;
			}
			throw SraUnsupportedError("the database has no SEQUENCE table (tables: " + names + ")");
		}
		prefix_ = "tbl/SEQUENCE/";
		if (ar.Find("md/cur") != nullptr) {
			db_metadata_ = ParseMetadata(ar.ReadFile("md/cur"));
		}
	} else if (!ar.IsDirectory("col")) {
		throw SraFormatError("the archive holds neither a database nor a table");
	}
	if (ar.Find(prefix_ + "md/cur") == nullptr) {
		throw SraUnsupportedError("the table has no md/cur metadata (a pre-release table)");
	}
	metadata_ = ParseMetadata(ar.ReadFile(prefix_ + "md/cur"));

	std::string col_dir = prefix_ + "col";
	for (const auto &name : ar.Children(col_dir)) {
		if (ar.IsDirectory(col_dir + "/" + name)) {
			physical_columns_.push_back(name);
		}
	}
	if (const MetadataNode *cols = metadata_.Find("col")) {
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
			sc.start_id = start_id ? start_id->ValueAsI64() : 1;
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
