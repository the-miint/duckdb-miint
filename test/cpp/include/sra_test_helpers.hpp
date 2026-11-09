#pragma once

#include "sra_archive.hpp"
#include "sra_column.hpp"

#include <fstream>
#include <iterator>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace miint_test {

// The bytes of a fixture. C++ tests run from the repository root, so data/ paths resolve.
inline std::vector<uint8_t> SlurpBytes(const std::string &path) {
	std::ifstream in(path, std::ios::binary);
	if (in.fail()) {
		throw std::runtime_error("cannot open " + path + " (run the tests from the repository root)");
	}
	return std::vector<uint8_t>(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// A fixture from data/sra (data/sra/README.md describes each one).
inline std::shared_ptr<const miint::sra::KarArchive> OpenSraFixture(const std::string &file) {
	return std::make_shared<const miint::sra::KarArchive>(
	    std::make_shared<miint::sra::SraFileSource>("data/sra/" + file));
}

inline std::shared_ptr<const miint::sra::KarArchive> SraArchiveFromBytes(std::vector<uint8_t> bytes) {
	return std::make_shared<const miint::sra::KarArchive>(
	    std::make_shared<miint::sra::SraMemorySource>(std::move(bytes)));
}

// A page map's rows expanded, for comparing small blobs with the values a test expects.
inline std::vector<uint32_t> RowLengths(const miint::sra::PageMap &pm) {
	std::vector<uint32_t> out;
	for (uint64_t i = 0; i < pm.row_count; i++) {
		uint32_t length;
		uint64_t offset;
		pm.Row(i, length, offset);
		out.push_back(length);
	}
	return out;
}
inline std::vector<uint64_t> RowOffsets(const miint::sra::PageMap &pm) {
	std::vector<uint64_t> out;
	for (uint64_t i = 0; i < pm.row_count; i++) {
		uint32_t length;
		uint64_t offset;
		pm.Row(i, length, offset);
		out.push_back(offset);
	}
	return out;
}
// A page map whose rows each have their own record, of the given lengths, end to end.
inline miint::sra::PageMap PackedRows(const std::vector<uint32_t> &lengths) {
	miint::sra::PageMap pm;
	uint64_t offset = 0;
	for (uint32_t l : lengths) {
		pm.runs.push_back({pm.row_count++, 1, l, offset, l});
		offset += l;
	}
	return pm;
}

// True if fn throws exactly an E whose message contains fragment. The SRA reader keeps damage
// (SraFormatError) apart from deliberate refusal (SraUnsupportedError), and its messages name
// what was found, so both the type and the fact are asserted.
template <typename E, typename Fn>
bool ThrowsTypeMentioning(Fn &&fn, const std::string &fragment) {
	try {
		fn();
	} catch (const E &e) {
		return std::string(e.what()).find(fragment) != std::string::npos;
	} catch (...) {
		return false;
	}
	return false;
}

} // namespace miint_test
