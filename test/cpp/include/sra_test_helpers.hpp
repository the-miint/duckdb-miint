#pragma once

#include "sra_archive.hpp"

#include <fstream>
#include <iterator>
#include <memory>
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
