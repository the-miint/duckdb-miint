// =============================================================================
// SylphDatabase.cpp — RAII C++ wrapper around the sylph database FFI.
//
// Translates the C/Rust thread-local-error contract into C++ exceptions:
// load failures throw std::runtime_error with the FFI's error message
// surfaced; the destructor unconditionally frees the handle. See
// SylphDatabase.hpp for the full lifetime / thread-safety contract.
// =============================================================================

#ifdef MIINT_HAS_SYLPH

#include "SylphDatabase.hpp"

#include "sylph.h"

#include <stdexcept>
#include <string>

namespace miint {

SylphDatabaseHandle::SylphDatabaseHandle(const std::string &path) : db_(nullptr) {
	db_ = sylph_database_load(path.c_str());
	if (db_ == nullptr) {
		// `sylph_database_load` populated the FFI's thread-local last-error.
		// Read it before any other sylph_* call to avoid clobber, then
		// surface as a C++ exception. Common messages: "failed to open
		// '<path>'" (missing file), "not a valid .syldb" (corrupt or older
		// sylph version), "contains zero genomes" (empty database).
		const char *err = sylph_get_last_error();
		std::string msg = err ? err : "<no sylph error message>";
		throw std::runtime_error("sylph_database_load failed: " + msg);
	}
}

SylphDatabaseHandle::~SylphDatabaseHandle() {
	// sylph_database_free is documented as NULL-tolerant; the explicit nullptr
	// check is defensive and lets us null out db_ for use-after-free debugging.
	if (db_ != nullptr) {
		sylph_database_free(db_);
		db_ = nullptr;
	}
}

size_t SylphDatabaseHandle::num_genomes() const {
	return sylph_database_num_genomes(db_);
}

}  // namespace miint

#endif  // MIINT_HAS_SYLPH
