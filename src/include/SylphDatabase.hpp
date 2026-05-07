#pragma once
#ifdef MIINT_HAS_SYLPH

#include <stddef.h>

#include <memory>
#include <string>

// Forward-declare the opaque sylph FFI type so callers don't need sylph.h.
// The full definition lives in the Rust crate and is not part of any C++
// translation unit's interface.
struct SylphDatabase;

namespace miint {

// RAII wrapper around the sylph_database_load / sylph_database_free FFI pair.
//
// Constructs by loading a `.syldb` (sylph's bincode-serialized
// `Vec<GenomeSketch>`) from disk; throws `std::runtime_error` with the FFI's
// thread-local error message on failure (typically "failed to open"
// for missing files, "not a valid .syldb" for corrupt or version-mismatched
// archives). The handle is freed automatically by the destructor.
//
// Thread safety: loading and freeing are NOT thread-safe (mirrors the FFI
// contract). Once loaded, multiple threads may share a const reference and
// pass `raw()` to `sylph_profile` concurrently — sylph treats the database
// as read-only during profiling.
//
// The wrapper is intentionally non-copyable and non-movable so the
// destructor's "free exactly once" invariant is trivial to verify by
// inspection.
class SylphDatabaseHandle {
public:
	// Throws std::runtime_error("sylph_database_load failed: <ffi error>")
	// on load failure.
	explicit SylphDatabaseHandle(const std::string &path);
	~SylphDatabaseHandle();

	SylphDatabaseHandle(const SylphDatabaseHandle &) = delete;
	SylphDatabaseHandle &operator=(const SylphDatabaseHandle &) = delete;
	SylphDatabaseHandle(SylphDatabaseHandle &&) = delete;
	SylphDatabaseHandle &operator=(SylphDatabaseHandle &&) = delete;

	// Number of reference genome sketches in the database. Constant after
	// construction.
	size_t num_genomes() const;

	// Borrow the raw FFI pointer (non-owning). Lifetime is tied to this
	// wrapper — do NOT call sylph_database_free on the returned pointer;
	// the destructor does that.
	const ::SylphDatabase *raw() const noexcept {
		return db_;
	}

private:
	::SylphDatabase *db_;
};

}  // namespace miint

#endif  // MIINT_HAS_SYLPH
