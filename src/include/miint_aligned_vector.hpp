#pragma once
//
// A std::vector whose storage starts on a fixed, known boundary.
//
// This exists for one reason, and it is a correctness reason rather than a
// performance one. Eigen vectorizes an expression by processing whole packets,
// and for a buffer whose address it does not know it PEELS a few leading scalars
// so the rest can use aligned loads. The number peeled is computed from the
// runtime pointer value, and peeled elements go down a different arithmetic path
// from packeted ones -- scalar `std::exp` against a packet polynomial, a
// different reduction order in a matrix product. So the same computation on the
// same data returns a different answer depending on WHERE THE ALLOCATOR PUT THE
// BUFFER.
//
// At SSE2 that never bites: packets are 16 bytes, malloc already guarantees 16,
// so nothing is ever peeled. At AVX2 (32) and AVX-512 (64) malloc's guarantee is
// not enough, and mmvec's run-to-run reproducibility -- the guarantee its own RNG
// transforms exist to provide -- quietly disappears. Measured: fitting the same
// data with the same seed twice in one process gave two different models, because
// the second fit reused the first one's freed blocks at a different offset.
//
// Aligning to the widest packet miint compiles for fixes it at the source: the
// peel count becomes a function of the shape alone, identical on every run. The
// alternative -- EIGEN_MAX_ALIGN_BYTES=0, which makes Eigen use unaligned access
// everywhere and peel nothing -- was measured at +75% (AVX2) and +113% (AVX-512),
// which puts AVX-512 behind the SSE2 baseline. Alignment is the cheap fix.

#include <cstddef>
#include <new>
#include <vector>

namespace miint {

//! The widest packet any miint variant is compiled for (AVX-512 = 64 bytes), and
//! therefore the alignment that makes every variant's peel count deterministic.
//! Also a cache line on every CPU miint targets.
inline constexpr std::size_t kSimdAlignment = 64;

//! Minimal allocator that hands back over-aligned storage.
//!
//! Uses the C++17 aligned operator new rather than posix_memalign or
//! std::aligned_alloc so it is portable to Emscripten and MSVC, both of which are
//! build targets here.
template <typename T, std::size_t Alignment = kSimdAlignment>
struct AlignedAllocator {
	using value_type = T;

	AlignedAllocator() noexcept = default;
	template <typename U>
	explicit constexpr AlignedAllocator(const AlignedAllocator<U, Alignment> &) noexcept {
	}

	template <typename U>
	struct rebind {
		using other = AlignedAllocator<U, Alignment>;
	};

	T *allocate(std::size_t n) {
		return static_cast<T *>(::operator new (n * sizeof(T), std::align_val_t {Alignment}));
	}

	void deallocate(T *p, std::size_t n) noexcept {
		::operator delete (p, n * sizeof(T), std::align_val_t {Alignment});
	}

	// Stateless, so any two instances can free each other's memory.
	bool operator==(const AlignedAllocator &) const noexcept {
		return true;
	}
	bool operator!=(const AlignedAllocator &) const noexcept {
		return false;
	}
};

//! Drop-in for std::vector<T> with SIMD-aligned storage. The full vector
//! interface is unchanged -- only the allocator differs -- so `.data()`,
//! `.resize()` and `operator[]` behave exactly as before.
template <typename T>
using AlignedVector = std::vector<T, AlignedAllocator<T>>;

} // namespace miint
