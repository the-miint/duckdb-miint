#pragma once

#include <cstddef>
#include <string>

namespace duckdb {
namespace miint {
namespace gpl_boundary {

/// Input shm segment that miint creates and the daemon reads.
///
/// Lifecycle: miint allocates → miint writes Arrow IPC bytes via `data()` →
/// miint sends `{shm_input, shm_input_size}` to the daemon → daemon mmaps
/// and reads → miint unlinks (via the destructor of this object).
///
/// The daemon does NOT call `fstat` on the inherited fd to derive size
/// (gpl-boundary commit `19306f6`); it trusts the explicit `shm_input_size`
/// the caller sent in the JSON request. This class exposes `size_bytes()`
/// so that field can be set authoritatively.
class InputShmRegion {
public:
	/// Allocate a new shm segment of `size_bytes` with a unique name. Throws
	/// `std::runtime_error` on shm_open / ftruncate / mmap failure.
	static InputShmRegion Create(std::size_t size_bytes);

	~InputShmRegion();

	InputShmRegion(const InputShmRegion &) = delete;
	InputShmRegion &operator=(const InputShmRegion &) = delete;
	InputShmRegion(InputShmRegion &&) noexcept;
	InputShmRegion &operator=(InputShmRegion &&) noexcept;

	/// POSIX shm name, e.g. `/miint-input-<uuid>`. Pass to the daemon's
	/// `BatchRequest.shm_input` JSON field.
	const std::string &name() const {
		return name_;
	}

	/// Number of bytes mapped. Pass to the daemon's `BatchRequest.shm_input_size`.
	std::size_t size_bytes() const {
		return size_;
	}

	/// Read/write pointer into the mmap'd shm segment.
	void *data() {
		return data_;
	}
	const void *data() const {
		return data_;
	}

private:
	InputShmRegion(std::string name, std::size_t size, void *data);
	std::string name_;
	std::size_t size_ = 0;
	void *data_ = nullptr;
};

/// Output shm segment that the daemon created and miint reads.
///
/// Lifecycle (per gpl-boundary `README.md`'s "Shared memory lifecycle" section):
///   - Daemon (gpl-boundary) creates the segment with a PID-based name and
///     reports `{shm_outputs: [{name, size, ...}]}` in the batch response.
///   - miint opens read-only with the explicit `size` from that response.
///   - **miint unlinks** (this destructor). The daemon does NOT call
///     `shm_unlink` on output segments; ownership of cleanup is explicitly
///     transferred to the consumer at the protocol level. (See
///     `gpl-boundary/README.md`: "Output shm: created by gpl-boundary
///     (PID-based name), read by miint, unlinked by miint".)
///
/// The size is supplied by the caller (from `ShmOutput.size`); the
/// implementation deliberately does NOT call `fstat` on the segment, mirroring
/// the gpl-boundary daemon-side discipline. Daemon segments are commonly
/// over-reserved relative to actual payload (see gpl-boundary's `ShmWriter`
/// 1 GiB reservation comment in `19306f6`'s commit message), so the segment's
/// physical size on disk is NOT a trustworthy proxy for payload size.
class OutputShmRegion {
public:
	/// Open an existing shm segment by name and mmap exactly `size_bytes`.
	/// Throws `std::runtime_error` on shm_open / mmap failure.
	static OutputShmRegion Open(const std::string &name, std::size_t size_bytes);

	~OutputShmRegion();

	OutputShmRegion(const OutputShmRegion &) = delete;
	OutputShmRegion &operator=(const OutputShmRegion &) = delete;
	OutputShmRegion(OutputShmRegion &&) noexcept;
	OutputShmRegion &operator=(OutputShmRegion &&) noexcept;

	const std::string &name() const {
		return name_;
	}
	std::size_t size_bytes() const {
		return size_;
	}
	const void *data() const {
		return data_;
	}

private:
	OutputShmRegion(std::string name, std::size_t size, void *data);
	std::string name_;
	std::size_t size_ = 0;
	void *data_ = nullptr;
};

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
