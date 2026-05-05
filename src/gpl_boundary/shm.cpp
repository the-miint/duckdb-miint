#include "gpl_boundary/shm.hpp"

#include <atomic>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace duckdb {
namespace miint {
namespace gpl_boundary {

namespace {

// Unique-name generator for input segments. Per-process counter + pid avoids
// collision between concurrent shm allocations in the same DuckDB connection.
// gpl-boundary just opens whatever name we give it, so the only requirement
// is uniqueness within /dev/shm/.
std::string make_input_name() {
	static std::atomic<uint64_t> counter {0};
	const uint64_t n = counter.fetch_add(1, std::memory_order_relaxed);
	const pid_t pid = ::getpid();
	std::string name;
	name.reserve(48);
	name += "/miint-input-";
	name += std::to_string(static_cast<uint64_t>(pid));
	name += "-";
	name += std::to_string(n);
	return name;
}

} // namespace

// =============================================================================
// InputShmRegion
// =============================================================================

InputShmRegion::InputShmRegion(std::string name, std::size_t size, void *data)
    : name_(std::move(name)), size_(size), data_(data) {
}

InputShmRegion InputShmRegion::Create(std::size_t size_bytes) {
	if (size_bytes == 0) {
		throw std::runtime_error("gpl_boundary: InputShmRegion size must be > 0");
	}

	// EEXIST retry: pid recycling can collide our (pid,counter) name space
	// with leftover segments from a previously SIGKILL'd miint process. Bump
	// the counter and try again. A small bounded loop is enough; if even
	// 64 consecutive names collide, something else is broken.
	std::string name;
	int fd = -1;
	for (int attempt = 0; attempt < 64; ++attempt) {
		name = make_input_name();
		fd = ::shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
		if (fd >= 0) {
			break;
		}
		if (errno != EEXIST) {
			throw std::runtime_error("gpl_boundary: shm_open(" + name + ") failed: " + std::string(::strerror(errno)));
		}
		// EEXIST: a stale segment is squatting on this name. Re-roll.
	}
	if (fd < 0) {
		throw std::runtime_error("gpl_boundary: shm_open exhausted retries on EEXIST — /dev/shm/ may "
		                         "be full of stale miint segments. Manual cleanup required.");
	}
	if (::ftruncate(fd, static_cast<off_t>(size_bytes)) != 0) {
		const int err = errno;
		::close(fd);
		::shm_unlink(name.c_str());
		throw std::runtime_error("gpl_boundary: ftruncate(" + name + ") failed: " + std::string(::strerror(err)));
	}
	void *p = ::mmap(nullptr, size_bytes, PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	::close(fd); // mmap holds its own reference; the fd is no longer needed
	if (p == MAP_FAILED) {
		const int err = errno;
		::shm_unlink(name.c_str());
		throw std::runtime_error("gpl_boundary: mmap(" + name + ") failed: " + std::string(::strerror(err)));
	}
	return InputShmRegion(name, size_bytes, p);
}

InputShmRegion::~InputShmRegion() {
	if (data_ && size_ > 0) {
		::munmap(data_, size_);
	}
	if (!name_.empty()) {
		// Best-effort: if the daemon already unlinked it (it shouldn't on the
		// input path, but) we don't care. ENOENT here is benign.
		::shm_unlink(name_.c_str());
	}
}

InputShmRegion::InputShmRegion(InputShmRegion &&other) noexcept
    : name_(std::move(other.name_)), size_(other.size_), data_(other.data_) {
	other.name_.clear();
	other.size_ = 0;
	other.data_ = nullptr;
}

InputShmRegion &InputShmRegion::operator=(InputShmRegion &&other) noexcept {
	if (this != &other) {
		this->~InputShmRegion();
		new (this) InputShmRegion(std::move(other));
	}
	return *this;
}

// =============================================================================
// OutputShmRegion
// =============================================================================

OutputShmRegion::OutputShmRegion(std::string name, std::size_t size, void *data)
    : name_(std::move(name)), size_(size), data_(data) {
}

OutputShmRegion OutputShmRegion::Open(const std::string &name, std::size_t size_bytes) {
	if (name.empty()) {
		throw std::runtime_error("gpl_boundary: OutputShmRegion name must be non-empty");
	}
	if (size_bytes == 0) {
		throw std::runtime_error("gpl_boundary: OutputShmRegion size must be > 0 (got 0). "
		                         "The daemon's ShmOutput.size is the authoritative payload "
		                         "length; passing 0 here would yield an empty mmap.");
	}
	const int fd = ::shm_open(name.c_str(), O_RDONLY, 0);
	if (fd < 0) {
		throw std::runtime_error("gpl_boundary: shm_open(" + name +
		                         ", O_RDONLY) failed: " + std::string(::strerror(errno)));
	}
	// We deliberately do NOT call fstat to derive size. gpl-boundary commit
	// 19306f6 made explicit size authoritative on both sides; the segment's
	// physical size on disk is unrelated to the payload size.
	void *p = ::mmap(nullptr, size_bytes, PROT_READ, MAP_SHARED, fd, 0);
	::close(fd);
	if (p == MAP_FAILED) {
		const int err = errno;
		throw std::runtime_error("gpl_boundary: mmap(" + name + ") failed: " + std::string(::strerror(err)));
	}
	return OutputShmRegion(name, size_bytes, p);
}

OutputShmRegion::~OutputShmRegion() {
	if (data_ && size_ > 0) {
		::munmap(data_, size_);
	}
	if (!name_.empty()) {
		::shm_unlink(name_.c_str());
	}
}

OutputShmRegion::OutputShmRegion(OutputShmRegion &&other) noexcept
    : name_(std::move(other.name_)), size_(other.size_), data_(other.data_) {
	other.name_.clear();
	other.size_ = 0;
	other.data_ = nullptr;
}

OutputShmRegion &OutputShmRegion::operator=(OutputShmRegion &&other) noexcept {
	if (this != &other) {
		this->~OutputShmRegion();
		new (this) OutputShmRegion(std::move(other));
	}
	return *this;
}

} // namespace gpl_boundary
} // namespace miint
} // namespace duckdb
