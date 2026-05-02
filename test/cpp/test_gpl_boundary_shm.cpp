#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gpl_boundary/shm.hpp"

#include <atomic>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

using namespace duckdb::miint::gpl_boundary;

// =============================================================================
// Cycle 2.1 — Input shm (write side, miint creates → daemon reads → miint unlinks)
// =============================================================================

TEST_CASE("InputShmRegion::Create returns a writable mmap of exactly size_bytes", "[gpl-boundary][shm]") {
	const size_t kSize = 4096;
	InputShmRegion region = InputShmRegion::Create(kSize);

	REQUIRE(region.size_bytes() == kSize);
	REQUIRE(region.data() != nullptr);

	// Write a recognizable pattern and read it back. If mmap is not actually
	// backed by the shm segment, the read-back will be junk.
	auto *p = static_cast<unsigned char *>(region.data());
	for (size_t i = 0; i < kSize; ++i) {
		p[i] = static_cast<unsigned char>(i & 0xFF);
	}
	for (size_t i = 0; i < kSize; ++i) {
		REQUIRE(p[i] == static_cast<unsigned char>(i & 0xFF));
	}
}

TEST_CASE("InputShmRegion::Create gives a name with the expected prefix", "[gpl-boundary][shm]") {
	InputShmRegion region = InputShmRegion::Create(64);
	const std::string name = region.name();
	INFO("region name: " << name);
	// gpl-boundary expects miint-side input segments to be uniquely named.
	// We use "/miint-input-<uuid>" — start with the slash, then the prefix.
	REQUIRE(name.size() > std::strlen("/miint-input-"));
	REQUIRE(name.compare(0, std::strlen("/miint-input-"), "/miint-input-") == 0);
}

TEST_CASE("InputShmRegion::Create produces unique names across instances", "[gpl-boundary][shm]") {
	InputShmRegion a = InputShmRegion::Create(64);
	InputShmRegion b = InputShmRegion::Create(64);
	REQUIRE(a.name() != b.name());
}

TEST_CASE("InputShmRegion is auto-unlinked on destruction", "[gpl-boundary][shm]") {
	std::string captured_name;
	{
		InputShmRegion region = InputShmRegion::Create(64);
		captured_name = region.name();

		// While alive, the shm segment should be openable from a second descriptor.
		const int probe = ::shm_open(captured_name.c_str(), O_RDONLY, 0);
		REQUIRE(probe >= 0);
		::close(probe);
	}
	// After destruction, opening should fail with ENOENT.
	const int probe2 = ::shm_open(captured_name.c_str(), O_RDONLY, 0);
	REQUIRE(probe2 == -1);
}

TEST_CASE("InputShmRegion is move-only and the moved-from instance does NOT unlink", "[gpl-boundary][shm]") {
	std::string name;
	{
		InputShmRegion source = InputShmRegion::Create(64);
		name = source.name();
		// Move into a new owner; the moved-from must not unlink at scope exit.
		InputShmRegion sink(std::move(source));
		// Source is now empty.
	}
	// `sink` was destroyed when its scope ended above (inner block) so the
	// segment is unlinked by `sink`. Verify.
	const int probe = ::shm_open(name.c_str(), O_RDONLY, 0);
	REQUIRE(probe == -1);
}

// =============================================================================
// Cycle 2.2 — Output shm (read side, daemon creates → miint reads → miint unlinks)
//
// Critical invariant from gpl-boundary 19306f6: explicit byte counts are the
// single source of truth on BOTH sides. Output-side mmap must NEVER call fstat
// to derive size — pass the size from the daemon's `ShmOutput.size`.
// =============================================================================

namespace {
// Helper: create an output-style shm segment OUTSIDE of OutputShmRegion (i.e.,
// stand in for the daemon). Returns the segment name + the data we wrote.
struct DaemonStubSegment {
	std::string name;
	size_t size;
};

// PID-based unique name generator for test fixtures. Avoids collisions when
// the test binary runs concurrently across CI shards or when a prior crashed
// run left stale segments in /dev/shm. Mirrors gpl-boundary's own
// `gb-<pid>-<idx>-<label>` naming convention.
std::string make_test_output_name(const char *label) {
	static std::atomic<unsigned> counter {0};
	const unsigned n = counter.fetch_add(1, std::memory_order_relaxed);
	std::string name = "/miint-test-output-";
	name += std::to_string(static_cast<unsigned long>(::getpid()));
	name += "-";
	name += std::to_string(n);
	name += "-";
	name += label;
	return name;
}

DaemonStubSegment fake_daemon_segment(const std::string &name, const std::string &payload) {
	// O_EXCL: fail loudly on a stale segment instead of silently truncating
	// over it. If we ever see EEXIST in CI it means a prior test run crashed
	// without cleanup AND reused the same pid+counter — treat it as a real
	// bug, not a transient.
	const int fd = ::shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
	if (fd < 0) {
		throw std::runtime_error("fake_daemon_segment shm_open failed for " + name + ": " +
		                         std::string(::strerror(errno)));
	}
	if (::ftruncate(fd, static_cast<off_t>(payload.size())) != 0) {
		::close(fd);
		::shm_unlink(name.c_str());
		throw std::runtime_error("fake_daemon_segment ftruncate failed");
	}
	void *p = ::mmap(nullptr, payload.size(), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	::close(fd);
	if (p == MAP_FAILED) {
		::shm_unlink(name.c_str());
		throw std::runtime_error("fake_daemon_segment mmap failed");
	}
	std::memcpy(p, payload.data(), payload.size());
	::munmap(p, payload.size());
	return {name, payload.size()};
}
} // namespace

TEST_CASE("OutputShmRegion::Open reads exactly the bytes the daemon wrote", "[gpl-boundary][shm]") {
	const std::string name = make_test_output_name("basic");
	const std::string payload = "hello arrow ipc, you arrived in shared memory";
	auto stub = fake_daemon_segment(name, payload);

	OutputShmRegion region = OutputShmRegion::Open(stub.name, stub.size);
	REQUIRE(region.size_bytes() == stub.size);
	REQUIRE(region.data() != nullptr);
	REQUIRE(std::memcmp(region.data(), payload.data(), payload.size()) == 0);
}

TEST_CASE("OutputShmRegion is auto-unlinked on destruction", "[gpl-boundary][shm]") {
	const std::string name = make_test_output_name("unlink");
	auto stub = fake_daemon_segment(name, "x");
	{
		OutputShmRegion region = OutputShmRegion::Open(stub.name, stub.size);
		(void)region;
	}
	// After destruction, opening should fail.
	const int probe = ::shm_open(name.c_str(), O_RDONLY, 0);
	if (probe >= 0) {
		::close(probe);
		::shm_unlink(name.c_str());
		FAIL("OutputShmRegion did not unlink on destruction");
	}
	REQUIRE(probe == -1);
}

TEST_CASE("OutputShmRegion::Open throws when the named segment does not exist", "[gpl-boundary][shm]") {
	REQUIRE_THROWS(OutputShmRegion::Open("/miint-test-output-definitely-not-exists-zxqv", 64));
}

TEST_CASE("OutputShmRegion::Open accepts size smaller than segment "
          "(daemon segments may be over-reserved)",
          "[gpl-boundary][shm]") {
	// gpl-boundary's ShmWriter sometimes reserves a 1 GiB segment but only
	// uses N bytes; the daemon reports the truthful N as ShmOutput.size.
	// The miint side must trust that N and never look at the segment's full size.
	const std::string name = make_test_output_name("undersize");
	auto stub = fake_daemon_segment(name, std::string(8192, 'A'));

	// Caller asks for only the first 100 bytes.
	OutputShmRegion region = OutputShmRegion::Open(stub.name, 100);
	REQUIRE(region.size_bytes() == 100);
	for (size_t i = 0; i < 100; ++i) {
		REQUIRE(static_cast<const unsigned char *>(region.data())[i] == 'A');
	}
}

TEST_CASE("InputShmRegion move-assignment transfers ownership exactly once", "[gpl-boundary][shm]") {
	std::string old_name;
	std::string new_name;
	{
		InputShmRegion sink = InputShmRegion::Create(64);
		old_name = sink.name();
		// Move-assign a different region into sink. The previous segment
		// (old_name) must be unlinked exactly once; the new segment's
		// ownership must transfer cleanly.
		InputShmRegion source = InputShmRegion::Create(128);
		new_name = source.name();
		REQUIRE(old_name != new_name);
		sink = std::move(source);
		REQUIRE(sink.size_bytes() == 128);
		REQUIRE(sink.name() == new_name);
		// old_name must be unlinked already (sink lost it on move-assign).
		const int probe = ::shm_open(old_name.c_str(), O_RDONLY, 0);
		if (probe >= 0) {
			::close(probe);
			::shm_unlink(old_name.c_str());
			FAIL("move-assign did not unlink the previous segment");
		}
		// new_name should still be alive while sink holds it.
		const int probe2 = ::shm_open(new_name.c_str(), O_RDONLY, 0);
		REQUIRE(probe2 >= 0);
		::close(probe2);
	}
	// After scope: sink is destroyed; new_name should be unlinked too.
	const int probe3 = ::shm_open(new_name.c_str(), O_RDONLY, 0);
	if (probe3 >= 0) {
		::close(probe3);
		::shm_unlink(new_name.c_str());
		FAIL("destructor of move-assigned region did not unlink the new segment");
	}
}

TEST_CASE("OutputShmRegion move-assignment transfers ownership exactly once", "[gpl-boundary][shm]") {
	const std::string name_a = make_test_output_name("mvassignA");
	const std::string name_b = make_test_output_name("mvassignB");
	auto stub_a = fake_daemon_segment(name_a, std::string(64, 'a'));
	auto stub_b = fake_daemon_segment(name_b, std::string(64, 'b'));

	{
		OutputShmRegion sink = OutputShmRegion::Open(stub_a.name, stub_a.size);
		OutputShmRegion source = OutputShmRegion::Open(stub_b.name, stub_b.size);
		sink = std::move(source);
		REQUIRE(sink.name() == name_b);
		REQUIRE(sink.size_bytes() == 64);
		// name_a (the old sink content) must be unlinked already.
		const int probe = ::shm_open(name_a.c_str(), O_RDONLY, 0);
		if (probe >= 0) {
			::close(probe);
			::shm_unlink(name_a.c_str());
			FAIL("output move-assign did not unlink the previous segment");
		}
	}
	// After scope: name_b must also be unlinked (sink owned it on exit).
	const int probe2 = ::shm_open(name_b.c_str(), O_RDONLY, 0);
	if (probe2 >= 0) {
		::close(probe2);
		::shm_unlink(name_b.c_str());
		FAIL("destructor of move-assigned output region did not unlink");
	}
}

TEST_CASE("InputShmRegion::Create rejects size 0", "[gpl-boundary][shm]") {
	REQUIRE_THROWS(InputShmRegion::Create(0));
}
