#include <catch2/catch_test_macros.hpp>

#include "shard_progress.hpp"

#include <string>

using duckdb::shard_progress::FormatShardDone;
using duckdb::shard_progress::FormatShardStart;
using duckdb::shard_progress::GroupThousands;

// These pin the EXACT user-facing format of the opt-in `progress` lines emitted
// by align_minimap2_sharded / align_bowtie2_sharded. The shaln --verbose layer
// (and users) read/grep these, so a format change must break a test rather than
// ship silently. Emit() (timestamp + stderr write) is impure and is covered
// end-to-end by shaln's --verbose stderr-capture test instead.

TEST_CASE("GroupThousands inserts comma separators", "[shard_progress]") {
	CHECK(GroupThousands(0) == "0");
	CHECK(GroupThousands(7) == "7");
	CHECK(GroupThousands(42) == "42");
	CHECK(GroupThousands(999) == "999");
	CHECK(GroupThousands(1000) == "1,000");
	CHECK(GroupThousands(48210) == "48,210");
	CHECK(GroupThousands(1000000) == "1,000,000");
}

TEST_CASE("FormatShardStart includes the read count when known", "[shard_progress]") {
	CHECK(FormatShardStart(1, 2, "shard_a", 50000) == "shard 1/2 'shard_a': index loaded, 50,000 reads");
}

TEST_CASE("FormatShardStart omits the read count when unknown (streamed)", "[shard_progress]") {
	// The bowtie2 path streams reads, so the per-shard count is unknown at load
	// (n_reads < 0) and must be omitted from the "started" line.
	CHECK(FormatShardStart(2, 5, "shard_b", -1) == "shard 2/5 'shard_b': index loaded");
}

TEST_CASE("FormatShardDone reports reads, alignments, and elapsed seconds", "[shard_progress]") {
	CHECK(FormatShardDone(1, 2, "shard_a", 50000, 48210, 2.34) ==
	      "shard 1/2 'shard_a': done - 50,000 reads, 48,210 alignments (2.34s)");
}

TEST_CASE("FormatShardDone rounds elapsed to two decimals", "[shard_progress]") {
	CHECK(FormatShardDone(3, 3, "s", 0, 0, 1.0) == "shard 3/3 's': done - 0 reads, 0 alignments (1.00s)");
}
