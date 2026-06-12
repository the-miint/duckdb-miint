#include <catch2/catch_test_macros.hpp>

#include "align_bowtie2_sharded.hpp"

#include <filesystem>
#include <fstream>
#include <string>
#include <unistd.h>

using duckdb::EffectiveShardThreads;
using duckdb::idx_t;

// EffectiveShardThreads decides bowtie2's `-p` (nthreads) for one shard's batch
// submit. The intent: as worker threads finish the tail of small shards and
// exit, the survivors should grow their thread count to fill the cores those
// workers freed — but never before the work is fully distributed (else they
// oversubscribe), and never below the user's configured floor. These tests
// encode WHY each rule exists; flip any rule and one assertion must fail.

TEST_CASE("EffectiveShardThreads: stays at base while shards still being handed out", "[bowtie2_sharded]") {
	// While work is still being distributed, bumping a worker would oversubscribe
	// the box against the workers about to claim the remaining shards. So the
	// gate (all_shards_claimed=false) pins everyone to base, regardless of how
	// few workers currently hold a shard.
	REQUIRE(EffectiveShardThreads(/*base*/ 4, /*db*/ 8, /*active*/ 1, /*all_claimed*/ false) == 4);
	REQUIRE(EffectiveShardThreads(4, 8, 2, false) == 4);
}

TEST_CASE("EffectiveShardThreads: sole survivor takes all cores once work is distributed", "[bowtie2_sharded]") {
	// The core win: a lone big shard left running after the tail drains should
	// use the whole box instead of idling half of it at the base `-p`.
	REQUIRE(EffectiveShardThreads(4, 8, 1, true) == 8);
}

TEST_CASE("EffectiveShardThreads: multiple survivors keep base (cores already full)", "[bowtie2_sharded]") {
	// Two survivors at base already saturate an 8-core box (2*4=8); bumping
	// either would oversubscribe. Fair share (8/2=4) equals base → no change.
	REQUIRE(EffectiveShardThreads(4, 8, 2, true) == 4);
}

TEST_CASE("EffectiveShardThreads: never drops below the configured floor", "[bowtie2_sharded]") {
	// More live workers than cores can back at base would give a fair share
	// below base (8/4=2), but base is the user's floor — honor it, don't shrink.
	REQUIRE(EffectiveShardThreads(4, 8, 4, true) == 4);
	// Floor also wins when the box has fewer cores than the configured `-p`
	// (e.g. SET threads=1 with max_threads_per_shard=4): keep the user's value.
	REQUIRE(EffectiveShardThreads(4, 1, 1, true) == 4);
}

TEST_CASE("EffectiveShardThreads: ramps with a finer base as survivors drop", "[bowtie2_sharded]") {
	// With base=2 on an 8-core box the ramp is visible at >1 survivor:
	REQUIRE(EffectiveShardThreads(2, 8, 4, true) == 2); // 8/4=2 == base
	REQUIRE(EffectiveShardThreads(2, 8, 2, true) == 4); // 8/2=4 > base
	REQUIRE(EffectiveShardThreads(2, 8, 1, true) == 8); // sole survivor
}

TEST_CASE("EffectiveShardThreads: zero active workers is a safe no-op", "[bowtie2_sharded]") {
	// A transient read where no worker holds a shard must not divide by zero.
	REQUIRE(EffectiveShardThreads(4, 8, 0, true) == 4);
}

// ShardIndexFiles drives prefetch: it must enumerate exactly the bowtie2 index
// files present for a prefix so we warm the right pages (and only those) into
// cache. The WHY: warming a wrong/empty set wastes the prefetch; warming a
// partial set is fine (the missing files just aren't there to read). These
// cases pin small vs large format, partial, and absent.
TEST_CASE("ShardIndexFiles enumerates the present bowtie2 index files", "[bowtie2_sharded]") {
	namespace fs = std::filesystem;
	const fs::path dir = fs::temp_directory_path() / ("miint-bt2-prefetch-" + std::to_string(::getpid()));
	fs::create_directories(dir);
	auto touch = [](const std::string &p) {
		std::ofstream(p) << "x";
	};

	// Small (.bt2) format: all four mandatory files present.
	const std::string small = (dir / "small").string();
	for (const char *e : {".1.bt2", ".2.bt2", ".rev.1.bt2", ".rev.2.bt2"}) {
		touch(small + e);
	}
	REQUIRE(duckdb::ShardIndexFiles(small).size() == 4);

	// Large (.bt2l) format: the variant used for big references.
	const std::string large = (dir / "large").string();
	for (const char *e : {".1.bt2l", ".2.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"}) {
		touch(large + e);
	}
	REQUIRE(duckdb::ShardIndexFiles(large).size() == 4);

	// Partial index: only the files that exist are listed (warming a subset is
	// harmless — the absent files simply aren't read).
	const std::string partial = (dir / "partial").string();
	touch(partial + ".1.bt2");
	touch(partial + ".rev.1.bt2");
	REQUIRE(duckdb::ShardIndexFiles(partial).size() == 2);

	// Absent index: nothing to warm.
	REQUIRE(duckdb::ShardIndexFiles((dir / "missing").string()).empty());

	fs::remove_all(dir);
}
