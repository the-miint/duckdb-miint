#include <catch2/catch_test_macros.hpp>

#include "align_bowtie2_sharded.hpp"

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

// -----------------------------------------------------------------------------
// Multi-thread / multi-shard correctness: verified by a local (untracked)
// harness, not a unit test. WHY there is no automated case here:
//
// The load-bearing concurrency invariant of the async pipeline is that no two
// daemon fingerprints (tool + config, where config carries bowtie2 `-p`) are
// ever outstanding on one connection at once — Session::AwaitNext correlates
// responses POSITIONALLY (FIFO), so a fingerprint change (shard switch OR a `-p`
// ramp step) must drain inflight() to 0 first. That logic only engages with
// >=2 worker threads each holding a real shard and a real gpl-boundary daemon
// fanning out on `-p`. Reproducing it in-process needs a live daemon AND
// on-disk bowtie2 indexes, neither of which exists in CI — a gated test would
// be inert everywhere it ran.
//
// So this is verified the same way the shipped `-p` ramp is: the DEPTH-INVARIANCE
// HARNESS at `scratch/bt2bench/depth_invariance_check.sh`. It aligns the same
// multi-shard read set three ways —
//     ground = threads=1 depth=1   (serial reference)
//     pipe1  = threads=N depth=1   (parallel, synchronous)
//     pipe2  = threads=N depth=2   (parallel, pipelined)
// — and asserts pipe1 and pipe2 are MULTISET-equal to ground (EXCEPT ALL both
// directions == 0) over the full alignment-identifying column set incl.
// `shard_name`, so a mis-attributed shard or a reordered/dropped row surfaces as
// a diff. Last run: 7,127,029 rows, 0 diff in every direction. sqllogictest
// can't stand in for it: it forces threads=1, collapsing the very interleaving
// this guards. The script lives under the untracked `scratch/` dir; the recipe
// to recreate it (and the bench data it needs) is in the bowtie2-sharded
// throughput notes (localdocs).
// -----------------------------------------------------------------------------
