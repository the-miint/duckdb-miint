#include <catch2/catch_test_macros.hpp>

#include "align_bowtie2_daemon_common.hpp"

using duckdb::bt2_daemon::ResolveBowtie2Nthreads;

// ResolveBowtie2Nthreads decides bowtie2's nthreads (`-p` for bowtie2-align,
// `--threads` for bowtie2-build) for the NON-sharded save_bowtie2_index /
// align_bowtie2 paths. The intent — and the regression these cases lock down —
// is that OMITTING `threads` must fall back to DuckDB's configured thread budget,
// not silently run on one core. The old code hard-defaulted to 1, so a 4-thread
// query built its index single-threaded (the qiita build_bowtie2_index report:
// gpl-boundary pegged at 100% of one core). An explicit value always wins so a
// user can still pin it. Flip either rule and one REQUIRE fails.
// (align_bowtie2_sharded has its own thread model — max_threads_per_shard ×
// active shards, see EffectiveShardThreads — and does not use this.)

TEST_CASE("ResolveBowtie2Nthreads: absent threads → DuckDB thread budget", "[bowtie2_nthreads]") {
	// The regression guard: no `threads` param on an 8-thread DuckDB must
	// build/align with 8, not 1 — exactly the case the hard-default-to-1 broke.
	REQUIRE(ResolveBowtie2Nthreads(/*threads_supplied*/ false, /*user_threads*/ 1, /*db_threads*/ 8) == 8);
	REQUIRE(ResolveBowtie2Nthreads(false, 1, 4) == 4);
}

TEST_CASE("ResolveBowtie2Nthreads: explicit threads wins over the DuckDB budget", "[bowtie2_nthreads]") {
	// A user who pins `threads := N` gets N regardless of the DuckDB thread count,
	// both below and above the budget.
	REQUIRE(ResolveBowtie2Nthreads(/*threads_supplied*/ true, /*user_threads*/ 2, /*db_threads*/ 8) == 2);
	REQUIRE(ResolveBowtie2Nthreads(true, 16, 8) == 16);
	// Explicit 1 is honored — a deliberate single-threaded build is not overridden
	// back up to the budget.
	REQUIRE(ResolveBowtie2Nthreads(true, 1, 8) == 1);
}

TEST_CASE("ResolveBowtie2Nthreads: degenerate db_threads floors at 1", "[bowtie2_nthreads]") {
	// TaskScheduler::NumberOfThreads() is >= 1 in practice, but a 0 must never
	// propagate to a 0/negative `-p`: floor it so the daemon gets a valid count.
	REQUIRE(ResolveBowtie2Nthreads(false, 1, 0) == 1);
}
