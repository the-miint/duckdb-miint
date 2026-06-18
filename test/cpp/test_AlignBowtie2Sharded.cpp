#include <catch2/catch_test_macros.hpp>

#include "align_bowtie2_sharded.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

using duckdb::EffectiveShardThreads;
using duckdb::idx_t;
using duckdb::InjectMemoryMappedDefault;

// InjectMemoryMappedDefault gives align_bowtie2_sharded its mm-off default:
// HPC telemetry measured sequential-fread index loads (memory_mapped=false) at
// 3.7x the throughput of --mm on a cold network FS, so sharded mode defaults the
// knob OFF — UNLESS the user supplied it, in which case their value (true or
// false) must survive untouched. The helper is templated over the map/value
// types (production: named_parameter_map_t + Value::BOOLEAN; test: a plain
// std::map<string,bool>) so this intent is exercised in the duckdb-free Catch2
// binary, the same reason EffectiveShardThreads/ShardIndexFiles are inline here.
// These cases encode WHY: absent must become false; a user value must be kept.
//
// NOT exercised here: the production named_parameter_map_t is case-INSENSITIVE
// (find/emplace fold case), whereas std::map is case-sensitive. That gap is
// immaterial — DuckDB lowercases unquoted named-parameter keys at parse time, so
// the key always arrives as "memory_mapped"; the case-insensitive map is just
// belt-and-braces for a quoted-identifier degenerate case, not logic this helper
// owns.
TEST_CASE("InjectMemoryMappedDefault: absent → defaults to false (mm-off)", "[bowtie2_sharded]") {
	std::map<std::string, bool> params;
	InjectMemoryMappedDefault(params, [](bool b) { return b; });
	REQUIRE(params.size() == 1);
	REQUIRE(params.count("memory_mapped") == 1);
	REQUIRE(params.at("memory_mapped") == false);
}

TEST_CASE("InjectMemoryMappedDefault: user memory_mapped=true is preserved", "[bowtie2_sharded]") {
	// A user who explicitly opts back into --mm must not be silently overridden
	// to mm-off — the default only fills the gap when the knob is absent.
	std::map<std::string, bool> params;
	params.emplace("memory_mapped", true);
	InjectMemoryMappedDefault(params, [](bool b) { return b; });
	REQUIRE(params.size() == 1);
	REQUIRE(params.at("memory_mapped") == true);
}

TEST_CASE("InjectMemoryMappedDefault: user memory_mapped=false is preserved (no duplicate)", "[bowtie2_sharded]") {
	// Same value as the default, but it must pass through the user-supplied path
	// (not be re-inserted), so the knob is never double-written.
	std::map<std::string, bool> params;
	params.emplace("memory_mapped", false);
	InjectMemoryMappedDefault(params, [](bool b) { return b; });
	REQUIRE(params.size() == 1);
	REQUIRE(params.at("memory_mapped") == false);
}

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

// FormatBatchTelemetry renders one per-batch record as a single TSV line. The
// WHOLE point of this telemetry is to split miint-side cost (fetch / encode /
// decode) from the opaque daemon round-trip (t_submit_ms) so we stop *inferring*
// where the per-shard seconds go. These tests pin: (1) the line is a single
// newline-terminated TSV row, (2) the data columns line up 1:1 with the header
// (so adding a field without updating the header — the classic schema-drift bug —
// fails the build), (3) the daemon black-box and miint-side phases are distinct
// columns, and (4) raw daemon `metrics` JSON rides in the trailing column so its
// embedded punctuation can't shift the earlier numeric columns.
static std::vector<std::string> SplitTabs(const std::string &line) {
	std::vector<std::string> out;
	std::string cur;
	for (char c : line) {
		if (c == '\n') {
			break; // stop at the first newline; the record is one line
		}
		if (c == '\t') {
			out.push_back(cur);
			cur.clear();
		} else {
			cur += c;
		}
	}
	out.push_back(cur);
	return out;
}

TEST_CASE("FormatBatchTelemetry: one newline-terminated TSV line, header-aligned", "[bowtie2_sharded]") {
	duckdb::BatchTelemetry r;
	r.wall_ms = 12.5;
	r.worker_id = 4242;
	r.shard = "shardA";
	r.batch_seq = 3;
	r.n_reads = 16384;
	r.input_bytes = 1000000;
	r.n_alignments = 50000;
	r.output_bytes = 2000000;
	r.t_open_stream_ms = 1.0;
	r.t_fetch_ms = 2.25;
	r.t_encode_ms = 3.0;
	r.t_submit_ms = 55000.0;
	r.t_decode_ms = 4.5;
	r.metrics = R"({"worker_majflt":17})";

	const std::string line = duckdb::FormatBatchTelemetry(r);

	// (1) Exactly one record: ends in a newline and contains no other newline.
	REQUIRE(line.back() == '\n');
	REQUIRE(std::count(line.begin(), line.end(), '\n') == 1);

	// (2) Data columns line up 1:1 with the header — a field added to the struct
	// but not the formatter (or vice-versa) breaks this.
	const auto cols = SplitTabs(line);
	const auto hdr = SplitTabs(duckdb::BatchTelemetryHeader());
	REQUIRE(cols.size() == hdr.size());
}

TEST_CASE("FormatBatchTelemetry: every column's value sits under its own header name (lockstep)", "[bowtie2_sharded]") {
	// Count-alignment alone misses a same-count reorder — e.g. n_reads <-> n_alignments
	// or input_bytes <-> output_bytes, the plausible copy-paste swaps. Pin every field
	// to a DISTINCT value and assert the value under each header NAME is the one that
	// field renders to, so any reorder/rename of either BatchTelemetryColumns or
	// FormatBatchTelemetry fails here instead of silently in an offline-parsed TSV.
	duckdb::BatchTelemetry r;
	r.wall_ms = 1.0;
	r.worker_id = 2;
	r.shard = "shard_v3";
	r.batch_seq = 4;
	r.n_reads = 5;
	r.input_bytes = 6;
	r.n_alignments = 7;
	r.output_bytes = 8;
	r.t_open_stream_ms = 9.0;
	r.t_fetch_ms = 10.0;
	r.t_encode_ms = 11.0;
	r.t_submit_ms = 12.0;
	r.t_decode_ms = 13.0;
	r.metrics = "m14";

	const auto cols = SplitTabs(duckdb::FormatBatchTelemetry(r));
	const auto hdr = SplitTabs(duckdb::BatchTelemetryHeader());
	REQUIRE(cols.size() == hdr.size());

	auto value_of = [&](const std::string &name) -> std::string {
		for (idx_t i = 0; i < hdr.size(); ++i) {
			if (hdr[i] == name) {
				return cols[i];
			}
		}
		FAIL("column not found in header: " << name);
		return {};
	};
	REQUIRE(value_of("wall_ms") == "1.000");
	REQUIRE(value_of("worker_id") == "2");
	REQUIRE(value_of("shard") == "shard_v3");
	REQUIRE(value_of("batch_seq") == "4");
	REQUIRE(value_of("n_reads") == "5");
	REQUIRE(value_of("input_bytes") == "6");
	REQUIRE(value_of("n_alignments") == "7");
	REQUIRE(value_of("output_bytes") == "8");
	REQUIRE(value_of("t_open_stream_ms") == "9.000");
	REQUIRE(value_of("t_fetch_ms") == "10.000");
	REQUIRE(value_of("t_encode_ms") == "11.000");
	REQUIRE(value_of("t_submit_ms") == "12.000");
	REQUIRE(value_of("t_decode_ms") == "13.000");
	REQUIRE(value_of("metrics") == "m14");
}

TEST_CASE("FormatBatchTelemetry: daemon round-trip is its own column, distinct from miint phases",
          "[bowtie2_sharded]") {
	// The split is the deliverable: t_submit_ms (the daemon black box) must be a
	// separate field from the miint-side fetch/encode/decode, else we're back to
	// inferring. Pin each phase to a distinct value and confirm all appear.
	duckdb::BatchTelemetry r;
	r.t_fetch_ms = 11.0;
	r.t_encode_ms = 22.0;
	r.t_submit_ms = 33.0;
	r.t_decode_ms = 44.0;
	const auto cols = SplitTabs(duckdb::FormatBatchTelemetry(r));
	const auto hdr = SplitTabs(duckdb::BatchTelemetryHeader());

	auto value_of = [&](const std::string &name) -> std::string {
		for (idx_t i = 0; i < hdr.size(); ++i) {
			if (hdr[i] == name) {
				return cols[i];
			}
		}
		FAIL("column not found in header: " << name);
		return {};
	};
	REQUIRE(value_of("t_fetch_ms") == "11.000");
	REQUIRE(value_of("t_encode_ms") == "22.000");
	REQUIRE(value_of("t_submit_ms") == "33.000");
	REQUIRE(value_of("t_decode_ms") == "44.000");
}

TEST_CASE("FormatBatchTelemetry: raw daemon metrics ride in the trailing column", "[bowtie2_sharded]") {
	// Folding the daemon's `metrics` object in as raw JSON keeps miint decoupled
	// from the (still-being-specced) metrics schema. It must be the LAST column so
	// its braces/commas can never be mistaken for additional TSV fields.
	duckdb::BatchTelemetry r;
	r.shard = "s";
	r.metrics = R"({"worker_majflt":17,"worker_reused":true})";
	const std::string line = duckdb::FormatBatchTelemetry(r);
	const auto hdr = SplitTabs(duckdb::BatchTelemetryHeader());
	REQUIRE(hdr.back() == "metrics");
	// The trailing column is the verbatim JSON, immediately before the newline.
	REQUIRE(line.find(r.metrics + "\n") != std::string::npos);

	// Absent metrics (old daemons) → empty trailing column, still header-aligned.
	duckdb::BatchTelemetry empty;
	const auto cols = SplitTabs(duckdb::FormatBatchTelemetry(empty));
	REQUIRE(cols.size() == hdr.size());
	REQUIRE(cols.back().empty());
}

TEST_CASE("FormatBatchTelemetry: control chars in free-form fields can't break the TSV row", "[bowtie2_sharded]") {
	// `shard` is an arbitrary user-supplied read_to_shard value; a tab or newline
	// in it (or in daemon metrics) must NOT split a column or truncate the line —
	// otherwise one bad shard name silently corrupts the offline-parsed telemetry.
	duckdb::BatchTelemetry r;
	r.shard = "weird\tshard\nname";
	r.metrics = "{\"k\":\"a\tb\"}";
	const std::string line = duckdb::FormatBatchTelemetry(r);
	const auto hdr = SplitTabs(duckdb::BatchTelemetryHeader());

	// Still exactly one line (the embedded newline became a space).
	REQUIRE(std::count(line.begin(), line.end(), '\n') == 1);
	REQUIRE(line.back() == '\n');
	// Still header-aligned (the embedded tabs became spaces, not new columns).
	REQUIRE(SplitTabs(line).size() == hdr.size());
}
