#include <catch2/catch_test_macros.hpp>

#include "align_bowtie2_sharded.hpp"

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <unistd.h>
#include <vector>

using duckdb::EstimateBatches;
using duckdb::idx_t;
using duckdb::InjectMemoryMappedDefault;
using duckdb::IsBigShard;

// InjectMemoryMappedDefault gives align_bowtie2_sharded its mm-off default:
// HPC telemetry measured sequential-fread index loads (memory_mapped=false) at
// 3.7x the throughput of --mm on a cold network FS, so sharded mode defaults the
// knob OFF — UNLESS the user supplied it, in which case their value (true or
// false) must survive untouched. The helper is templated over the map/value
// types (production: named_parameter_map_t + Value::BOOLEAN; test: a plain
// std::map<string,bool>) so this intent is exercised in the duckdb-free Catch2
// binary, the same reason IsBigShard / ShardIndexFiles are inline here.
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

// EstimateBatches / IsBigShard drive the two-pool scheduler's partition. A shard
// is submitted to the daemon in chunks of `submit_batch_reads` reads; the number
// of submits it takes is ceil(read_count / submit_batch_reads). IsBigShard calls
// a shard "big" iff it needs >= big_shard_min_batches submits — enough reads that
// it is compute-bound and worth handing the full-thread (`-p db_threads`) big
// lane; fewer ⇒ a load-bound "tail" shard that runs single-threaded (`-p1`)
// alongside others. These tests pin the breakpoint math because mis-partitioning
// is the one way the scheduler can regress without changing the alignment output.
TEST_CASE("EstimateBatches: ceil division, exact and off-by-one boundaries", "[bowtie2_sharded]") {
	// Empty / sub-batch shards still cost one submit (ceil), except a literally
	// zero-read shard which costs none.
	REQUIRE(EstimateBatches(0, 16384) == 0);
	REQUIRE(EstimateBatches(1, 16384) == 1);
	REQUIRE(EstimateBatches(16384, 16384) == 1); // exactly one batch
	REQUIRE(EstimateBatches(16385, 16384) == 2); // one over ⇒ a second submit
	REQUIRE(EstimateBatches(32768, 16384) == 2); // exactly two
	REQUIRE(EstimateBatches(32769, 16384) == 3); // one over ⇒ a third
}

TEST_CASE("EstimateBatches: zero batch size is a safe no-op (no divide-by-zero)", "[bowtie2_sharded]") {
	// submit_batch_reads is bind-validated to >= 1, but the pure fn must not trap
	// if ever called with 0 (defensive — a div-by-zero here would crash a worker).
	REQUIRE(EstimateBatches(1000, 0) == 0);
}

TEST_CASE("IsBigShard: breakpoint is `>= min_batches` submits", "[bowtie2_sharded]") {
	// Default min_batches=2: a single-batch shard (read_count <= submit_batch_reads)
	// is tail; the first read that forces a second submit makes it big. This is the
	// exact line that sent 837/1000 single-batch shards to the load-bound tail.
	REQUIRE(IsBigShard(/*reads*/ 16384, /*sbr*/ 16384, /*min*/ 2) == false); // 1 batch  -> tail
	REQUIRE(IsBigShard(16385, 16384, 2) == true);                            // 2 batches -> big
	REQUIRE(IsBigShard(1, 16384, 2) == false);                               // tiny shard -> tail
}

TEST_CASE("IsBigShard: min_batches=1 makes every non-empty shard big", "[bowtie2_sharded]") {
	// The floor case: with min_batches=1 even a one-read shard clears the bar, so
	// the whole set goes to the big lane (sequential, full-thread). Pins that the
	// breakpoint is inclusive (`>=`), not strict.
	REQUIRE(IsBigShard(1, 16384, 1) == true);
	REQUIRE(IsBigShard(16384, 16384, 1) == true);
	// A literally zero-read shard (defensive: read_count is >=1 in practice via the
	// GROUP BY) has 0 batches, so it never clears even min_batches=1 ⇒ tail.
	REQUIRE(IsBigShard(0, 16384, 1) == false);
}

TEST_CASE("IsBigShard: higher min_batches pushes the breakpoint out", "[bowtie2_sharded]") {
	// min_batches=3 requires a third submit to qualify as big: two full batches is
	// still tail, the read that forces the third flips it. Confirms the threshold
	// tracks min_batches rather than being hard-coded at 2.
	REQUIRE(IsBigShard(32768, 16384, 3) == false); // 2 batches -> tail
	REQUIRE(IsBigShard(32769, 16384, 3) == true);  // 3 batches -> big
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
