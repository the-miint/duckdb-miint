# DuckDB engine notes (developer reference)

Non-obvious DuckDB engine behavior that affects how we build and tune this extension.
Only surprising, easily-missed, or consequential items are recorded here — not a general
DuckDB tutorial.

**Sources:** the DuckDB internals docs (`duckdb.org/docs/current/internals/`) and the
"Design and Implementation of DuckDB Internals" course (Torsten Grust, U Tübingen,
`github.com/DBatUTuebingen/DiDi`, combined deck at `blobs.duckdb.org/slides/DiDi.pdf`).
Version-specific claims below are as of DuckDB v1.4/v1.5 — re-verify when we bump.

## Storage

**New database files default to storage version 64 (= v1.0.0), even on v1.5.**
This is the single most surprising item here. Newer compression algorithms are gated
behind newer storage versions, so a database written with defaults silently forgoes
them. Opt in explicitly:

```sql
ATTACH 'file.db' (STORAGE_VERSION 'latest');   -- or 'v1.4.0', etc.
SET storage_compatibility_version = 'latest';  -- CLI/other clients
```
```bash
duckdb -storage-version latest my.db
```

Check what a file actually is: `SELECT database_name, tags FROM duckdb_databases();`
Version map: 68 = v1.5.x, 67 = v1.4.x, 66 = v1.3.x, 65 = v1.2.x, 64 = v0.9.x–v1.1.x.

The first 20 bytes of a `.db` file are `uint64` checksum + `DUCK` magic + `uint64`
version — readable without DuckDB, which is handy for test fixtures and CI assertions.

**Compression is off for in-memory databases.** It is applied to persistent databases
only unless you `ATTACH` with `COMPRESS`. In-memory benchmarks therefore overstate the
footprint of sequence data relative to what an on-disk database would use.

**DuckDB's format is not automatically smaller than Parquet.** The docs' own rule of
thumb: 100 GB of CSV → ~25 GB of DuckDB, but 100 GB of Parquet → ~120 GB of DuckDB.
Relevant whenever we advise users on where to land miint output.

Available codecs: Constant, RLE, Bit Packing, FOR, Dictionary, FSST, ALP, Chimp, Patas,
Zstd. FSST is the one that matters for sequence and read-ID strings; ALP for floats.

**Row groups default to 122,880 rows** (`ATTACH ... (ROW_GROUP_SIZE n)` to change). A row
group is simultaneously the unit of parallelism, of compression, and of zonemap pruning.

**Two processes must not share a `temp_directory`.** Spill files are named
`duckdb_temp_storage_<size>-<index>.tmp`, with no process identifier, and a database's
teardown deletes every file beginning `duckdb_temp_` in its temp directory, or the whole
directory if it created it (`duckdb/src/storage/temporary_file_manager.cpp`). An in-memory
database defaults `temp_directory` to `.tmp` under the working directory
(`DBConfig::SetDefaultTempDirectory`), so two concurrent processes started in the same
directory write the same files at the same offsets. Measured on v1.5.4 with a forced spill
and a shared working directory: most runs segfaulted, and some reported `Invalid unicode
(byte sequence mismatch) detected in segment statistics update`, which was another process's
bytes rather than invalid input. A per-process `temp_directory` removed every failure. miint
never sets `temp_directory`; a caller running concurrent processes sets one per process.

## Execution format (vectors)

**Input vectors are not necessarily flat.** The storage layer emits CONSTANT vectors
(from constant compression), DICTIONARY vectors (from dictionary compression), and
SEQUENCE vectors (for row IDs) directly into execution. Calling `FlatVector::GetData()`
on an *input* vector is a silent correctness bug, not a performance issue. Always go
through `ToUnifiedFormat()` / `UnifiedVectorFormat`, or use the `UnaryExecutor` /
`BinaryExecutor` helpers which do it for you.

Our code already follows this — `FlatVector::GetData()` appears only on output vectors.
Preserve the convention in new scalar functions.

**`STANDARD_VECTOR_SIZE` is 2048 by default but is a compile-time setting.** Never
hardcode 2048; reference the macro.

**`string_t` inlines only ≤ 12 bytes** (4-byte length + 4-byte prefix + 8-byte pointer,
union'd with 4-byte length + 12 inline bytes). Consequences for us:
- Every read sequence and quality string is heap-allocated with pointer chasing; none are
  ever inlined.
- The 4-byte prefix exists as a comparison early-out, and our data defeats it: read IDs
  sharing an instrument/run prefix (`A00123:...`) collide on all 4 prefix bytes, so every
  comparison follows the pointer anyway. Sorting or joining on read ID is more expensive
  than the same operation on a synthetic integer key.

**A NULL `string_t` slot has no defined contents.** Check validity before `GetString()` or
`GetData<string_t>`. Observed: the slot under a NULL read back from a materialized table
yields `""`, while an inline constant-NULL vector yields bytes that trip an unrelated error.
A test that feeds NULL only one way can pass against the bug, so cover both a NULL stored in
a table and a literal `NULL`.

A **morsel is 122,880 rows (60 data chunks)** — the granularity a thread claims from a
source.

## Arrow export (ArrowAppender)

**`ArrowBuffer` capacity is always a power of two, and it grows by `realloc`.**
`ArrowBuffer::reserve` does `NextPowerOfTwo(bytes)` before allocating
(`duckdb/src/include/duckdb/common/arrow/arrow_buffer.hpp`). So a finished Arrow batch
holding *N* bytes of payload occupies `NextPowerOfTwo(N)` — up to 2× the payload, and
`realloc` may need both blocks live while it moves. Two consequences:

- **Never let an Arrow batch's byte size scale with the input.** A batch sized by *rows*
  has an unbounded byte size, so the slack is an unbounded fraction of the corpus. Size
  by bytes against a power-of-two ceiling instead: a batch that stops just short of 2^k
  has capacity exactly 2^k, so the slack is one row rather than one doubling. That is
  what `RYPE_ARROW_BATCH_BYTES` in `src/include/rype_input_stream.hpp` does, and it is
  most of the fix for the-miint/Qiita#459.
- The buffers use plain `malloc`/`realloc`, not DuckDB's allocator, so none of this is
  visible to `memory_limit` or to `duckdb_memory()`. Peak RSS on an Arrow-exporting path
  is not bounded by `SET memory_limit` — measured directly: 10.2 GB peak at a 2 GB limit
  vs 11.0 GB at 24 GB, for the same query.

**The appender's `initial_capacity` argument is a row count that also sizes the data
buffer.** `ArrowVarcharData::Initialize` reserves `(capacity + 1) * sizeof(int64_t)` for
offsets *and* `capacity` bytes for the payload, each rounded up to a power of two. Passing
a large row hint to an appender that will be closed early therefore over-reserves the
offsets buffer by more than the batch ever uses.

## Pipelines and parallelism

Execution is **push-based**, not the Volcano pull model (changed in 0.3.0, Oct 2021).
Sources and sinks are parallelism-aware; the operators between them need not be. Our
table functions are *sources*, which is why they own partitioning themselves rather than
inheriting it — see the File Reading Pattern and `MaxThreads()` in `architecture.md`.

The pipeline driver already **buffers rows when a FILTER or PROBE emits very few**
(DuckDB avoids pushing chunks under ~64 rows). Don't add our own micro-batching to solve
that particular problem; it is handled upstream.

Sinks (pipeline breakers) run Sink → Combine → Finalize: `HASH_GROUP_BY`, the build side
of `HASH_JOIN`, `ORDER_BY`, `UNGROUPED_AGGREGATE`.

`PERFECT_HASH_GROUP_BY` replaces `HASH_GROUP_BY` when the grouping key provably has few
distinct values — which is the common case for our `sample_id`, SAM flag, and strand
groupings.

DuckDB's sort-key normalization is exposed at SQL level as `create_sort_key()`.

### Worker threads get a small stack — never recurse per input element

DuckDB worker threads are created with a **544 KB stack**; the main thread has the usual
8 MB. Whether a given table function body runs on a worker or on the main thread depends
on how the scheduler happens to place the pipeline, so an algorithm whose stack depth
scales with input size **crashes nondeterministically on identical input** — and it
crashes with a signal (SIGBUS/SIGSEGV), not a catchable exception, so the query does not
fail with an error: a multi-statement script simply dies partway, leaving a
partially-populated database and no message. That makes it far worse than an error.

Consequently, anything whose depth is driven by data — tree/graph traversal, nested
parsing, divide-and-conquer — must use an explicit heap-allocated work stack. `Newick`
parsing and serialization (`src/NewickTree.cpp`) and all `NewickTree` traversals
(`postorder()`, `preorder()`) are iterative for exactly this reason; the parser was
originally recursive and crashed on ordinary fragment-insertion phylogenies, which nest
~8,600 levels deep (issue #249).

A Catch2 test does **not** catch this by default: it runs on the main thread's 8 MB
stack. Tests for depth-sensitive code must pin a small stack explicitly via
`pthread_attr_setstacksize` (`std::thread` gives no control over stack size, and on Linux
inherits the 8 MB default) — see the `[depth]`-tagged cases in
`test/cpp/test_NewickParser.cpp`.

### Progress reporting

The progress bar estimates completion with an unscented Kalman filter
(`duckdb/src/common/progress_bar/unscented_kalman_filter.cpp`). An operator that reports no
progress for several seconds can drive the filter's innovation covariance toward zero; its
inverse is then infinite and the process dies with SIGFPE. Observed on an INSERT operator
whose Finalize makes several HTTPS round trips, and absent with the progress bar disabled.
`GetSinkProgress` returning `ProgressData::SetInvalid()` keeps such an operator out of the
estimate (`src/include/ena_object_insert_op.hpp`). A table function without a
`table_scan_progress` callback already reports invalid progress
(`duckdb/src/execution/operator/scan/physical_table_scan.cpp`), so a table function that
blocks in `InitGlobal` needs nothing added. An iterative fit that stops early has no known
total, so its iteration count is not a progress fraction.

## Optimizer

**Passes run once each, in a fixed order, never to fixpoint.** v1.5 has 30+ passes with a
total budget around 1 ms (10–50 µs per pass). The design explicitly trades missed
optimizations for predictable planning cost, and the pass-ordering problem is real: an
opportunity exposed by a late pass is simply not revisited.

Practical consequence: the optimizer will not clean up after messy generated SQL. Macros
and generated queries in this repo should emit reasonable plans on their own rather than
relying on the optimizer to rescue them.

Inspection and control:
```sql
SELECT * FROM duckdb_optimizers();
SET disabled_optimizers = 'filter_pullup,join_order';
PRAGMA disable_optimizer;
PRAGMA explain_output = 'optimized_only';
```

**Some optimizations are invisible to plain `EXPLAIN`.** `join_filter_pushdown` (the
build side's min/max, or an `IN` list when there are ≤ 50 distinct values, used to
pre-filter the probe side) and `row_group_pruner` are runtime effects — only
`EXPLAIN ANALYZE` reveals them. A plan can look worse than it actually runs.

**A cast on the column side of a predicate forfeits pushdown.** `WHERE id::VARCHAR IN (…)`
on a BIGINT column is a filter over an expression: `EXPLAIN` shows
`(CAST(id AS VARCHAR) IN (...))` rather than a table filter on `id`, and every row is
evaluated. Measured on a 5.2M-row table, a 1,000-id slice cost 86× the instructions of the
native predicate. On a VARCHAR column the cast is elided, so the cost appears only once ids
are integers. Compare in the column's type and cast the literals instead:
`IdInPredicate` in `src/unifrac_function_common.cpp` does this for integer and UUID ids,
whose text round trip is exact, and keeps the cast for FLOAT, DOUBLE and DECIMAL ids, whose
text form is not.

**Validation in a table macro goes in one CASE, evaluated before any filter can drop a row.**
Two engine behaviours decide this. The terms of `cond AND check` need not run in the order
written: `filter_pushdown` rebuilds the conjunction, and over both a table function and a
native table has put `check` first (observed on DuckDB v1.5.5), so a cast or `error()` that
`cond` was meant to guard runs on rows `cond` excludes (`reorder_filter`, below, is not what
moves it). A CASE tests its WHENs in order. And a check written in the SELECT list is pruned
when the query does not select that column, while a WHERE, HAVING or QUALIFY clause is
evaluated whichever columns the query selects. The plain form is therefore `WHERE CASE WHEN
<excluded> THEN FALSE WHEN <bad> THEN error(...) ELSE TRUE END` (or the same CASE in HAVING or
QUALIFY), with any cast placed after the WHEN that proves it safe.

The plain form sees only the rows that reach its clause. A query that stops early, such as one
with `LIMIT`, does not read the rest of the file. And the same pushdown can put a caller's own
filter ahead of the macro's WHERE, so a malformed record that filter excludes raises nothing:
an ungated `read_gff` did under `WHERE type IN (...)`. The gated form closes that gap. Compute
the CASE once as a column in a subquery (`... END AS keep`), wrap every output column in
`CASE WHEN keep THEN ... END`, and filter `WHERE keep`. In `EXPLAIN`, a caller's filter on a
gated column stays above the projection that computes `keep`, so the check runs on every row
the scan reads. `WHERE keep` is still needed: it drops the rows the CASE excludes, which would
otherwise come back as rows of NULLs, and it runs the check for a query that reads no column.
`read_gff` uses the gated form. On a 2-million-line, 189 MB GFF (8 threads, the file in the
page cache, five interleaved runs each) its median times stayed within 2 ms of the ungated
macro's: about 70 ms for a count, a GROUP BY or an IN filter, and about 310 ms with the
attributes parsed. The other macros in `src/include/miint_macros.hpp` use the plain form.

Passes worth knowing about for our workloads:
- `late_materialization` — keeps rows narrow through the heavy operators, then joins wide
  payload columns back by row ID. This is exactly the shape our wide records
  (sequence + quality + tags) want, and it happens automatically.
- `sum_rewriter` — `sum(e + k)` → `sum(e) + k * count(e)`, integers only (FP `+` is not
  associative).
- `window_self_join` — only for aggregate window functions whose frame is specified
  exclusively by `PARTITION BY`.
- `reorder_filter` — cheap conjuncts first, with text expressions charged a +200 cost
  penalty. Conjuncts that can throw are not reordered (that would hide a side effect), but
  `filter_pushdown` may already have reordered them; see the validation paragraph above.
- `in_clause` — a large `IN (…)` becomes a MARK join against a `COLUMN_DATA_SCAN`.

## Indexes

**Zonemaps are free and automatic** — min/max per column per row group, always present,
used for predicate pushdown into `SEQ_SCAN`. Their effectiveness depends entirely on
*physical row order*: an unordered column produces wide, overlapping spans that prune
nothing.

This is a real free win for us. Coordinate-sorted BAM input written out in input order
gives `reference_name` / `position` predicates useful row-group pruning at no cost —
another reason to preserve input ordering where we already do.

**ART indexes are used far more narrowly than PostgreSQL habits suggest.** As of v1.4,
only a single-column predicate against a single-column index is served by an index scan,
and only when selective. All of the following are *not* used:

| Predicate | Index | Used? |
|---|---|---|
| `c1 <op> val` | `t(c1)` | yes, if selective |
| `c1 <op> val` | `t(c1,c2)` composite | no |
| `c1 <op> v1 AND c2 <op> v2` | `t(c1,c2)` composite | no |
| `c1 <op> v1 AND c2 <op> v2` | `t(c1)` and/or `t(c2)` | no |
| `c1 <op> v1 OR c2 <op> v2` | `t(c1)` and `t(c2)` | no |
| `c1 LIKE 'prefix%'` | `t(c1)` | no |

**ART memory cannot be evicted.** As of v1.4 it is accounted by the buffer manager but
pinned for the session. A `CREATE INDEX` on a large miint table is a permanent memory
cost, so indexes should not be offered as a general-purpose speedup in our docs.

## Allocator

Since **v1.5.3 jemalloc ships as a third-party library, not an extension**. It is
statically linked and cannot be loaded at runtime. Linux builds include it by default;
macOS requires `BUILD_JEMALLOC=1`; Windows has no jemalloc.

Two tuning knobs worth knowing:
- `DUCKDB_JE_MALLOC_CONF` — DuckDB's rename of jemalloc's `MALLOC_CONF`, renamed to avoid
  clashing with the host process's own allocator config.
- `SET allocator_background_threads = true` — off by default; purges allocations
  asynchronously instead of on the foreground thread. The docs call out allocation-heavy
  workloads on many-core CPUs.

**Caveat that limits both:** DuckDB's jemalloc symbols are prefixed (`duckdb_je_`), so
these knobs reach only allocations that actually route through it. Extension STL and most
embedded libraries use the system allocator, so `allocator_background_threads` is not a
general fix for our allocation-bound paths (tree build/shear, alignment). minimap2 is the
deliberate exception — see the jemalloc integration notes in `embedded-tools.md`, which
also cover why the loadable-extension build falls back to system malloc. Never link a
second allocator into the process.

**Freed memory isn't returned to the OS until a thread idles.** With
`allocator_background_threads` off (the default), the *only* place DuckDB automatically
flushes a thread's tcache and purges its jemalloc arena is
`TaskScheduler::ExecuteForever`'s idle-timeout path (`duckdb/src/parallel/task_scheduler.cpp`)
— triggered after a worker sits idle on the task queue for 500ms, and once more at thread
exit. A thread that frees a large amount of memory and then either keeps working (never
idles) or was never a `TaskScheduler` worker to begin with (e.g. a table function's
`InitGlobal`, which runs on the query's calling thread) never hits that path: DuckDB's own
accounting correctly shows the memory as released (`duckdb_memory()`), but the physical
pages stay resident — invisible to `memory_limit` and to an OS/cgroup ceiling. A table
function with a large allocate-then-free burst outside the normal per-batch execution loop
needs to flush explicitly. Prefer `Allocator::ThreadFlush(background_threads_setting,
/*threshold=*/0, /*thread_count=*/1)` — it mirrors the scheduler's own forced flush and
purges only the calling thread's arena — over `Allocator::FlushAll()`, which purges every
arena in the process and is meant for genuinely global events (buffer-pool eviction
pressure, DB shutdown). See `FlushThisThreadsFreedMemory` in `align_common.hpp` for a
concrete use: `align_minimap2` calls it after the corpus-snapshot materialization, and both
`align_minimap2` and `align_minimap2_sharded` hand it to `Minimap2PartCursor`, which calls it
after every point where a multi-part index transition may have freed a part.
