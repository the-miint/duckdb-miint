# DuckDB engine notes (developer reference)

Non-obvious DuckDB engine behavior that affects how we build and tune this extension.
Only surprising, easily-missed, or consequential items are recorded here — not a general
DuckDB tutorial.

**Sources:** the DuckDB internals docs (`duckdb.org/docs/current/internals/`) and the
"Design and Implementation of DuckDB Internals" course (Torsten Grust, U Tübingen,
`github.com/DBatUTuebingen/DiDi`, combined deck at `blobs.duckdb.org/slides/DiDi.pdf`).
Version-specific claims below are as of DuckDB v1.4/v1.5 unless they name 2.0 — re-verify when
we bump.

## Storage

**On DuckDB 1.x, new database files default to storage version 64 (= v1.0.0), even on
v1.5.** (2.0 changed this; see below.) This is the single most surprising item here.
Newer compression algorithms are gated behind newer storage versions, so a database
written with defaults silently forgoes them. Opt in explicitly:

```sql
ATTACH 'file.db' (STORAGE_VERSION 'latest');   -- or 'v1.4.0', etc.
SET storage_compatibility_version = 'latest';  -- CLI/other clients
```
```bash
duckdb -storage-version latest my.db
```

Check what a file actually is: `SELECT database_name, tags FROM duckdb_databases();`
Version map: 69 = v2.0.x, 68 = v1.5.x, 67 = v1.4.x, 66 = v1.3.x, 65 = v1.2.x, 64 = v0.9.x–v1.1.x.

The first 20 bytes of a `.db` file are `uint64` checksum + `DUCK` magic + `uint64`
version — readable without DuckDB, which is handy for test fixtures and CI assertions.
Since DuckDB 2.0 that field reads **999** (`StorageVersion::DEPRECATED`) on files
written at the v2.0 storage version; the real version moved into the database header.

**DuckDB 2.0 flips the default to the latest storage version (v2.0.0), and DuckDB 1.x
cannot open those files** ("Trying to read a database file with version number 999").
Opening and writing an *existing* older file keeps its version, so 1.x can still read
it. To create a file 1.x can read, pin the version:

```sql
ATTACH 'file.db' (STORAGE_VERSION 'v1.5.0');
```
```bash
duckdb -storage-version v1.5.0 my.db
```

`ATTACH 'old.db' (STORAGE_VERSION 'latest')` upgrades an existing file in place, and
1.x can no longer open it afterwards. All verified on DuckDB 2.0 (`96df096fb00`)
against 1.5.5.

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

**On DuckDB 2.0, `NULL AS x` in a CTAS makes a column of type `"NULL"`.** On 1.5.5 the same
column was INTEGER. A function that checks a column's type then refuses it, as
`rype_classify`, `rype_log_ratio` and `rype_filter` refuse a `sequence2` made this way;
`NULL::VARCHAR AS x` keeps a usable type. Test fixtures that build an all-NULL column cast it.

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
- Most read sequences and quality strings are longer, so they are on the heap, reached through
  the pointer. Not all: a read, a trimmed read, a barcode or a UMI can be 12 bytes or fewer, and
  is then inline. A `const char *` from `GetData()` must then come from the `string_t` in the
  vector, not from a copy of it, which takes the bytes with it (the FASTA/FASTQ writers'
  `SequenceChunkReader` returns references for this reason).
- The 4-byte prefix exists as a comparison early-out, and our data defeats it: read IDs
  sharing an instrument/run prefix (`A00123:...`) collide on all 4 prefix bytes, so every
  comparison follows the pointer anyway. Sorting or joining on read ID is more expensive
  than the same operation on a synthetic integer key.

**A NULL `string_t` slot has no defined contents.** Check validity before `GetString()` or
`GetData<string_t>`. Observed: the slot under a NULL read back from a materialized table
yields `""`, while an inline constant-NULL vector yields bytes that trip an unrelated error.
A test that feeds NULL only one way can pass against the bug, so cover both a NULL stored in
a table and a literal `NULL`.

**An all-valid validity buffer is not free.** A `ValidityMask` with no buffer already means
every row is valid; `SetAllValid()` allocates a buffer and fills it with ones. The storage
append (`validity_uncompressed.cpp`) takes its fast path only when there is no buffer, so a
column handed an all-ones buffer is appended row by row. On a list column this counts the
child too: an all-ones score mask on `read_fastx`'s quality lists sends every score down the
row-by-row path, about 3 billion of them for 10M pairs of 150-base reads. Leave a fresh output mask
alone: DuckDB resets each output chunk's masks before the source fills it
(`PipelineExecutor` calls `Reset()`, and `VectorCache::ResetFromCache` drops the buffers).

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

### Interrupts and `max_execution_time` (DuckDB 2.0)

**`IsInterrupted()` does not see `max_execution_time`.** It reads only the interrupt flag that
Ctrl-C or `Interrupt()` sets. The deadline is enforced by `ClientContext::InterruptCheck()`,
which reads the clock on every 256th call on a thread (a per-thread count, kept across
queries), by `BlockingInterruptCheck()`, which reads it on every call, and between the
statements of a multi-statement query; each throws "Query exceeded maximum execution time"
(`duckdb/src/main/client_context.cpp`, `client_context_execution.cpp`). DuckDB's pipeline
executor calls `InterruptCheck()` for every chunk, which does nothing for an `InitGlobal` or
`Execute` that runs for an hour before it returns. So a long step checks for itself:
`InterruptCheck()` in a hot loop, `BlockingInterruptCheck()` between coarse steps (a batch, a
wave, one library call) and in a wait. Either may be called from a worker thread of the
function's own that is joined before the function returns: they read an atomic, the clock, and
a deadline that is set before the query runs and cleared after it ends.

A sleep is a wait too. Every sleep in the HTTP clients (a retry's backoff, a rate limit, BLAST's
waits for a search) goes through `SleepUnlessInterrupted` (`src/include/http_client_common.hpp`),
which calls `BlockingInterruptCheck()` after every slice of at most 100 ms, the last one
included, so a cancel during the wait also stops the request the client would send next.
DuckDB's own HTTP layer does not check. `HTTPUtil::RunRequestWithRetry`
(`duckdb/src/main/http/http_util.cpp`) retries a failed idempotent request up to `http_retries`
times (default 3), and one answered 429 or 503 five times more, sleeping up to 20 s between tries
(`http_retry_budget.cpp`). With the defaults, a GET to a host that keeps answering 503 holds the
query for 44 to 88 s, by the jittered schedule (68 s measured, #414), before the client sees the
status, whatever Ctrl-C does. A POST is never retried there. The clients take their HTTP
parameters from the database (`InitializeParameters(db, url)`), so `SET GLOBAL http_retries`
changes theirs and a session `SET` does not. A GET whose body streams to a handler has more to it:
[A GET that streams its body](#a-get-that-streams-its-body-duckdb-20-httpfs).

### Waiting and background work (DuckDB 2.0)

Verified on DuckDB 2.0 (`d54cbbe`). `rype_filter` uses all of it: `src/rype_filter.cpp`, with
its state machine in `src/include/rype_filter_schedule.hpp`.

**A table function can wait without holding a thread.** It parks: it sets
`async_result = AsyncResultType::BLOCKED`, keeps a copy of `*input.interrupt_state`, and returns.
Calling `Callback()` on that copy later puts the task back on the queue. The rules:

- Park only when `results_execution_mode == AsyncResultsExecutionMode::TASK_EXECUTOR` and
  `interrupt_state` is set. A synchronous scan has no task to wake, and neither has a function
  another operator drives: a positional join, or a LATERAL `in_out_function`, which gets
  neither (`physical_tableinout_function.cpp`). Such a function needs a second way to wait,
  a timed condition-variable wait that calls `BlockingInterruptCheck()` each time round.
- Register the waiter last, under the lock every wake-up takes, so no wake-up falls between
  the check and the registration. Wake by swapping the waiter list out under the lock, then
  call `Callback()` after unlocking.
- Never call back one registration twice, and never one the calling thread made in this call.
  `Executor::RescheduleTask` spins until the task has returned from its call. A thread that
  has just started background jobs therefore looks again instead of parking: starting a job
  can fail, and the failure wakes every parked thread, its own registration included.
- Set `async_result` on every return. An empty chunk left at the default (`IMPLICIT`) reads
  as FINISHED, and the scan ends with work still to emit.

**Background jobs go on the async pool.** A `TaskExecutor(context, TaskSchedulerType::ASYNC)`
(`duckdb/src/include/duckdb/parallel/task_executor.hpp`) runs `BaseExecutorTask`s that no
pipeline task waits on. In detail:

- Idle regular workers take async tasks too (`TaskScheduler::ExecuteForever`,
  `duckdb/src/parallel/task_scheduler.cpp`), so `async_threads` does not cap how many run.
- `NumberOfAsyncThreads()` is 0 in a WASM build without pthreads, and `SET async_threads = 0`
  makes it 0 anywhere. Keep an inline path that runs the job on the calling thread.
- A `ScheduleTask` that fails retires the task, which calls its `Cancel()`, and rethrows.
  `Cancel()` must change whatever the waiters watch, on every path, or they wait for a job
  that will never run (`FilterJobTask::Cancel`).
- `CancelAndDrain()` waits for running jobs to finish. Call it before freeing anything the jobs
  share, and give them that state by `shared_ptr`, so no job depends on the global state's
  lifetime (`RypeFilterTableFunction::GlobalState::~GlobalState`). The wait means a job's own
  sleeps, waits and pipe reads delay teardown unless something can cut them short.
- A library error kept in a thread-local, such as `rype_get_last_error()`, is read inside the
  job that made the call.

### COPY TO

Verified on DuckDB 2.0 (`d54cbbe`).

- **The format's `execution_mode` picks the operator.**
  - With none, the COPY is REGULAR: `PhysicalCopyToFile` is a sequential sink, so the pipeline,
    source included, is one task whatever the source's `MaxThreads()` (`plan_copy_to_file.cpp`,
    `pipeline.cpp`).
  - PARALLEL lets every thread prepare batches.
  - BATCH (`PhysicalBatchCopyToFile`) restores source order from batch indexes.
  - What miint's sequence writers (those that plan through `SequenceCopyPlan`) return is in
    [architecture.md](architecture.md#record-writers-on-duckdbs-batch-api), and so is how they turn
    order off by default: in their `plan` callback, on the bound `LogicalCopyToFile`, since the
    callback below cannot.
- **"Order kept" is what the callback is told.** It is false for a source that declares
  `NO_ORDER`; true for an `ORDER BY` (`FIXED_ORDER`) whatever the setting; otherwise the
  `preserve_insertion_order` setting decides. The COPY's `PRESERVE_ORDER` overrides all of these,
  and `PARTITION_BY`, `PER_THREAD_OUTPUT`, `FILE_SIZE_BYTES` and `ORDER_BY` force it off
  (`PhysicalPlanGenerator::PreserveInsertionOrder`, called from `plan_copy_to_file.cpp`). `EXPLAIN` does not show which mode
  ran, and the COPY's `PhysicalOperator` log entries carry no thread id; `PhysicalCopyToFile::parallel` in the plan
  does show it (`test_CopySequenceOrder.cpp` plans a COPY and reads it).
- **`ExtractPlan` then `PhysicalPlanGenerator::Plan` resolves a plan's columns twice.** `ExtractPlan`
  already runs `ColumnBindingResolver`, and `Plan` runs it again; most plans survive, but a top-N
  fails with an internal error (`Failed to bind column reference ... inequal num bindings/types`).
  To build a physical plan in a test, parse, plan and optimize the statement yourself inside
  `RunFunctionInTransaction`, then call `Plan`. v2.0's `PreparedStatement` no longer exposes its
  plan.
- **An order-kept COPY in REGULAR mode still runs its source on one thread.** That covers
  miint's sequence writers, and CSV or Parquet over a source without a batch index
  ([Output order](#output-order-costs-threads-and-memory)); over one with a batch index, CSV and
  Parquet run in BATCH mode. A table function with heavy
  per-row work does that work on its own jobs and only emits from `Execute` (`rype_filter`
  maps reads on async jobs). Time such a function through the sink its users write to: a CTAS
  of a source that declares `NO_ORDER` is a parallel sink, so it does not show the problem.
- **Batches.** DuckDB collects rows into a `ColumnDataCollection`, then calls the format's
  `prepare_batch` and `flush_batch`. A format without them gets a legacy adapter: its
  `copy_to_sink` and `copy_to_combine` run once per batch, with a fresh local state, under the
  per-file lock (`FlushLegacyCopyBatch` in `physical_copy_to_file.cpp`).
  - In REGULAR and PARALLEL, DuckDB appends whole input chunks to a batch and checks the caps
    after each append (`PhysicalCopyToFile::Sink`). With neither `ROW_GROUP_SIZE` nor `ROW_GROUP_SIZE_BYTES`,
    every chunk is a batch of its own (`CopyFunctionBatchAnalyzer::AnyBatchQualifies`).
  - A batch is flushed at whichever cap it reaches first: `ROW_GROUP_SIZE` rows or
    `ROW_GROUP_SIZE_BYTES`. The row cap never cuts a chunk, so a 2,048-row chunk stays whole
    under `ROW_GROUP_SIZE 64`.
  - Only `ROW_GROUP_SIZE_BYTES` cuts a chunk (`CopyBatchSlicer::Slice`): the batch takes rows up
    to and including the one that reaches the cap. The exception is the rest of a chunk: it is
    left whole unless the batch reaches 9/8 of the cap before the chunk's last row, so a batch
    can stay under 9/8 of the cap until its last row, which can take it further.
    `'8MB'` is 8,000,000 bytes, as `ColumnDataCollection::SizeInBytes` counts them: the rows the
    query hands the format, not what the format writes. SAM's `SEQUENCE_DATA` reads, looked up
    while a record is built, are not in them, so a capped batch can make records of any size.
  - BATCH mode is different: it holds each source batch whole before cutting it, which for a
    table or Parquet row group of genomes is most of the input.
- **Defaults.** A format can default the rows (`desired_batch_size`) but not the bytes: the
  `copy_default_batch_size_bytes_t` typedef has no field. miint's sequence writers add both
  through a `plan` callback instead (`SequenceCopyPlan`). `plan` replaces binding, so it must rebind as a different
  format, or it recurses. DuckDB checks options against the format named in the statement before
  calling `plan`. `PARTITION_PATH` survives the rebind: the rebound statement binds it as any
  COPY would (miint's sequence writers write correct partition directories with it, as DuckDB's
  JSON writer, which also rewrites its COPY through `plan`, does).
- **Concurrency.** In PARALLEL and partitioned COPYs, `prepare_batch` runs outside any lock and
  is given the first file's global state, not necessarily the destination's
  (`PhysicalCopyToFile::PrepareBatch`). `flush_batch` runs under the destination's per-file
  lock.
- **A failed COPY leaves no output.** Each file and directory a COPY creates is registered with
  its `CopyOutputLifecycle`, which removes them unless the COPY succeeds
  (`copy_output_lifecycle.cpp`); this covers `PARTITION_BY` and `PER_THREAD_OUTPUT` too. A COPY
  to one local file that already exists writes `tmp_<name>` beside it and moves it over the
  target only on success (`ResolveUseTmpFile` in `bind_copy.cpp`, `plan_copy_to_file.cpp`), so the
  old file survives a failure. The format is given the `tmp_` path, so its errors name that.
  So a format need not clean up after its own error, and an error raised after some batches were
  written leaves nothing partial.
- **Observing batches.** `CALL enable_logging('PhysicalOperator')` (in the shell, add
  `storage = 'memory'`), then read `duckdb_logs_parsed('PhysicalOperator')`. There is one
  `PhysicalCopyToFile`/`FlushBatch` row per batch, with `info.file`, `info.rows`, `info.size`
  and `info.reason` (`BATCH_SIZE`, `BATCH_SIZE_BYTES`, `FORCED_FLUSH`), but no thread id.
  `test/sql/copy_fastx_batches.test` uses it.

### Output order costs threads and memory

Verified on DuckDB 2.0 (`d54cbbe`). Under the default `preserve_insertion_order`:

- **A source that declares no order and no batch index runs its consumer on one thread.** A
  table function left at `INSERTION_ORDER` without `get_partition_data` gets a
  `SOURCE_ORDERED` result sink (`physical_result_collector.cpp`), which is not a parallel sink
  (`PhysicalResultSink::ParallelSink`), and a CTAS gets a single-threaded insert
  (`plan_create_table.cpp`). Its own threads then go unused. A function whose rows come out
  in no fixed order sets `order_preservation_type = OrderPreservationType::NO_ORDER`.
- **A batch index lets the consumer run on every thread and still keep the order.**
  - **Setting `get_partition_data` is all it takes.** DuckDB does not check what the request
    asks for (`PhysicalTableScan::SupportsPartitioning`, a FIXME).
  - **The consumers.** `PhysicalPlanGenerator::UseBatchIndex` admits a plan when there is more
    than one thread and every source has a batch index, whatever the sources' own thread counts.
    With order kept, these then use the index:
    - CREATE TABLE AS and INSERT (`PhysicalBatchInsert`; `plan_create_table.cpp`,
      `plan_insert.cpp`);
    - COPY to a format with a batch mode, CSV and Parquet (`plan_copy_to_file.cpp`);
    - a query's result, materialized or streamed (`ResultOrdering::BATCH_INDEX_ORDERED`,
      `physical_result_collector.cpp`);
    - a materialized CTE (`plan_cte.cpp`);
    - a LIMIT, with a constant limit and offset of at most 10,000 rows, over anything but a bare
      scan, which includes a table function's own (`UseBatchLimit`, `plan_limit.cpp`): so
      `FROM read_fastx(glob) LIMIT 10` gets the streaming limit;
    - a `USING SAMPLE n ROWS (system)`; the default row-count sample is a reservoir sample, which
      uses no batch index (`plan_sample.cpp`).

    `EXPLAIN` names the batch operators "Batch Create Table As", "Batch Insert" and
    "Batch Copy To File".
  - `read_fastx` declares one through `ScanBatchIndex`
    ([architecture.md](architecture.md#output-order)).
- **The batch-index contract:**
  - `get_partition_data` runs after each non-empty chunk, so the local state must still describe
    the chunk `Execute` just returned (`PipelineExecutor::NextBatch`).
  - A source's index is below `PipelineBuildState::BATCH_INCREMENT - 2` (10^13 - 2)
    (`PipelineExecutor::ToPipelinePartitionData`).
  - On a thread, indexes never decrease (`PipelineExecutor::NextBatch`). A thread starts from the
    lowest index any thread holds when it starts, and no new index may fall below the lowest one
    any thread still holds (`Pipeline::UpdateBatchIndex`). Indexes may repeat (the chunks then
    form one batch), and gaps are fine.
  - Two threads never share an index: the sinks reject it (`PhysicalBatchInsert::AddCollection`,
    `PhysicalBatchCopyToFile`).

  A source that breaks the contract fails the query with an `InternalException`, which in
  DuckDB 2.0 invalidates the database. A source whose indexes rise with the work units it hands
  out keeps it; `ScanBatchIndex` does so for readers that claim whole files.
- **Order-keeping sinks hold every batch until enough follow it.** `PhysicalBatchInsert`
  keeps each batch unmerged, at a cost of its own even when it holds a few rows, until
  consecutive batches add up to a row group. A batch index per small work unit made 60k rows
  in 7-row batches run out of memory under a 64 MB limit at 8 threads. Each thread holds whole
  batches, so their memory also follows a batch's bytes: a reader of large records caps a
  chunk's bytes, as `read_fastx`'s `max_batch_bytes` does.
- **A reader that gives each thread a whole file makes the sink hold most of the input.**
  - Batches from later files wait until the file with the lowest index is done. For a few large
    files that is nearly everything.
  - The batch sinks' `BatchMemoryManager` reserves memory for them, at most a quarter of the
    operator memory limit at a time (`BatchMemoryManager::SetMemorySize`). Past its reservation,
    it blocks every task whose batch index is above the lowest until that index advances, trading
    speed for memory. RSS counts more than the reservation (the readers' buffers, the output), so
    it can pass `memory_limit`.
  - Measured with `read_fastx` on DuckDB 2.0, 8 threads pinned to 8 cores, medians of 3 runs,
    peak RSS from `/usr/bin/time`, for
    `COPY (SELECT read_id, sequence1 FROM read_fastx(<eight 75 MB FASTQ files>)) TO '/dev/null' (FORMAT CSV)`:
    1.1 GB and 0.10 s under an 8 GB `memory_limit`; 425 MB and 0.61 s under 256 MB; and before
    `read_fastx` had a batch index, on one thread, 102 MB and 0.57 s.
  - DuckDB's own readers split a file across threads, which keeps that window small.
- **A one-thread source gains nothing from a batch index, and it costs a CREATE TABLE AS.**
  A one-thread source with a batch index still gets `PhysicalBatchInsert` (the consumers, above).
  - From a one-file `read_fastx` of 604 MB of short reads (8 threads pinned to 8 cores, medians
    of 3 to 5 runs), it merged one 2,048-row batch per chunk into row groups: 1.86 s, against
    1.16 s before `read_fastx` had a batch index, about 60% slower (and 1.12 s for the regular
    insert of the same build with order off).
  - A table function's `Bind` can withhold the index per scan, because the bound copy of the
    function is writable (`input.table_function.get_partition_data = nullptr`);
    `ScanBatchIndex::KeepForSeveralFiles` does so for a one-file scan.
- **A thread that never emits holds the minimum batch index**, so no later batch can be
  flushed past it. A task that waits on background work without emitting stalls an
  order-keeping sink; finish the task (return FINISHED) once it has nothing left to claim.
- Keep batch indexes coarse, or declare `NO_ORDER` when the order means nothing.

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
another reason to preserve input ordering where we already do. (Writing to Parquet keeps
input order by default; miint's sequence writers do not, unless the query
sorts or says `PRESERVE_ORDER true`.)

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

## Exceptions

**`what()` is DuckDB's serialized form of the error.** `duckdb::Exception` derives from
`std::runtime_error`, and its constructor stores the type, the message and any extra fields
(an HTTP error's status, headers and response body) as JSON, `{"exception_type":"IO",
"exception_message":"..."}`, in `what()` (`Exception::ToJSON`,
`duckdb/src/common/exception.cpp`). DuckDB parses it back when it reports the error, so users
see the JSON only when extension code builds a new message from `what()`.
`ErrorData(e).RawMessage()` gives the message alone and passes a `what()` that is not JSON
through unchanged, except `std::bad_alloc`'s, which becomes "Allocation failure"
(`duckdb/src/common/error_data.cpp`); `Message()` adds the type prefix
("IO Error: ...") and, for an internal error, DuckDB's assertion notice and stack trace. A
catch of `std::exception` or `std::runtime_error` can receive a DuckDB exception; a catch of
`std::invalid_argument` or another `std::logic_error` cannot.

## HTTP

### A GET that streams its body (DuckDB 2.0, httpfs)

Given a `response_handler` and a `content_handler`, `GetRequestInfo` hands a GET's body over in
chunks instead of holding it in `response->body`
(`duckdb/src/include/duckdb/main/http/http_util.hpp`). Verified against httpfs `77392a7`
(`src/http/httpfs_curl_client.cpp`, `src/http/httpfs_httplib_client.cpp`), two things are
easy to miss:

- **A retry sends the body again from its first byte.** `HTTPUtil::RunRequestWithRetry`
  ([Interrupts and `max_execution_time`](#interrupts-and-max_execution_time-duckdb-20)) retries a
  GET whose transfer fails partway, since a dropped connection is a request error, and the
  content handler then receives the new attempt's body from the start. The response handler
  runs once per attempt: before a 2xx's body, and after the body of an error status, which
  httpfs keeps in `response->body`. A sink starts over there.
- **A handler must not throw.** The curl client holds an exception from a handler and rethrows
  it after the transfer; the httplib clients do not catch it. Hold it in a
  `std::exception_ptr`, return false, and rethrow once the request returns.

The curl client sets no limit on a transfer's total time, but aborts one that stays below 1 KB/s
for `http_timeout` seconds (`CURLOPT_LOW_SPEED_LIMIT`, `CURLOPT_LOW_SPEED_TIME`), as "Timeout was
reached".

`miint::HttpGetStreaming` (`src/http_client_common.cpp`) handles both, passes only a 2xx's body
on, and checks for an interrupt on every chunk. Use it for any GET whose body is large, such as a
download to a file. Its retries wrap `RunRequestWithRetry`'s
([Interrupts and `max_execution_time`](#interrupts-and-max_execution_time-duckdb-20)).

**To read only the start of a remote file, open it through `FileSystem`, not a GET.** httpfs then
asks for ranges: a plain file's first `Read` of 64 KB fetched 64 KB. A compressed one is read in
pieces of `REMOTE_READ_IN_BUF_SIZE`, 16 MiB, or the whole file if smaller
(`duckdb/src/common/compressed_file_system.cpp`), whatever the reader asks for, and no setting
changes it. `read_csv` fetches more: for the header alone, `LIMIT 1`, `DESCRIBE` and `sniff_csv`
each fetched 32 to 40 MB of a 54 MB gzip file and all of a 29 MB plain one, with `sample_size`,
`buffer_size`, `threads` and the external file cache making no difference (measured 2026-11-10,
counting response bytes with `enable_logging('HTTP')`). Stopping a streamed GET early is no
substitute: a content handler returns false to stop, and the curl client returns the response as
the 200 it was, but the httplib client (deprecated) reports "Canceled", a request error, which
`RunRequestWithRetry` retries up to `http_retries` times. The GTDB readers read headers through
`FileSystem` (`ReadHeader`, `src/read_gtdb.cpp`).
