---
title: Per-sample pattern
description: Partitioning a relation by sample_id and running a pipeline per sample.
---

A family of miint table functions accept a `sample_id` named parameter that partitions an input relation by a column and runs the function's pipeline once per distinct sample value. The output gains a prepended sample column whose type matches the source column.

Functions currently using the pattern: `massql`, `woltka_ogu`, `deblur`, `align_mafft`, `detect_chimera_uchime_denovo`, `detect_chimera_uchime`.

## The helper

`src/include/per_sample_table_function.hpp` exposes a small, opinionated API. Callers use it for bind-time sample discovery and exec-time atomic claim; everything else — the non-sample branch, data loading, result buffering, chunk emission — stays in the caller.

```cpp
struct PerSampleBindInfo {
    string sample_id_col;
    LogicalType sample_id_type;
    vector<Value> sample_values;
};

void DiscoverSamples(Connection &conn, const string &source, const string &sample_id_col,
                     const std::vector<std::string> &reserved_lowercase_output_names,
                     const string &fn_label, PerSampleBindInfo &out);

struct PerSampleGlobalState : GlobalTableFunctionState {
    atomic<idx_t> next_sample_idx{0};
    idx_t max_threads = 1;
    idx_t MaxThreads() const override { return max_threads; }
};

void InitPerSampleGlobal(ClientContext &context, PerSampleGlobalState &gstate,
                         idx_t num_samples, idx_t max_threads_hint = 0);

bool ClaimNextSample(PerSampleGlobalState &gstate, idx_t num_samples, idx_t &out_idx);
```

### Error strings

All helper-thrown errors are prefixed with `<fn_label>: `. The substrings that downstream tests have historically asserted on are preserved verbatim (`sample_id column '<X>' not found`, `sample_id column name must not be empty`, `NULL values in sample_id column '<X>'`, `sample_id column '<X>' has no non-NULL values`, `collides with an output column`). When adding a new per-sample function, reuse the same substrings in tests — they will all come out of the helper.

### Thread-safety knob

`max_threads_hint` controls the concurrency cap:

- `0` (default): use the DuckDB scheduler's thread count, clamped to `num_samples`. Use when the per-sample backend is re-entrant.
- `1`: force serial execution. Use when the backend has a process-wide mutex (MAFFT) or is explicitly not thread-safe (vsearch's `VsearchChimeraWrapper`). Spawning more threads would just queue them behind the lock — worse, not better.
- `N > 1`: cap to `min(N, num_samples, db_threads)`. Use if a backend has limited internal parallelism.

## Adding per-sample support to a new function

Three steps.

### 1. Named parameter + bind-time discovery

```cpp
tf.named_parameters["sample_id"] = LogicalType::VARCHAR;
```

```cpp
auto sample_it = input.named_parameters.find("sample_id");
if (sample_it != input.named_parameters.end()) {
    data->has_sample_id = true;
    data->sample_info.sample_id_col = sample_it->second.GetValue<string>();
}

if (data->has_sample_id) {
    auto &db = DatabaseInstance::GetDatabase(context);
    Connection conn(db);
    DiscoverSamples(conn, data->input_table, data->sample_info.sample_id_col,
                    /*reserved=*/ {/* lowercase output column names */},
                    "<fn_name>", data->sample_info);

    names.push_back(data->sample_info.sample_id_col);
    return_types.push_back(data->sample_info.sample_id_type);
}
```

Store `PerSampleBindInfo sample_info` in the Data struct — not its three fields individually; that's the historical anti-pattern this helper was carved out of.

### 2. Global / local state

Inherit the GlobalState from `PerSampleGlobalState`; initialize it at `InitGlobal`:

```cpp
if (data.has_sample_id) {
    InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size(),
                        /*max_threads_hint=*/0 /* or 1 for non-re-entrant backends */);
} else {
    gstate->max_threads = 1;
    // ... existing non-sample setup ...
}
```

Provide a `LocalState` with a `unique_ptr<Connection> conn` (for the per-thread data load) and whatever per-sample result buffer the function needs. The Connection is only allocated in the sample_id path.

### 3. Per-sample Execute loop

Template used by all six existing callers:

```cpp
while (true) {
    if (lstate.current_row < lstate.results.size()) {
        EmitRows(/* ... */);  // emit a chunk from the per-sample buffer
        return;
    }
    idx_t sample_idx;
    if (!ClaimNextSample(gstate, data.sample_info.sample_values.size(), sample_idx)) {
        output.SetCardinality(0);
        return;
    }
    lstate.sample_value = data.sample_info.sample_values[sample_idx];
    // Load this sample's rows and run the algorithm. Use the per-thread lstate.conn.
    lstate.results = RunAlgorithmForSample(*lstate.conn, data, lstate.sample_value);
    lstate.current_row = 0;
}
```

When emitting, the sample column at output position 0 should be filled via a ConstantVector reference — constant across the chunk, zero-copy:

```cpp
output.data[0].Reference(lstate.sample_value);
```

Avoid the per-row `SetValue` anti-pattern; a 2048-row chunk would otherwise do 2048 string-interning operations for the sample column alone.

## Per-sample data loading

Two loader patterns are in use:

**SQL transpile** (`massql`, `woltka_ogu`). The function's per-sample SQL is wrapped by `CREATE OR REPLACE TEMP VIEW __<fn>_per_sample AS SELECT * FROM <src> WHERE <col> = <literal>` on the thread's Connection, then the existing SQL runs against the view. The TEMP VIEW is scoped to the connection, so parallel threads don't collide despite the shared name.

**In-memory load** (`deblur`, `align_mafft`, `uchime_denovo`, `uchime_ref`). The function pulls rows directly:
- `deblur` and `uchime_denovo` run custom `SELECT … WHERE <sample_predicate> ORDER BY <ordering>` queries.
- `align_mafft` uses `LoadSingleEndSequences(Connection&, table, fn, strict, where_sql)` — an overload that accepts a caller-owned Connection and an arbitrary WHERE clause. Existing callers that don't want a predicate keep the original `ClientContext&` overload.
- `uchime_ref` streams queries lazily via `QuerySequenceStream(Connection&, table, schema, batch_size)` — an overload that uses the caller's Connection so the TEMP VIEW and the stream share a session (the view isn't visible otherwise).

All per-sample predicates use `CAST(<col> AS VARCHAR) = CAST(<literal> AS VARCHAR)`. This is correct for VARCHAR, integer, and date sample IDs. **DECIMAL sample IDs are not supported**: DuckDB's VARCHAR cast preserves trailing zeros (`3.500`) while `Value::ToSQLString()` strips them (`3.5`), so equality silently misses rows. If DECIMAL sample IDs become a requirement, switch every call site to type-aware literal binding.

## Reserved output-column names

Each caller passes a list of lowercase output column names that the `sample_id` column must not collide with (case-insensitive). This prevents the output from having two columns with the same name.

Current lists:
- `deblur`: `{read_id, sequence, abundance}`
- `align_mafft`: `{sequence_index, read_id, aligned_sequence, original_length, aligned_length}`
- `uchime_denovo` and `uchime_ref`: `GetUchimeOutputNames()` (the 18 uchimeout columns).
- `woltka_ogu`: `{feature_id, value}`
- `massql`: `{}` (no fixed output schema — the schema is inferred from sample[0]).

Input-column collisions (e.g. `sample_id := 'sequence1'`) are intentionally *not* rejected. The user picks their own foot, but there's no ambiguity in the output.

## Testing

Each new per-sample function gets a `test/sql/<fn>_sample_id.test` covering:
1. Missing column → `sample_id column '<name>' not found`
2. Empty string → `sample_id column name must not be empty`
3. NULL values → `NULL values in sample_id column '<name>'`
4. Empty relation filter → `has no non-NULL values`
5. Output-column collision (case-insensitive) → `collides with an output column`
6. Happy path with two samples — verify isolation and output schema
7. Integer-typed sample_id round-trips
8. Single-sample equivalence vs. the non-sample path
9. `PRAGMA threads=4` parallelism cycle — exercises the atomic claim

`massql_sample_id.test` and `deblur_sample_id.test` are canonical references for layout.
