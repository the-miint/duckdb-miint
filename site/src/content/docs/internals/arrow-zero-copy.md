---
title: Arrow zero-copy
description: Consuming Arrow C Data Interface streams without copying.
---

When consuming Arrow C Data Interface streams (e.g., from RYpe or other FFI libraries), use DuckDB's built-in Arrow-to-DuckDB conversion instead of manually copying data. This gets zero-copy for fixed-width types (integers, floats) including inside List children.

## Key Components

- `ArrowTableFunction::PopulateArrowTableSchema()` — parses Arrow schema into `ArrowTableSchema` with type metadata needed by the converter
- `ArrowToDuckDBConversion::SetValidityMask()` — copies the null bitmap from Arrow to DuckDB
- `ArrowToDuckDBConversion::ColumnArrowToDuckDB()` — converts one column, dispatching by type:
  - Fixed-width types (UBIGINT, DOUBLE, etc.) → `DirectConversion()` → `FlatVector::SetData()` (zero-copy pointer swap)
  - List types → converts offset buffer, then recursively converts child (child gets zero-copy)
  - Strings → copies (different layout between Arrow and DuckDB)

## Lifetime Management for Zero-Copy

Arrow buffers must stay alive while DuckDB Vectors reference them. The pattern:

1. Wrap each Arrow batch in `shared_ptr<ArrowArrayWrapper>` (not raw `ArrowArray`)
2. Set `ArrowArrayScanState::owned_data = shared_ptr` before calling `ColumnArrowToDuckDB`
3. DuckDB stores the `shared_ptr` in `ArrowAuxiliaryData` on the Vector's buffer, preventing premature deallocation

```cpp
// In GlobalState: use shared_ptr, NOT raw ArrowArray
shared_ptr<ArrowArrayWrapper> current_chunk;

// Fetching a new batch:
auto wrapper = make_shared_ptr<ArrowArrayWrapper>();
stream.get_next(&stream, &wrapper->arrow_array);
current_chunk = std::move(wrapper);

// In Execute, for each Arrow column that maps directly to a DuckDB column:
auto &array_state = lstate.GetState(col_idx);
array_state.owned_data = gstate.current_chunk;  // ref-count keeps batch alive
ArrowToDuckDBConversion::SetValidityMask(output.data[col], array, batch_offset, size, batch.offset, -1);
ArrowToDuckDBConversion::ColumnArrowToDuckDB(output.data[col], array, batch_offset, array_state, size, arrow_type);
```

## When Manual Access is Needed

For columns that require transformation (e.g., id→read_id lookup), direct buffer access is fine for *that column* — use `reinterpret_cast` on `array.buffers[1]` as in `src/rype_classify.cpp`.

## Anti-pattern

Do **NOT** manually iterate Arrow List child data with element-by-element copy. Use `ColumnArrowToDuckDB`, which handles List→LIST conversion with zero-copy on the child values.

## Headers

```cpp
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/function/table/arrow.hpp"
```

## Reference Implementations

- `src/rype_extract.cpp` — zero-copy Arrow List<UInt64> → DuckDB LIST(UBIGINT)
- `src/rype_classify.cpp` — manual Arrow access (all columns need transformation)
- `duckdb/src/function/table/arrow_conversion.cpp` — DuckDB's internal conversion engine
- `duckdb/src/main/capi/arrow-c.cpp` (`duckdb_data_chunk_from_arrow`) — C API usage example
