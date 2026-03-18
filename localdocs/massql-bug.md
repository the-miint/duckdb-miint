# Bug: Python MassQL discards qualifier operator (`<`, `=`, `>` all produce identical behavior)

## Summary

The MassQL parser (`msql_parser.py`) discards the comparison operator on qualifier expressions. All three of `INTENSITYPERCENT=50`, `INTENSITYPERCENT<50`, and `INTENSITYPERCENT>50` produce identical results — the value is always treated as a minimum threshold (`i_norm > value/100`).

## Affected version

`massql==2025.12.10` (latest as of 2026-03-12)

## Minimal reproduction

```python
from massql import msql_engine

# Create a tiny mzML with one MS1 spectrum containing two peaks:
#   mz=100 intensity=1000 (i_norm=0.2, i.e. 20%)
#   mz=200 intensity=5000 (i_norm=1.0, i.e. 100%)
mzml = "data/mzml/basic_3spectra.mzML"  # or any mzML with peaks at varying intensities

# Peak at mz=100 has i_norm=0.2 (20%).
# INTENSITYPERCENT<50 should MATCH (20% < 50%), but returns 0 rows:
df_lt = msql_engine.process_query(
    "QUERY scannum(MS1DATA) WHERE MS1MZ=100:INTENSITYPERCENT<50", mzml
)
print(f"<50 on 20% peak: {len(df_lt)} rows (expected 1)")  # prints 0

# Peak at mz=200 has i_norm=1.0 (100%).
# INTENSITYPERCENT<50 should NOT match (100% is not < 50%), but returns 1 row:
df_lt2 = msql_engine.process_query(
    "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT<50", mzml
)
print(f"<50 on 100% peak: {len(df_lt2)} rows (expected 0)")  # prints 1

# For comparison, >50 on the same peaks:
df_gt = msql_engine.process_query(
    "QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT>50", mzml
)
print(f">50 on 100% peak: {len(df_gt)} rows (expected 1)")  # prints 1 (correct)

df_gt2 = msql_engine.process_query(
    "QUERY scannum(MS1DATA) WHERE MS1MZ=100:INTENSITYPERCENT>50", mzml
)
print(f">50 on 20% peak: {len(df_gt2)} rows (expected 0)")  # prints 0 (correct)
```

**Expected output:**
```
<50 on 20% peak: 1 rows (expected 1)
<50 on 100% peak: 0 rows (expected 0)
>50 on 100% peak: 1 rows (expected 1)
>50 on 20% peak: 0 rows (expected 0)
```

**Actual output:**
```
<50 on 20% peak: 0 rows (expected 1)
<50 on 100% peak: 1 rows (expected 0)
>50 on 100% peak: 1 rows (expected 1)
>50 on 20% peak: 0 rows (expected 0)
```

The `<` results are inverted — they behave identically to `>`.

## Root cause

**`msql_parser.py` lines 141-158**, `MassQLToJSON.qualifier()`:

The method receives three items for a qualifier like `INTENSITYPERCENT<50`:
- `items[0]` = `"qualifierintensitypercent"` (field name)
- `items[1]` = `Tree('lessthan', [])` (the operator)
- `items[2]` = `50.0` (the value)

It stores `items[0]` and `items[-1]` but **skips `items[1]`** entirely. The operator is never written to the parsed JSON. All three operators produce the same parsed output:

```json
{"qualifierintensitypercent": {"name": "qualifierintensitypercent", "value": 50.0}}
```

**`msql_engine_filters.py` lines 29-63**, `_get_minintensity()`:

Takes the qualifier value and returns it as `min_percent_intensity = value / 100`. This is then used in a hardcoded `i_norm > min_intpercent` comparison (e.g. line 468). There is no code path for `<` or `=`.

## Impact

- The `<` operator is silently treated as `>`, producing inverted results.
- The `=` operator is silently treated as `>`, producing different results (should be `>=`).
- Only `>` works correctly.
- The EBNF grammar (`msql.ebnf` lines 35-46) correctly defines `lessthan` and `greaterthan` as distinct tokens; the grammar is not the problem.
- Affects `INTENSITYPERCENT`, `INTENSITYVALUE`, and `INTENSITYTICPERCENT` — any qualifier where operator direction matters.
- No existing test in the MassQL test suite exercises `<` or `=` operators, so the bug has not been caught.
