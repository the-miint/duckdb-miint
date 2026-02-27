# TODO: align_minimap2 should handle sequence2 column in subject table

## Problem

`align_minimap2` rejects subject tables that have a `sequence2` column, even
when it's all NULL. Since `read_fastx` always produces `sequence2` in its
schema, users can't do `CREATE TABLE _subjects AS SELECT * FROM read_fastx(...)`
and pass it directly — they must manually select only `read_id, sequence1`.

The CLI works around this with:
```sql
CREATE OR REPLACE VIEW _subjects AS
  SELECT read_id, sequence1 FROM read_fastx('...')
```

## Suggested fix

In `align_minimap2`'s bind function, instead of erroring when `sequence2`
exists, either:
1. Ignore `sequence2` if it's present (just don't read it), or
2. Only error if `sequence2` has non-NULL values (requires a probe query)

Option 1 is simpler and matches user expectations — the function only needs
`read_id` and `sequence1` from subjects.
