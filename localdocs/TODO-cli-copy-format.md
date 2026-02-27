# TODO: Smarter COPY FORMAT FASTQ/FASTA in CLI

## Current behavior

`miint convert parquet -o output.fastq` uses `_select_for_format()` to hardcode
single-end columns (`read_id, comment, sequence1, qual1` for FASTQ; drop `qual1`
for FASTA). This avoids the `INTERLEAVE parameter required for paired-end data`
error that triggers when `sequence2`/`qual2` columns are present (they always
exist in `read_fastx` output, even when NULL).

## What should be smarter

1. **Paired-end round-trip**: If the parquet actually has non-NULL `sequence2`,
   the current CLI silently drops it. Should detect paired-end data and either:
   - Add `INTERLEAVE true` to produce interleaved output
   - Use `{ORIENTATION}` placeholder to write split R1/R2 files
   - Require an explicit `--interleave` or `--split` flag

2. **Column detection**: Could probe the parquet schema (or sample data) to
   decide single-end vs paired-end automatically:
   ```python
   has_pe = con.execute(
       "SELECT COUNT(*) FROM read_parquet($p) WHERE sequence2 IS NOT NULL LIMIT 1",
       {"p": path}
   ).fetchone()[0] > 0
   ```

3. **`sequence_index` preservation**: `ID_AS_SEQUENCE_INDEX` option exists in
   COPY FASTQ/FASTA but the CLI doesn't expose it. If the parquet has a
   `sequence_index` column and the user wants it as the read ID, we need a flag.

4. **Comment round-trip**: `INCLUDE_COMMENT` defaults to false in COPY FASTQ/FASTA.
   The CLI should probably default to true when the parquet has a non-NULL
   `comment` column.

## Priority

Low — single-end is the common case for now.
