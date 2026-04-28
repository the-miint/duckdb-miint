Scalar functions take one or more typed inputs and return a single value
per row. They are the building blocks of `WHERE`, `SELECT`, and aggregate
expressions. miint's scalar functions cluster around two broad concerns:

- **SAM flag bit tests** — `alignment_is_paired`, `alignment_is_unmapped`,
  etc., for filtering alignment rows by their SAM flag bits.
- **CIGAR-derived alignment quality** — `alignment_seq_identity`,
  `alignment_query_length`, `alignment_query_coverage`, for computing
  per-alignment metrics from CIGAR strings (and optionally NM/MD tags).

All scalars compose freely with the rest of SQL — you can use them in
`WHERE`, `GROUP BY`, window functions, CTEs, and subqueries.
