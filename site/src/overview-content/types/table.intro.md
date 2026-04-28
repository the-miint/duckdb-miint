Table functions return a relation: one or more columns of typed rows that
DuckDB treats as a virtual table. Most miint table functions either **read**
a file format (`read_alignments`, `read_fastx`, `read_biom`, …) or **run a
tool** that produces tabular output (`align_minimap2`, `cluster_sequences`,
…). They are invoked in the `FROM` clause:

```sql
SELECT read_id, mapq
FROM read_alignments('sample.bam')
WHERE mapq >= 30;
```

Many table functions accept **named parameters** that control output schema
or behavior — e.g. `include_seq_qual=true` to include sequence and quality
columns, or `reference_lengths='my_refs'` to supply reference metadata for
headerless SAM. Defaults for named parameters are listed on each function's
reference page.
