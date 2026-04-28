Alignment-IO table functions read SAM and BAM files (text and binary
respectively, including gzip-compressed SAM and remote `https://` paths)
into DuckDB-native columnar relations that downstream SQL can filter,
join, and aggregate.

The canonical entry point is
[`read_alignments`](./read_alignments/), which:

- Accepts a single path, an array of paths, a glob pattern, or stdin.
- Auto-detects SAM vs BAM from the file content (no extension required).
- Optionally adds sequence and quality columns (`include_seq_qual=true`).
- Optionally adds a per-row source filename column (`include_filepath=true`).
- Supports headerless SAM via a `reference_lengths` named parameter that
  points at a DuckDB table or view of `(name, length)` rows.

`read_sam` is a backward-compatibility alias of `read_alignments` —
predates the rename and is preserved so existing pipelines don't break.
The two are identical in every other respect; new code should prefer
`read_alignments`.

## Composing with the rest of miint

Alignment-IO output flows naturally into the
[SAM flag scalars](../sam-flags/) and the
[CIGAR-derived alignment-quality scalars](../alignment-quality/):

```sql
SELECT read_id, mapq, cigar
FROM read_alignments('sample.bam', include_seq_qual=false)
WHERE alignment_is_primary(flags)
  AND mapq >= 30
  AND alignment_query_coverage(cigar) >= 0.9;
```

For aligning sequences (rather than reading existing alignments) see the
`align_*` table functions (`align_minimap2`, `align_bowtie2`, `align_mafft`,
`align_sortmerna`) — those will land in their own categories as the
`RegisterDocumented*` rollout reaches them.
