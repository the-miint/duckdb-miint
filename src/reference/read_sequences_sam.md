#### `read_sequences_sam(filename)`

<!-- manual:description -->
Read SAM and BAM sequences in the same column layout as `read_fastx`.

<!-- manual:example -->
```sql
SELECT read_id, sequence1, qual1 FROM read_sequences_sam('reads.bam');
```

<!-- manual:example -->
```sql
SELECT * FROM read_sequences_sam('samples/*.bam', include_filepath = true);
```

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `filename` | `VARCHAR` or `VARCHAR[]` | Required | Path, glob, list of literal paths, or `-` / `/dev/stdin` as the only input. |
| `include_filepath` | `BOOLEAN` | `false` | Include the source path. |

| Column | Type | Meaning |
| --- | --- | --- |
| `sequence_index` | `BIGINT` | 1-based record index within each file. |
| `read_id` | `VARCHAR` | Query name (QNAME). |
| `comment`, `sequence2` | `VARCHAR` | Always `NULL`. |
| `sequence1` | `VARCHAR` | Sequence as stored in the alignment record. |
| `qual1` | `UTINYINT[]` | Phred scores, or `NULL` if absent. |
| `qual2` | `UTINYINT[]` | Always `NULL`. |
| `filepath` | `VARCHAR` | Source path, returned with `include_filepath = true`. |

Each alignment record produces one row. Mates remain separate records even when they share a query name. Primary and unmapped records must carry SEQ. Use `read_alignments` for flags, coordinates, or alignment tags.

<!-- manual:url /reading/#read_sequences_samfilename -->
