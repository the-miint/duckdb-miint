#### `read_alignment_header(path)`

<!-- manual:description -->
Read the reference names and lengths from a local SAM or BAM header in header order.

<!-- manual:example -->
```sql
SELECT tid, reference, length FROM read_alignment_header('reads.bam') ORDER BY tid;
```

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `path` | `VARCHAR` | Required | Single local SAM or BAM path. |

| Column | Type | Meaning |
| --- | --- | --- |
| `tid` | `INTEGER` | 0-based index in the `@SQ` list, which is the BAM record's on-disk `refID`. |
| `reference` | `VARCHAR` | Reference name (`SN:`). |
| `length` | `BIGINT` | Reference length (`LN:`), which can exceed 2³¹. |

The function reads only the header. An intact header remains readable if a record is malformed. A header with no `@SQ` lines returns zero rows, as is common for unaligned BAM files written with `COPY … (FORMAT UBAM)`. Use `read_alignments` to stream alignment records from remote storage.

BAM coordinate sort uses `tid` as its reference sort key. `read_alignments` exposes the reference names resolved from the file's own header. To verify the [SAM and BAM writer's reference order](writing.md#sam-and-bam), compare each `tid` to the index produced by sorting the reference names.

```sql
SELECT COUNT(*) = 0 AS header_is_name_sorted FROM (
  SELECT tid, ROW_NUMBER() OVER (ORDER BY reference) - 1 AS expected_tid
  FROM read_alignment_header('sorted.bam')
) WHERE tid != expected_tid;
```

<!-- manual:url /reading/#read_alignment_headerpath -->
