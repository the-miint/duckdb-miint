SAM flag scalars test individual bits of the 16-bit `flags` column produced
by [`read_alignments`](/reference/table-functions/alignment-io/read_alignments/).
Bit definitions follow the [SAM specification](https://samtools.github.io/hts-specs/SAMv1.pdf)
(section 1.4):

| Bit    | Canonical                     | HTSlib alias        |
|-------:|-------------------------------|---------------------|
| `0x1`  | `alignment_is_paired`         | `is_paired`         |
| `0x2`  | `alignment_is_proper_pair`    | `is_proper_pair`    |
| `0x4`  | `alignment_is_unmapped`       | `is_unmapped`       |
| `0x8`  | `alignment_is_mate_unmapped`  | `is_munmap`         |
| `0x10` | `alignment_is_reverse`        | `is_reverse`        |
| `0x20` | `alignment_is_mate_reverse`   | `is_mreverse`       |
| `0x40` | `alignment_is_read1`          | `is_read1`          |
| `0x80` | `alignment_is_read2`          | `is_read2`          |
| `0x100`| `alignment_is_secondary`      | `is_secondary`      |
| `0x200`| `alignment_is_qc_failed`      | `is_qcfail`         |
| `0x400`| `alignment_is_duplicate`      | `is_dup`            |
| `0x800`| `alignment_is_supplementary`  | `is_supplementary`  |

`alignment_is_primary` is the only entry without a SAM bit of its own —
it returns true when **both** `alignment_is_secondary` (0x100) and
`alignment_is_supplementary` (0x800) are false. Most analyses want
`alignment_is_primary(flags)` rather than negating two other tests.

The HTSlib aliases mirror the predicate names available on
[`samtools view --filter`](http://www.htslib.org/doc/samtools-view.html);
they exist so queries written against `samtools` cheat sheets translate
unchanged.

## Composing flag tests

Flag predicates take a `USMALLINT` and return `BOOLEAN`, so they compose
with normal SQL booleans:

```sql
-- Properly-paired primary alignments only
SELECT *
FROM read_alignments('sample.bam')
WHERE alignment_is_proper_pair(flags)
  AND alignment_is_primary(flags)
  AND NOT alignment_is_qc_failed(flags);
```
