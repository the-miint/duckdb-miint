CIGAR-derived alignment-quality scalars compute per-row metrics from a
read's CIGAR string and (optionally) its NM/MD edit-distance tags. They
are useful as `WHERE` filters or `SELECT` projections downstream of
[`read_alignments`](/reference/table-functions/alignment-io/read_alignments/).

## Identity vs length vs coverage

Three families of metrics, each answering a different question:

- **`alignment_seq_identity(cigar, nm, md, type)`** — *how similar is the
  alignment?* Returns sequence identity in [0, 1]. Four `type` formulas
  trade off correctness vs. tag availability — see the function's
  reference page for the full table; default `'gap_compressed'` matches
  Heng Li's recommendation in
  [his blog post on the definition of sequence identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

- **`alignment_query_length(cigar [, include_hard_clips])`** — *how long
  was the read?* Returns the BIGINT query length implied by the CIGAR.
  With `include_hard_clips=false` (matches HTSlib's `bam_cigar2qlen`)
  the result is the on-disk read length minus hard-clipped bases.

- **`alignment_query_coverage(cigar [, type])`** — *what fraction of the
  read aligned?* Returns DOUBLE in [0, 1]. `type='aligned'` (default)
  counts only `M`/`=`/`X` ops in the numerator; `type='mapped'` also
  counts insertions, so only soft/hard clips reduce coverage.

`alignment_query_length` and `alignment_query_coverage` need only the
CIGAR string. `alignment_seq_identity` needs the NM and/or MD tag values
unless its `type='cigar'` overload is used (`=`/`X` extended CIGAR ops
required).

## Filtering recipe

```sql
-- High-quality, primary, well-aligned reads
SELECT *
FROM read_alignments('sample.bam')
WHERE alignment_is_primary(flags)
  AND mapq >= 30
  AND alignment_query_coverage(cigar) >= 0.9
  AND alignment_seq_identity(cigar, tag_nm, tag_md) >= 0.95;
```
