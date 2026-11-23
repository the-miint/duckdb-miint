### PAF

<!-- manual:description -->
Read Pairwise Mapping Format (PAF) alignment records with 1-based half-open query and target intervals.

<!-- manual:example -->
```sql
SELECT * FROM 'alignments.paf';
SELECT * FROM 'alignments.paf.gz';
```

DuckDB uses `read_paf` for these [file suffixes](#reading-a-file-by-its-name). Call the function to pass options or open a file with a different suffix.

#### Read several PAF files

<!-- manual:example -->
```sql
SELECT * FROM read_paf('samples/*.paf.gz', include_filepath = true);
```

#### PAF parameters

| Parameter | Type | Default | Purpose |
| --- | --- | --- | --- |
| `filename` | `VARCHAR` or `VARCHAR[]` | Required | Path, glob, or list of literal paths. |
| `include_filepath` | `BOOLEAN` | `false` | Include the source path. |

#### PAF output columns

Each row is one PAF record.

| Column | Type | Meaning |
| --- | --- | --- |
| `read_id` | `VARCHAR` | Query sequence name. |
| `query_length` | `BIGINT` | Query sequence length. |
| `query_start` | `BIGINT` | 1-based start of the aligned query block. |
| `query_stop` | `BIGINT` | Exclusive end of the aligned query block. |
| `strand` | `VARCHAR` | `+` or `-`, or `*` for an unmapped record. |
| `reference` | `VARCHAR` | Target sequence name, or `*` for an unmapped record. |
| `reference_length` | `BIGINT` | Target sequence length. |
| `position` | `BIGINT` | 1-based start of the aligned target block. |
| `stop_position` | `BIGINT` | Exclusive end of the aligned target block. |
| `residue_matches` | `BIGINT` | Number of matching bases. |
| `block_length` | `BIGINT` | Alignment block length. |
| `mapq` | `UTINYINT` | Mapping quality, with 255 indicating a missing value. |
| `tags` | `MAP(VARCHAR, VARCHAR)` | Optional fields without a dedicated column, keyed by tag name. |
| `filepath` | `VARCHAR` | Source path, returned with `include_filepath = true`. |

#### PAF coordinates and unmapped records

Both intervals use MIINT's [1-based half-open convention](#coordinate-conventions). The reader adds 1 to both ends of each mapped interval. The target interval matches `read_alignments`, so `stop_position - position` is the aligned target length, and `query_stop - query_start` is the aligned query length. The query interval is on the [read as sequenced](#reads), including on the reverse strand.

An unmapped record, such as one produced by `minimap2 --paf-no-hit`, has `*` for `strand` and `reference` and 0 for all four interval endpoints.

#### Optional PAF tags

Minimap2's PAF tags have typed columns. An absent tag is `NULL`. Minimap2 writes the `SA` tag only in SAM.

| Column | Type | Tag |
|---|---|---|
| `tag_tp` | VARCHAR | `tp:A`, alignment type (`P` primary, `S` secondary, `I` / `i` inversion) |
| `tag_cm` | BIGINT | `cm:i`, minimizers on the chain |
| `tag_s1` | BIGINT | `s1:i`, chaining score |
| `tag_s2` | BIGINT | `s2:i`, chaining score of the best secondary chain |
| `tag_nm` | BIGINT | `NM:i`, mismatches and gaps in the alignment |
| `tag_md` | VARCHAR | `MD:Z`, for regenerating the reference sequence of the alignment |
| `tag_as` | BIGINT | `AS:i`, DP alignment score |
| `tag_ms` | BIGINT | `ms:i`, DP score of the best-scoring segment |
| `tag_nn` | BIGINT | `nn:i`, ambiguous bases in the alignment |
| `tag_ts` | VARCHAR | `ts:A`, transcript strand (splice mode) |
| `tag_cg` | VARCHAR | `cg:Z`, CIGAR |
| `tag_cs` | VARCHAR | `cs:Z`, difference string |
| `tag_dv` | DOUBLE | `dv:f`, approximate per-base divergence |
| `tag_de` | DOUBLE | `de:f`, gap-compressed per-base divergence |
| `tag_rl` | BIGINT | `rl:i`, length of query regions with repetitive seeds |
| `tag_zd` | BIGINT | `zd:i`, Z-drop break (bit 1 left, bit 2 right) |

The reader converts `NM:i:3` to the integer `tag_nm = 3`. For tags without a dedicated column, it validates the type letter against the value and stores only the text value. For example, `zz:i:7` becomes `tags['zz'] = '7'`. Cast map values when you need numbers.

Float tags (`f`) accept `nan`, `inf`, and `infinity` in any case with an optional sign. Values outside the finite `DOUBLE` range overflow to infinity. The reader also accepts untyped `key:value` fields. When a key has a dedicated column, the reader converts the value to that column's type.

#### PAF paths, compression, and parallel reads

A glob's matches are sorted alphabetically. Lists contain literal paths. The reader decompresses gzip files, including multi-member (bgzip) files, and opens remote paths such as `s3://` and `https://` through DuckDB's file system. Files can be read in parallel, with one thread per file. Empty files return zero rows. PAF input must be a file path, glob, or list of paths.

#### Malformed PAF records

A malformed record fails the scan. The error identifies the file, line, and defect.

- fewer than 12 fields, or a blank line
- an empty query or target name
- a fractional, negative, or non-numeric value in one of the first 12 fields that requires an integer, or a mapping quality above 255
- an interval that runs backwards or past its sequence's length
- a strand other than `+` or `-`, except `*` on an unmapped record
- `*` as only one of strand and target name, or an unmapped record with a non-zero query or target start or end
- a query or target end too large to take the `+ 1`
- a malformed tag, or a tag that appears twice on one line
- one of minimap2's tags with a type letter other than the one in the table, or a value that does not fit its column (an integer beyond 64 bits)

#### Filter and summarize PAF alignments

```sql
-- Primary alignments of reads with at most 2% divergence
SELECT read_id, reference, position, stop_position
FROM read_paf('reads.paf.gz')
WHERE tag_tp = 'P' AND tag_de <= 0.02;

-- Fraction of each read covered by its alignments, pooling fragments on the read
SELECT read_id, query_length,
       list_sum([iv.stop - iv.start FOR iv IN compress_intervals(query_start, query_stop)]) / query_length
           AS fraction_aligned
FROM read_paf('reads.paf')
WHERE reference <> '*'
GROUP BY read_id, query_length;

-- Unaligned reads
SELECT read_id FROM read_paf('reads.paf') WHERE reference = '*';
```

<!-- manual:url /reading/#paf -->
