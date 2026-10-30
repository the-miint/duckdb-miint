# Multi-part minimap2 index fixture

`shard_a.mmi` and `shard_b.mmi` are described in `test/sql/align_minimap2_sharded_single_read.test`. This file covers the multi-part fixture:

- `multipart_fixture.mmi` is a two-part minimap2 index: `part1_ref` in its first part, `part2_ref` in its second.
- `multipart_fixture_part1.mmi` and `multipart_fixture_part2.mmi` are its two parts, each a single-part index on its own, so that a multi-part run can be compared with two single-index runs.

`save_minimap2_index` always builds one part. The minimap2 CLI's `-I <batch>` builds several by writing one `mm_idx_dump` block per part into the same file, and each single-part file is already a complete block, so concatenating two gives the same layout without the CLI:

```sql
CREATE TABLE p1 AS SELECT * FROM (VALUES
    ('part1_ref', 'ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC')
) AS t(read_id, sequence1);
CREATE TABLE p2 AS SELECT * FROM (VALUES
    ('part2_ref', 'TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA')
) AS t(read_id, sequence1);
SELECT success FROM save_minimap2_index('p1', 'data/shards/multipart_fixture_part1.mmi', k := 5);
SELECT success FROM save_minimap2_index('p2', 'data/shards/multipart_fixture_part2.mmi', k := 5);
```

```bash
cat data/shards/multipart_fixture_part1.mmi data/shards/multipart_fixture_part2.mmi > data/shards/multipart_fixture.mmi
```

The queries in `align_minimap2_multipart.test` that are meant to match, such as `query_part1` and `query_part2`, are prefixes of these two references. If you rebuild the fixture from other references, change those queries with it.
