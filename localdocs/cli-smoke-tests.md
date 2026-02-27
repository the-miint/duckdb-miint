# Python CLI Smoke Tests

Manual smoke tests for the `miint` Python CLI. Run from the repo root.
Requires the `duckdb-144` conda environment and a built extension.

## Setup

```bash
# Install CLI in editable mode
cd python && conda run -n duckdb-144 pip install -e . && cd ..

# Set extension path
EXT=./build/release/extension/miint/miint.duckdb_extension
```

## Convert commands

```bash
# sequence (single-end)
miint --extension-path $EXT convert sequence \
  -1 data/fastq/small_a.fq -o /tmp/miint_rt.parquet

# sequence (paired-end)
miint --extension-path $EXT convert sequence \
  -1 data/fastq/small_a_r1.fq -2 data/fastq/small_a_r2.fq \
  -o /tmp/miint_paired.parquet

# alignment (SAM)
miint --extension-path $EXT convert alignment \
  -i data/sam/foo_has_header.sam -o /tmp/miint_aln.parquet

# alignment (BAM)
miint --extension-path $EXT convert alignment \
  -i data/sam/foo_with_tags.bam -o /tmp/miint_bam.parquet

# biom
miint --extension-path $EXT convert biom \
  -i data/biom/test.biom -o /tmp/miint_biom.parquet
```

## Convert parquet (round-trip)

```bash
# parquet -> FASTQ
miint --extension-path $EXT convert parquet \
  -i /tmp/miint_rt.parquet -o /tmp/miint_rt.fastq
diff data/fastq/small_a.fq /tmp/miint_rt.fastq  # should match

# parquet -> FASTQ.gz
miint --extension-path $EXT convert parquet \
  -i /tmp/miint_rt.parquet -o /tmp/miint_rt.fastq.gz

# parquet -> FASTA
miint --extension-path $EXT convert parquet \
  -i /tmp/miint_rt.parquet -o /tmp/miint_rt.fasta

# parquet -> SAM (requires --lengths)
echo -e "name\tlength\nchr1\t1000\nchr2\t2000" > /tmp/miint_test_lengths.tsv
miint --extension-path $EXT convert parquet \
  -i /tmp/miint_aln.parquet --lengths /tmp/miint_test_lengths.tsv \
  -o /tmp/miint_rt.sam
```

## Align commands

Requires external test data. Create a 1000-sequence subset first to keep
the test fast (full files take 10+ minutes single-threaded):

```bash
DB=../rype-local-items/perf-data/wol2-genomes/G000157935.fasta.gz
QUERY=../rype-local-items/perf-assessment/query-files/long_read.fastq.gz

# Create subset
conda run -n duckdb-144 python -c "
import duckdb
con = duckdb.connect(config={'allow_unsigned_extensions': 'true'})
con.execute(\"LOAD '$EXT'\")
con.execute(\"COPY (SELECT * FROM read_fastx('$QUERY') LIMIT 1000) TO '/tmp/miint_subset_1k.parquet' (FORMAT PARQUET, parquet_version 'v2', compression 'zstd')\")
"

# align minimap2 (FASTA subject)
miint --extension-path $EXT align minimap2 \
  -1 /tmp/miint_subset_1k.parquet -d $DB \
  -o /tmp/miint_test_align.parquet
```

## Transform commands

Requires alignment parquet and a genome lengths file:

```bash
# Create lengths TSV from the reference FASTA
conda run -n duckdb-144 python -c "
import duckdb
con = duckdb.connect(config={'allow_unsigned_extensions': 'true'})
con.execute(\"LOAD '$EXT'\")
con.execute(\"COPY (SELECT read_id AS name, LENGTH(sequence1) AS length FROM read_fastx('$DB')) TO '/tmp/miint_genome_lengths.tsv' (FORMAT CSV, DELIMITER '\t', HEADER true)\")
"

# genome-coverage
miint --extension-path $EXT transform genome-coverage \
  -i /tmp/miint_test_align.parquet \
  --lengths /tmp/miint_genome_lengths.tsv \
  -o /tmp/miint_test_cov.parquet

# woltka-ogu (no filters)
miint --extension-path $EXT transform woltka-ogu \
  -i /tmp/miint_test_align.parquet \
  -o /tmp/miint_test_ogu.parquet

# woltka-ogu (sequence identity filter)
miint --extension-path $EXT transform woltka-ogu \
  -i /tmp/miint_test_align.parquet \
  --alignment-seq-identity 0.95 \
  -o /tmp/miint_test_ogu_filtered.parquet

# woltka-ogu (all filters)
miint --extension-path $EXT transform woltka-ogu \
  -i /tmp/miint_test_align.parquet \
  --alignment-seq-identity 0.95 \
  --lengths /tmp/miint_genome_lengths.tsv \
  --genome-coverage 0.01 \
  -o /tmp/miint_test_ogu_all.parquet
```

## Error cases

These should all exit non-zero with clean error messages (no tracebacks):

```bash
# SAM output without --lengths
miint --extension-path $EXT convert parquet \
  -i /tmp/miint_aln.parquet -o /tmp/test.sam
# Expected: "error: --lengths is required for SAM output"

# Missing input file
miint --extension-path $EXT convert sequence \
  -1 /tmp/nonexistent.fq -o /tmp/out.parquet
# Expected: "error: --r1 not found: /tmp/nonexistent.fq"

# Bad input extension
miint --extension-path $EXT align minimap2 \
  -1 /tmp/miint_rt.sam -d /tmp/miint_rt.fasta -o /tmp/out.parquet
# Expected: "error: Unrecognized input extension..."

# --genome-coverage without --lengths
miint --extension-path $EXT transform woltka-ogu \
  -i /tmp/miint_test_align.parquet --genome-coverage 0.01 -o /tmp/out.parquet
# Expected: "error: --lengths is required when --genome-coverage is specified"
```

## Not yet tested

- `align minimap2-sharded` (needs RYpe index + shard directory)
- `align minimap2` with `.mmi` index
- `convert parquet` to BIOM or BAM format
