# MIINT Documentation

MIINT provides a range of bioinformatic capabilities. The materials included here are coupled where possible with detailed explanation and examples.

## Passing relations by name

Many MIINT functions take a table or view **by name**, as a quoted string literal rather than a relation
argument — `align_minimap2('reads', subject_table := 'refs')`. Those names are resolved against the
catalog, so a `CREATE TEMP TABLE` or `CREATE TEMP VIEW` works anywhere an ordinary table does. This
matters most under a read-only catalog, where a temporary relation may be the only one you can create.

Two exceptions, both with a note in the relevant function's documentation:

- [`massql()`](massql.md) and the per-sample (`sample_id := ...`) paths of
  [`woltka_ogu`](profiling.md), [`sylph_profile`](profiling.md) and
  [`detect_chimera_uchime`](chimera.md) cannot read a temporary relation. Materialize the source as an
  ordinary table first.
- A relation created inside an explicit `BEGIN` that has not been committed yet is not visible to any of
  these functions, temporary or not.

A subquery cannot be used in place of the literal — DuckDB rejects it outright with `Binder Error: Table
function cannot contain subqueries`. To choose the relation at runtime, use `SET VARIABLE` plus
`getvariable()`:

```sql
SET VARIABLE reads = 'my_temp_reads';
SELECT * FROM align_minimap2(getvariable('reads'), subject_table := 'refs');
```

### Arrow relations

A relation registered from Arrow works by name like any other, because DuckDB's registration creates a
temporary view. This includes a `RecordBatchReader` streamed from an external source such as Arrow
Flight, so a streaming query set can be handed to a function directly rather than spilled to Parquet
first:

```python
con.register('reads', flight_reader)          # pyarrow Table or RecordBatchReader
con.execute("SELECT * FROM align_minimap2('reads', subject_table := 'refs')")
```

Two constraints come with it:

- **A `RecordBatchReader` is consumed once.** A second reference to the same registered name returns zero
  rows, silently — that is Arrow's semantics, not something MIINT can detect. `align_minimap2` (both
  modes) and `align_minimap2_sharded` read their relation exactly once and are safe;
  `align_bowtie2_sharded` re-reads per shard and will raise an error rather than return a partial result.
  Other relation-taking functions have **not** been audited for this yet — notably the per-sample
  (`sample_id := ...`) paths, which read the relation once per sample — so prefer a materialized table
  there. See `docs/internals/reading-tables-views.md` for the audit status.
- **Do not register a reader that is itself backed by a DuckDB query on the same connection.** In DuckDB
  1.5.4 `con.execute(...).arrow()` returns exactly that — a lazy `RecordBatchReader`, not a materialized
  table — and passing it to a MIINT function **deadlocks the process** with no error and no timeout. The
  cause is a lock-ordering hazard inside DuckDB's own Arrow scan, and nothing in MIINT can intercept it.
  Materialize first: `con.execute(...).arrow().read_all()`, or `CREATE TEMP TABLE ... AS SELECT`.

## Reading & writing

- [Reading files](reading.md) - File formats with read support in MIINT (FASTA/FASTQ, SAM/BAM, SFF, BIOM, mzML/mzXML, GFF, jplace, Newick).
- [Writing files](writing.md) - File formats with write support in MIINT.

## Public data repositories

- [EBI/ENA](insdc_ena.md) - Interaction with EBI/ENA including submission and fetching of data.
- [NCBI](insdc_ncbi.md) - Interaction with NCBI including fetching of data and the NCBI taxonomy (tree, name classes, lineages, retired/deleted-taxid handling).

## Alignment

- [Reference based alignment](alignment_reference.md) - High throughput reference based read alignment.
- [Pairwise alignment](alignment_pairwise.md) - Pairwise alignment methods.
- [Multiple sequence alignment](alignment_multiple.md) - Multiple sequence alignment methods.
- [Alignment analysis](alignment_analysis.md) - Operations for the analysis of alignment data.

## Amplicon & community analysis

- [Quality control](qc.md) - Adapter / polyG / polyX / quality trimming and per-read filtering.
- [Denoising](denoising.md) - Denoise sequencing data into sub-OTU variants.
- [Chimera checking](chimera.md) - Reference-based and de novo chimera detection.
- [Sequence clustering](clustering.md) - Cluster sequences into OTUs at an identity threshold.
- [Sequence search](search.md) - Search sequences against a reference database.
- [Taxonomic classification](classification.md) - Classify sequences by minimizer/k-mer content (RYpe).
- [Profiling & feature tables](profiling.md) - Abundance profiling (sylph) and OGU feature tables (Woltka).
- [Absolute quantification](absolute_quantification.md) - Turn compositional read counts into absolute quantities: cell counts from synDNA spike-ins (Zaramela et al. 2022), and copies of each ORF's ssRNA.

## Phylogeny & diversity

- [Phylogeny estimation](phylogeny.md) - Phylogenetic estimation, tree manipulation (shear, resolve multifurcations/placements), and comparative methods (independent contrasts, ancestral state reconstruction: Brownian-motion, parsimony, and Mk maximum likelihood).
- [Alpha and beta diversity](diversity.md) - Methods to compute and analyze alpha and beta diversity, including non-phylogenetic community distances, ordination, and sample clustering.
- [Community simulation](simulation.md) - Generate synthetic gradient / clustered OTU tables with known ground truth for benchmarking resemblance and ordination methods.
- [Multi-omics integration](multiomics.md) - MMvec: learn which metabolites co-occur with which microbes from paired count tables, without the spurious associations correlation produces on compositional data.

## Mass spectrometry

- [MassQL](massql.md) - Operating on mass spectrometry data through MassQL.

## Other

- [Utilities](utilities.md) - Sequence helpers, versions/diagnostics, and optional-tool installation.
