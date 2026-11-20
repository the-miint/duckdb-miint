## Writing Prodigal's FASTA

```sql
prodigal_proteins(relation, [per_genome := false])
prodigal_genes(relation, contigs, [genome_id := 'col'], [per_genome := false])
```

<!-- manual:url https://the-miint.github.io/duckdb-miint/gene_prediction/#writing-prodigals-fasta -->
<!-- manual:description -->
Convert `predict_genes` rows into Prodigal protein FASTA records. `per_genome` defaults to false.

The table macros make Prodigal's `-a` (proteins) and `-d` (gene nucleotides) FASTA records from `predict_genes`' rows for [`COPY ... (FORMAT FASTA)`](writing.md#fasta) to write. On the contigs in `data/prodigal/`, single-genome and metagenomic, the files hold the Prodigal CLI's records byte for byte. The records come out in no set order, as a FASTA file's records are a set; the CLI writes them in contig order.

- `relation` (VARCHAR, positional): a table or view of `predict_genes`' output, in its order (contig order).
- `contigs` (VARCHAR, `prodigal_genes` only): a table or view of the contigs `predict_genes` was given, `read_id` and `sequence1`. A gene's bases come from the contig with its `read_id`, which must be there once, with `read_id` of the genes' type. A gene whose contig is missing, named twice, NULL, too short for it, or holds a character outside ASCII `A`..`z` (which `predict_genes` refuses), and a gene that is empty (`stop_position` not after `position`), are refused.
- `genome_id` (VARCHAR, `prodigal_genes` only): the `contigs` column of genomes, as for `predict_genes`. Given, a gene's contig is the one with its genome and `read_id`, so genomes may share contig names (`contig_1`, ...).
- `per_genome` (BOOLEAN, default false): contigs are numbered (the `n` in `ID=n_g`) from 1 across the relation, as [`COPY ... (FORMAT GFF, DIALECT 'prodigal')`](writing.md#gff) numbers them in one file; true numbers each genome's from 1, as one CLI run per genome does.

They return `genome_id`, `read_id` (`<contig>_<gene>`), `comment` (`# begin # end # 1|-1 # ID=...;partial=..;start_type=..;rbs_motif=..;rbs_spacer=..;gc_cont=..`) and `sequence1`. Prodigal wraps `-a` at 60 characters and `-d` at 70, so write them with `INCLUDE_COMMENT true` and `LINE_WIDTH 60` or `70`. `-d`'s bases are the contig's as Prodigal reads them: `a`/`c`/`g`/`t` in either case are `ACGT`, every other character (other letters, and ``[\]^_` ``) is `N`, and a `-` gene is the reverse complement. A contig's row without a gene is numbered but has no record; a genome's own row is skipped.

Given a `genes` table or view from `predict_genes`, preview the protein records.

<!-- manual:example -->
```sql
SELECT * FROM prodigal_proteins('genes');
```

```sql
CREATE TABLE contigs AS SELECT * FROM read_fastx('genome.fna', include_filepath := true);
CREATE TABLE genes AS SELECT * FROM predict_genes('contigs', genome_id := 'filepath');

-- prodigal -a
COPY (FROM prodigal_proteins('genes')) TO 'genes.faa' (FORMAT FASTA, INCLUDE_COMMENT true, LINE_WIDTH 60);

-- prodigal -d
COPY (FROM prodigal_genes('genes', 'contigs')) TO 'genes.fna' (FORMAT FASTA, INCLUDE_COMMENT true, LINE_WIDTH 70);

-- A file per genome, each as the CLI would write it for that genome alone
COPY (FROM prodigal_proteins('genes', per_genome := true)) TO 'faa_by_bin'
(FORMAT FASTA, INCLUDE_COMMENT true, LINE_WIDTH 60, PARTITION_BY (genome_id), FILE_EXTENSION 'faa');
```
