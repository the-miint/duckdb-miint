# Gene prediction

`predict_genes` calls protein-coding genes in prokaryotic contigs, one row per gene. It runs Prodigal 2.6.5, built from [the-miint's fork](https://github.com/the-miint/prodigal) of [Prodigal](https://github.com/hyattpd/Prodigal). Prodigal is GPL-licensed, so it runs in a separate process, [gpl-boundary](internals/embedded-tools.md#gpl-boundary-process-isolation-host-for-gpl-tools), which miint starts for the call, rather than inside miint. Its genes, proteins and per-contig model data equal the Prodigal CLI's, field by field, for `-p single`, `-p meta`, `-c`, `-g 4` and `-m` on the contigs in `data/prodigal/`. To write them as Prodigal's GFF, see [`prodigal_gff_lines`](writing.md#gff).

## predict_genes

```sql
predict_genes(relation, genome_id := 'col', [method := 'prodigal'], [order_by := 'sequence_index'],
              [meta_mode := false], [trans_table := 11], [closed_ends := false], [mask_regions := false],
              [training_path := 'model.trn'])
predict_genes(chunks, mapping_table := 'table', ...)
```

Requires gpl-boundary 0.5.0 or later (`SELECT install_gpl_boundary();`, or `install_gpl_boundary(true)` to replace an older one). Its `prodigal` and `prodigal-train` must report Prodigal 2.6.5 in `gpl-boundary --version`, the version `prodigal_gff_lines` names: `predict_genes` refuses another, since its output would be labeled with the wrong version. Not available on Windows or in WASM.

**Parameters:**

- `relation` (VARCHAR, positional): a table or view of contigs, one per row: `read_id` (VARCHAR, BIGINT or UUID; the contig's name) and `sequence1` (VARCHAR). `read_fastx`'s output is one. Or a table or view of chunked contigs: see [Chunked input](#chunked-input).
- `genome_id` (VARCHAR, **required**): the column that groups contigs into genomes (VARCHAR, BIGINT or UUID). Each distinct non-NULL value is a genome; rows where it is NULL are skipped. `filepath` from `read_fastx(..., include_filepath := true)` makes each FASTA file a genome.
- `method` (VARCHAR, default `'prodigal'`): the gene caller. Prodigal is the only one.
- `order_by` (VARCHAR, default `sequence_index`): an integer column that orders a genome's contigs, as they are given to Prodigal (ties by `read_id`). The default is the column `read_fastx` adds, so a genome's contigs are in file order, as the Prodigal CLI reads them; a relation without it names another column. Training on a genome uses its contigs in this order.
- `meta_mode` (BOOLEAN, default false): false is Prodigal's single-genome mode (`-p single`, the CLI's default): a model is trained on each genome, then its genes are called with it. true is metagenomic mode (`-p meta`): each contig is called with the best of Prodigal's 50 pre-trained models, with no training.
- `trans_table` (INTEGER, default 11): the NCBI translation table a model is trained with, one of those gpl-boundary's `prodigal-train` lists (`gpl-boundary --describe prodigal-train`; for Prodigal 2.6.5, 4, 5, 9, 11, 13, 15, 21, 25, 26 and 32), checked when the call runs. Metagenomic mode takes 11 only (its models choose 11 or 4 per contig).
- `closed_ends` (BOOLEAN, default false): genes may not run off a contig's ends (`-c`).
- `mask_regions` (BOOLEAN, default false): runs of N are treated as gaps, which genes do not cross (`-m`).
- `training_path` (VARCHAR, single-genome mode only): a Prodigal model file (the CLI's `-t`). If it exists, every genome's genes are called with it, and no genome is trained; `trans_table` is then refused (the model holds it), and `closed_ends` and `mask_regions` should be what it was trained with (the model does not store them). If it does not exist, the relation must hold exactly one genome, and the file's directory must exist and be writable; a model is trained on it, written there, and used. A file Prodigal cannot read as a model fails the call.

### Chunked input

Genomes stored in blocks, a row per block of a contig's sequence, are read as they are. The relation has `chunk_data` (VARCHAR or BLOB), `chunk_index` (an integer, 0, 1, 2, ... within a contig) and one key column naming the contig, `read_id` or `feature_idx` (VARCHAR, BIGINT or UUID). A contig's sequence is its rows' `chunk_data` joined in `chunk_index` order. `mapping_table` (VARCHAR, required) names a table or view of the key column and `genome_id` (VARCHAR, BIGINT or UUID) that puts each contig in its genome. In place of `genome_id` and `order_by`:

- A genome's contigs are given to Prodigal in key order.
- `read_id` in the output is the key, in its type.
- Chunks whose key `mapping_table` does not list, or maps to a NULL `genome_id`, are skipped. A key listed twice, a contig whose `chunk_index` values miss or repeat one, and a NULL `chunk_data` fail the call.
- Every contig is joined from its chunks by one `GROUP BY` over the whole relation before the first genome is called, not a genome at a time.

```sql
-- chunks(feature_idx BIGINT, chunk_index INTEGER, chunk_data VARCHAR); bins(feature_idx BIGINT, genome_id VARCHAR)
SELECT * FROM predict_genes('chunks', mapping_table := 'bins');
```

**Output:** one row per gene; one row for each contig without genes, with its contig columns and NULL gene columns; and one row for each genome that failed training or has no genes at all, with only `genome_id`, `status` and `error`.

| Column | Type | Description |
|---|---|---|
| `genome_id` | the input's | The gene's genome |
| `read_id` | the input's | Its contig |
| `gene_index` | INTEGER | 1-based number of the gene on its contig (the `<n>` of the CLI's `<contig>_<n>` protein names) |
| `position` | BIGINT | 1-based start of the gene's leftmost base, on either strand |
| `stop_position` | BIGINT | 1-based **half-open** end (the CLI's end + 1), stop codon included |
| `strand` | VARCHAR | `+` or `-` |
| `partial_left`, `partial_right` | BOOLEAN | The gene runs off the contig's left or right end (the CLI's `partial=`) |
| `start_type` | VARCHAR | `ATG`, `GTG`, `TTG`, or `Edge` for a gene with no start codon |
| `rbs_motif`, `rbs_spacer` | VARCHAR | The ribosome binding site found, NULL where the CLI prints `None` |
| `score` | DOUBLE | `cscore + sscore` (the CLI's `score=`) |
| `cscore`, `sscore`, `rscore`, `uscore`, `tscore` | DOUBLE | Coding, start, RBS, upstream and start-type scores |
| `confidence` | DOUBLE | The CLI's `conf=`, 50 to 100 |
| `gene_gc_content` | DOUBLE | The gene's GC fraction (the CLI's `gc_cont=`) |
| `protein` | VARCHAR | The translation (the CLI's `-a`): `*` for the stop codon, none where the gene runs off an end |
| `contig_length` | BIGINT | The contig's length |
| `contig_gc_content` | DOUBLE | The contig's GC fraction |
| `translation_table` | INTEGER | The table the contig was called with |
| `uses_sd` | BOOLEAN | Whether its model scores Shine-Dalgarno motifs |
| `meta_bin` | INTEGER | Metagenomic mode: the pre-trained model chosen for the contig (0..49); NULL in single-genome mode |
| `model_desc` | VARCHAR | Metagenomic mode: that model's description (the CLI's `model=`); NULL in single-genome mode |
| `model_gc_content` | DOUBLE | The GC fraction of the model the contig was called with (the CLI's Model Data `gc_cont`, which is in percent) |
| `status` | VARCHAR | `ok` on a gene's row; `no_genes` on a contig's row (`read_id` set) and on a genome's own row (`read_id` NULL); `training_failed` on a genome's own row |
| `error` | VARCHAR | For `training_failed`, Prodigal's reason; NULL otherwise |

**Genomes that fail:** in single-genome mode, Prodigal refuses to train on a genome that is too short (the CLI wants at least 20,000 bases). That genome gets one row, `status = 'training_failed'`, with `error` giving the reason and every other column but `genome_id` NULL, and the other genomes are called as usual. Each contig without genes gets a row with `status = 'no_genes'` (so [`prodigal_gff_lines`](writing.md#gff) can number contigs as the CLI does), and a genome with no genes at all also gets a row of its own, `read_id` NULL. Count genes with `WHERE status = 'ok'`. Anything else that goes wrong (gpl-boundary exiting, say) fails the query. Rerunning failed genomes in metagenomic mode is what CheckM does:

```sql
CREATE TABLE contigs AS SELECT * FROM read_fastx('bins/*.fna', include_filepath := true);
CREATE TABLE genes AS SELECT * FROM predict_genes('contigs', genome_id := 'filepath');

CREATE VIEW failed AS
SELECT c.* FROM contigs c JOIN genes g ON c.filepath = g.genome_id WHERE g.status = 'training_failed';
INSERT INTO genes SELECT * FROM predict_genes('failed', genome_id := 'filepath', meta_mode := true);
```

**Errors:** a contig with a NULL `read_id`, a NULL or empty `sequence1`, a non-NULL `sequence2`, or a character outside `A`..`z` (which Prodigal would drop, moving every coordinate after it) fails the call, naming the contig and its genome.

**Order and resources:** genomes come out in `genome_id` order, each genome's genes and gene-less contigs in `order_by` order of their contigs, then the genome's own row if it has one. One genome is in memory at a time (in metagenomic mode, up to 64 MiB of bases of one), and genomes are called one after another in one gpl-boundary process. In single-genome mode, each genome's model is written to a temporary directory (`$TMPDIR/miint-prodigal-*`), removed when the query ends.

```sql
-- Single-genome mode, one genome per FASTA file
SELECT genome_id, count(*) AS genes
FROM predict_genes('contigs', genome_id := 'filepath')
WHERE status = 'ok'
GROUP BY genome_id;

-- A metagenome's contigs
SELECT * FROM predict_genes('assembly', genome_id := 'sample_id', meta_mode := true);

-- Train once on a reference genome, then call other genomes with its model
SELECT * FROM predict_genes('reference', genome_id := 'filepath', training_path := 'ref.trn');
SELECT * FROM predict_genes('others', genome_id := 'filepath', training_path := 'ref.trn');
```
