# Prodigal fixtures

Every file is gzipped (`gzip -n -9`). The checksums below are of the decompressed files.

From [the-miint/prodigal](https://github.com/the-miint/prodigal) at `8d7022715c06830db2c162015b2787748744d735`, the
commit gpl-boundary v0.5.0 builds Prodigal 2.6.5 from:

| File | Source path | What it is |
|---|---|---|
| `anthus_aco.fas.gz` | `anthus_aco.fas` | 22 contigs of 960 bp (CRLF line endings) |
| `ref_single.gff.gz`, `ref_single.proteins.gz` | `testdata/ground_truth/` | `prodigal -p single -f gff` / `-a` on `anthus_aco.fas` |
| `ref_meta.gff.gz`, `ref_meta.proteins.gz` | `testdata/ground_truth/` | `prodigal -p meta -f gff` / `-a` on `anthus_aco.fas` |

The reference files were made by upstream Prodigal 2.6.3, so their GFF names it (`Prodigal_v2.6.3`).
The fork's `docs/testing.md` states its library and CLI reproduce them; gpl-boundary's
`src/tools/prodigal.rs` checks its genes and contigs against them column by column.

The GFFs have had their carriage returns removed (`tr -d '\r'`). Prodigal copies each FASTA header whole into the
`# Sequence Data` comment, so on `anthus_aco.fas`'s CRLF headers it writes `seqhdr="61430_aco\r"`, which `read_gff`
and `read_csv` refuse ([#396](https://github.com/the-miint/duckdb-miint/issues/396)). The upstream files' SHA-256,
before that:

```
66afb33984d90fae9415a211cdbd257154004f3d19d30e9b0e5c17049fee28ed  ref_meta.gff
44ef5d1a95c2594e82d3224cc3d183e4a955e26884df4e7d9f8858381024e830  ref_single.gff
```

## Made here

Made with the same upstream Prodigal 2.6.3 (bioconda's `prodigal=2.6.3`), which reproduces `ref_single.gff` byte for
byte, to check the single-genome options against the CLI. The GFFs' carriage returns are removed as above.

| File | Command |
|---|---|
| `ref_single_c.gff.gz`, `ref_single_c.proteins.gz` | `prodigal -i anthus_aco.fas -p single -c -f gff -o ref_single_c.gff -a ref_single_c.proteins` |
| `ref_single_g4.gff.gz`, `ref_single_g4.proteins.gz` | `prodigal -i anthus_aco.fas -p single -g 4 -f gff -o ref_single_g4.gff -a ref_single_g4.proteins` |
| `anthus_masked.fas.gz` | `anthus_aco.fas` with bases 401..460 of each contig replaced by `N` (LF line endings); on `anthus_aco.fas` itself `-m` changes nothing |
| `ref_masked.gff.gz` | `prodigal -i anthus_masked.fas -p single -f gff -o ref_masked.gff` |
| `ref_masked_m.gff.gz`, `ref_masked_m.proteins.gz` | `prodigal -i anthus_masked.fas -p single -m -f gff -o ref_masked_m.gff -a ref_masked_m.proteins` |
| `ref_single.genes.gz` | `prodigal -i anthus_aco.fas -p single -d ref_single.genes` (this Prodigal reproduces `ref_single.proteins` byte for byte) |
| `ref_meta.genes.gz` | `prodigal -i anthus_aco.fas -p meta -d ref_meta.genes` (this Prodigal reproduces `ref_meta.proteins` byte for byte) |

`predict_genes`' own output on `anthus_aco.fas`, for `copy_gff_prodigal.test` to run without gpl-boundary
(`predict_genes.test` checks it is still `predict_genes`' output):

| File | Query |
|---|---|
| `anthus_single_genes.csv.gz` | `predict_genes(..., genome_id := 'genome')`, in the file's contig order |
| `anthus_meta_genes.csv.gz` | `predict_genes(..., genome_id := 'genome', meta_mode := true)`, in the file's contig order |

```sql
CREATE TABLE anthus AS
SELECT 'anthus' AS genome, sequence_index, read_id, sequence1 FROM read_fastx('anthus_aco.fas.gz');
COPY (SELECT g.* FROM predict_genes('anthus', genome_id := 'genome') g JOIN anthus USING (read_id)
      ORDER BY sequence_index, gene_index) TO 'anthus_single_genes.csv.gz';
```

## SHA-256 (decompressed)

```
42904913b6ebc0d40b68a953bc0b2a60d2dd4ffac6424260c90e960ac03641f8  anthus_aco.fas
75cb3638fd39c2f16d281e301cdf368a5b5968c8df6d7448ad7ffcad4c9d8ef6  anthus_masked.fas
766a066c4d4f3cb118343980bc910b2eda8cd550d6795bc2161d7e5c4ff3a365  anthus_meta_genes.csv
c41ab7dd508065edd384a288bd2570beff9eaded73ebd2ef829f7bf11e498c88  anthus_single_genes.csv
ed4d5259945d09194b2ef5aae4606bf713df27256a7a0d10145da1000a972a41  ref_masked_m.gff
615fcd75a2dde0d51658585428d1e3b5b7301799d24685c241c81e3e7f3e3529  ref_masked_m.proteins
6a51e610470c41fffef1a50f4c7e1413a560953daadf265873bbf9d365f5653c  ref_masked.gff
bfd3e41dea798f524a7b4665d8548e244e1cabeeafe4abf9d77a6084ac91348c  ref_meta.genes
86d1f55762754b579d8acf6f28f81204da9ce633b61637c4ce8fecec9c49bc29  ref_meta.gff
a662d74385723c7771892c82c5a08620f4b9d86b679cd8d397a0c8123a20a4b6  ref_meta.proteins
5391437b9aa8e96bb24fe3c103282ddca5a4bc1c5a109369ad532a37fefe263c  ref_single_c.gff
361684b4cef4efca8a16d843f53d9546440ba656765c1ade66d1eb5a2be0f858  ref_single_c.proteins
581d4c6e73ca0afa69b986f212c64a771fc16c6f355378f24030d56aefa3e14c  ref_single_g4.gff
9c080f5febcbb784047b635e1f281777656aa1e23c910a29e1863fb9af9e661a  ref_single_g4.proteins
36c262e4bd3533de768ea68c1f67aca62166412e0d83554570026d010a81b299  ref_single.genes
c7d51315fd61e7505b7d5a328f8b55eb2b71824fa04c19b7567d30be2b0a5eba  ref_single.gff
1d803c1d922cff73a1e98c8879ff392e65930908d8ca6c48213f4b5dbcafa1dc  ref_single.proteins
```
