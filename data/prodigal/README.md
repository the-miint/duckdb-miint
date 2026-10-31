# Prodigal fixtures

Copied unchanged from [the-miint/prodigal](https://github.com/the-miint/prodigal) at
`8d7022715c06830db2c162015b2787748744d735`, the commit gpl-boundary v0.5.0 builds Prodigal 2.6.5 from.

| File | Source path | What it is |
|---|---|---|
| `anthus_aco.fas` | `anthus_aco.fas` | 22 contigs of 960 bp (CRLF line endings) |
| `ref_single.gff`, `ref_single.proteins` | `testdata/ground_truth/` | `prodigal -p single -f gff` / `-a` on `anthus_aco.fas` |
| `ref_meta.gff`, `ref_meta.proteins` | `testdata/ground_truth/` | `prodigal -p meta -f gff` / `-a` on `anthus_aco.fas` |

The reference files were made by upstream Prodigal 2.6.3, so their GFF names it (`Prodigal_v2.6.3`).
The fork's `docs/testing.md` states its library and CLI reproduce them; gpl-boundary's
`src/tools/prodigal.rs` checks its genes and contigs against them column by column.

SHA-256:

```
42904913b6ebc0d40b68a953bc0b2a60d2dd4ffac6424260c90e960ac03641f8  anthus_aco.fas
66afb33984d90fae9415a211cdbd257154004f3d19d30e9b0e5c17049fee28ed  ref_meta.gff
a662d74385723c7771892c82c5a08620f4b9d86b679cd8d397a0c8123a20a4b6  ref_meta.proteins
44ef5d1a95c2594e82d3224cc3d183e4a955e26884df4e7d9f8858381024e830  ref_single.gff
1d803c1d922cff73a1e98c8879ff392e65930908d8ca6c48213f4b5dbcafa1dc  ref_single.proteins
```

## Made here

Not from the fork: made with the same upstream Prodigal 2.6.3 (bioconda's `prodigal=2.6.3`), which reproduces
`ref_single.gff` byte for byte, to check the single-genome options against the CLI.

| File | Command |
|---|---|
| `ref_single_c.gff`, `ref_single_c.proteins` | `prodigal -i anthus_aco.fas -p single -c -f gff -o ref_single_c.gff -a ref_single_c.proteins` |
| `ref_single_g4.gff`, `ref_single_g4.proteins` | `prodigal -i anthus_aco.fas -p single -g 4 -f gff -o ref_single_g4.gff -a ref_single_g4.proteins` |
| `anthus_masked.fas` | `anthus_aco.fas` with bases 401..460 of each contig replaced by `N` (LF line endings); on `anthus_aco.fas` itself `-m` changes nothing |
| `ref_masked.gff` | `prodigal -i anthus_masked.fas -p single -f gff -o ref_masked.gff` |
| `ref_masked_m.gff`, `ref_masked_m.proteins` | `prodigal -i anthus_masked.fas -p single -m -f gff -o ref_masked_m.gff -a ref_masked_m.proteins` |

```
75cb3638fd39c2f16d281e301cdf368a5b5968c8df6d7448ad7ffcad4c9d8ef6  anthus_masked.fas
6a51e610470c41fffef1a50f4c7e1413a560953daadf265873bbf9d365f5653c  ref_masked.gff
ed4d5259945d09194b2ef5aae4606bf713df27256a7a0d10145da1000a972a41  ref_masked_m.gff
615fcd75a2dde0d51658585428d1e3b5b7301799d24685c241c81e3e7f3e3529  ref_masked_m.proteins
804948ba99dc92c32c9ff71afcec79fb1d8af461551b729051c2d4de4a3d2497  ref_single_c.gff
361684b4cef4efca8a16d843f53d9546440ba656765c1ade66d1eb5a2be0f858  ref_single_c.proteins
61e2ec1a6335d06d6c2afe43d53506834559e67b483d4746ad9298e393fd0a09  ref_single_g4.gff
9c080f5febcbb784047b635e1f281777656aa1e23c910a29e1863fb9af9e661a  ref_single_g4.proteins
```
