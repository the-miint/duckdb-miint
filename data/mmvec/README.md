# MMvec test fixtures

Input data for the `mmvec` implementation's unit tests and sqllogictests. Expected
values are **not** here — they are carved into `test/cpp/mmvec_oracle.hpp`, which also
records the oracle environment and per-tier tolerances.

Every count fixture is a long-form COO feature table with the schema
`(sample_id VARCHAR, feature_id VARCHAR, value DOUBLE)` — the same schema
[`read_biom`](../../docs/reading.md) emits and every diversity function consumes. Zero
cells are omitted (the sparse feature-table contract). All 14 count fixtures were
verified to contain no all-zero rows, no all-zero columns, no NULLs, no stored zeros and
no non-finite values; within each pair the `sample_id` sets are identical.

## Ordering rule (matters for correctness, not just cosmetics)

Consumers must sort feature ids **lexicographically**. MMvec pins the *first* Y feature
as the reference category with a fixed zero logit, and because the Gaussian priors
penalize parameters (which are not softmax-shift-invariant), **the MAP solution depends
on which feature that is**. Measured at tight convergence by permuting Y and realigning:
the choice of reference moves ranks by up to 31% of their magnitude on `toy` and lands on
a genuinely different (lower) optimum. Permuting the *remaining* Y columns is harmless.

Every oracle value in `test/cpp/mmvec_oracle.hpp` was generated under this rule. For the
`soils` fixture the rule is a no-op — its file order and lexicographic order happen to
agree on the first metabolite — but for `cf` it selects a different reference than the
BIOM's storage order would.

## Fixtures

| File | Shape | nnz | Source |
|---|---|---|---|
| `toy_x` / `toy_y` | 5 × 8 / 5 × 6 | 27 / 20 | the `mmvec` docstring example in scikit-bio 0.7.3 |
| `toy_x_test` / `toy_y_test` | 3 × 8 / 3 × 6 | 19 / 14 | the held-out samples from that docstring's `predict`/`score` examples |
| `synth_a_x` / `synth_a_y` | 8 × 4 / 8 × 6 | 32 / 48 | `random_multimodal(4, 6, 8, 2, seed=7)` |
| `synth_b_x` / `synth_b_y` | 50 × 8 / 50 × 10 | 285 / 485 | `random_multimodal(8, 10, 50, 3, seed=42)` |
| `synth_c_x` / `synth_c_y` | 150 × 8 / 150 × 8 | 1186 / 1200 | `random_multimodal(8, 8, 150, 2, x_noise_std=2, x_total=1000, y_total=10000, seed=1)` |
| `soils_x` / `soils_y` | 19 × 466 / 19 × 85 | 4534 / 1609 | desert biocrust wetting study |
| `cf_x` / `cf_y` | 172 × 2720 / 172 × 462 | 9810 / 18429 | cystic fibrosis sputum study |
| `cf_microbe_meta` | 2720 rows | — | `(feature_id, taxon)`; 39 rows match `Pseudomonas` |
| `cf_metabolite_meta` | 462 rows | — | `(feature_id, expert_annotation)`; 20 rows annotated |

`random_multimodal` is scikit-bio's private test simulator
(`skbio.stats.ordination._mmvec.random_multimodal`). It is deliberately **not** ported to
C++ — its outputs are committed here instead, which avoids having to replicate NumPy's
PCG64 stream and the simulator's unused regression-coefficient draw, which never reaches
the generated data but does advance the RNG stream (an artifact inherited from upstream
`biocore/mmvec`, so not removable without changing every seeded result).

### Fixture checksums

```
c933402112a83e3a1e3f713414da54b7a3ef681049907286d235207f6f3cff64  toy_x.parquet
bbf4c117245fac5b06d6ee96ed5cf1e3634973b9e51280eefb260c9e4b851327  toy_y.parquet
8b48633254f23f431341fed1f17cbec1a6db025c50680e25b6593b427beb32e1  toy_x_test.parquet
7554aeb44b8f69be09df2a69e18ddd1b06ae9c5855fde01275f0bc3e5c6c6e47  toy_y_test.parquet
c172daa4900d34affed8161f71997c6472c029fcbe196adc525cc5a35384c870  synth_a_x.parquet
106b74217a5f30210919cf723d6e74abd9faefb8e1e5c6ecb8fbe9767ef54afd  synth_a_y.parquet
2514eff1481830e6bd9330e7e1f9ccd53375b80ddf07ea010d4fa3a3b6719ce8  synth_b_x.parquet
6114c028bf151b948057bbe93553a63c86afa677768422a93820528280c64600  synth_b_y.parquet
e21c0a3ed5b1c52e7e496b41e4134812e84224810ebddc709a3d41bb53f2342e  synth_c_x.parquet
c4880929eda9c420de5ab4995765025d16ae7feaa9f5e01014b2e528ee508021  synth_c_y.parquet
f219b6e3afca95e5ddc64c5cd566793939a5e2f3900d4d52608bb3a505f819cc  soils_x.parquet
276f4c607599cdb79eb96f27139d2cdcdc081f1181d934bc049388b90145ed93  soils_y.parquet
d38c813fd0b508bdc446714b75245d02a9ac2899534c8444665abd59efd0ffaf  cf_x.parquet
fd5e2d58a05c7f047e3e789407f3a443a2edc9a4030f5ef2744fa3d89a3ed0ea  cf_y.parquet
49586ba574b8f8e0730d818b262c5c602047b99cc35d4cb8775091394973492e  cf_microbe_meta.parquet
24853e3ffdb84d1af1ae1491a80eaef5b388171381f7d8e6288fd1a985e8f8a5  cf_metabolite_meta.parquet
```

## Derivation of the real-data fixtures

### soils

Unpivoted from the two TSVs committed in scikit-bio 0.7.3
(`skbio/stats/ordination/tests/data/soils/`), which are stored features × samples and are
transposed by scikit-bio's own test. Verified against scikit-bio's committed sidecars:

```
69f018c4a907d8548d21bdca9f73f7392bd54ccf0c8d1e738ff13185912a025f  microbes.tsv
802e9076be6e14073e2b3730964a1a20c0c18c27214e40617d9feff7e052c0ae  metabolites.tsv
```

### cf

scikit-bio ships only `.sha256` sidecars for the CF tables; the TSVs themselves were
removed in scikit-bio PR #2448 and are not distributed. The fixtures here are therefore
derived from the upstream BIOM files, which scikit-bio's first MMvec commit
(`3b445f3b`) did commit at byte-identical sizes:

| Upstream file (`biocore/mmvec` @ `master:examples/cf/`) | sha256 |
|---|---|
| `otus_nt.biom` (1,617,749 B, HDF5) | `13445d4d1c963737a4de62c753d42bd9ab6f53743f74bf28aed330a973b7cbeb` |
| `lcms_nt.biom` (220,253 B, HDF5) | `97a332c0b893a21d604a6da6bf459a17166e4eff9ffec5f085f9bbcb5c272805` |
| `microbe-metadata.txt` (1,959,705 B) | `ade57e616119dee81e4d8544f2ec3540dfa979102723962a4ad34dcd60b72342` |
| `metabolite-metadata.txt` (224,918 B) | `58281f7568a51cc4055adfb423bc591c70b760e7a8af8f4f285a09d7b374b205` |

Derivation, via `read_biom` plus a sample intersection — exactly the steps scikit-bio's
original CF test performed (intersect common samples, drop zero-sum columns, drop
zero-sum rows):

1. Read both BIOMs. As stored they are 636 samples × 7558 microbes and 180 samples × 462
   metabolites.
2. Restrict both to the **172 samples present in both**.
3. Features that become all-zero simply have no COO rows, leaving **2720 microbes** and
   **462 metabolites**.
4. Project the metadata to the columns the tests use, restricted to surviving features.
   The metabolite metadata's id column is named `sampleid` upstream even though it holds
   metabolite ids; it is renamed `feature_id` here to avoid a genuine footgun.

That reproduces scikit-bio's documented shape (172 × 2720 × 462) and, independently, its
documented metadata counts (`n=39` *Pseudomonas*, `n=20` expert-annotated) — three
numbers agreeing, so this is their data, not an approximation of it.

## Provenance and citations

If you use these datasets, cite the studies they come from.

**Method.** Morton, J. T.; Aksenov, A. A.; Nothias, L. F.; Foulds, J. R.; Quinn, R. A.;
Badri, M. H.; …; Knight, R. (2019) "Learning representations of microbe–metabolite
interactions", *Nature Methods* 16(12), 1306–1314.
doi: [10.1038/s41592-019-0616-3](https://doi.org/10.1038/s41592-019-0616-3)

**soils.** Swenson, T. L.; Karaoz, U.; Swenson, J. M.; Bowen, B. P.; Northen, T. R.
(2018) "Linking soil biology and chemistry in biological soil crust using isolate
exometabolomics", *Nature Communications* 9(1), 19.
doi: [10.1038/s41467-017-02356-9](https://doi.org/10.1038/s41467-017-02356-9)

19 samples over 4 time points after a laboratory wetting event × 4 biocrust successional
stages; microbes by rplO genotyping, metabolites by LC/MS.

**cf.** Quinn, R. A.; Comstock, W.; Zhang, T.; Morton, J. T.; da Silva, R.; Tran, A.;
…; Dorrestein, P. C. (2018) "Niche partitioning of a pathogenic microbiome driven by
chemical gradients", *Science Advances* 4(9), eaau1908.
doi: [10.1126/sciadv.aau1908](https://doi.org/10.1126/sciadv.aau1908)

172 samples, 16S rRNA V4 ASVs (feature ids are the ASV sequences) and LC-MS/MS
metabolites.

Both datasets and the reference implementation are BSD-3-Clause: the data via
[`biocore/mmvec`](https://github.com/biocore/mmvec) and
[`scikit-bio`](https://github.com/scikit-bio/scikit-bio). See
[`THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md) for the full record.
