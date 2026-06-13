# Sylph test fixtures

Small, deterministic fixtures used by `test/sql/sylph_profile.test` and
`test/cpp/test_Sylph*.cpp`. All inputs come from upstream sylph's own
`test_files/` directory and are redistributable under sylph's MIT/Apache-2.0
license.

## Files

| File | Source | Purpose |
|---|---|---|
| `tiny_refs/e.coli-EC590.fasta.gz` | sylph 0.9.0 `test_files/` | Reference genome 1 |
| `tiny_refs/e.coli-K12.fasta.gz`   | sylph 0.9.0 `test_files/` | Reference genome 2 (the matching one) |
| `tiny_refs/e.coli-o157.fasta.gz`  | sylph 0.9.0 `test_files/` | Reference genome 3 |
| `tiny_reads_R1.fq.gz`             | gzip of sylph 0.9.0 `test_files/k12_R1.fq` | Synthetic paired reads from K12 (R1) |
| `tiny_reads_R2.fq.gz`             | gzip of sylph 0.9.0 `test_files/k12_R2.fq` | Synthetic paired reads from K12 (R2) |
| `tiny.syldb`                      | regenerated (see below) | Sylph reference database over the 3 refs |
| `expected_profile.tsv`            | regenerated (see below) | Golden output of `sylph profile` |
| `tiny_oracle.submodule.sha`       | `git rev-parse HEAD` of the embedded sylph fork that produced the above | Drift detector |

## Expected behaviour

`tiny_reads_R{1,2}.fq.gz` are sylph's bundled synthetic K12 paired reads.
Profiling against the 3-genome syldb is expected to identify K12 only, with
~100 % taxonomic and sequence abundance. EC590 and o157 should fall below
the default `min_ani=0.95` threshold despite being closely related E. coli
strains — this validates that sylph's containment-ANI cutoff works as advertised.

## Regenerating

```bash
# Build the embedded fork (Phase 1+; not yet wired):
bash build.sh

# Regenerate fixtures:
bash tools/sylph_oracle.sh

# Pin the new oracle:
( cd ext/sylph && git rev-parse HEAD ) > data/sylph/tiny_oracle.submodule.sha
```

`run_tests.sh` cross-checks `tiny_oracle.submodule.sha` against the embedded
fork's actual HEAD and refuses to run regression tests against a stale oracle.
Same pattern as `data/sortmerna/real_oracle.submodule.sha`.

## Why these inputs

We reuse sylph's own bundled test data instead of synthesizing our own because:

1. They already test the full single-genome-mostly-matches code path that we
   need to assert.
2. They're small (~5 MB total compressed) and redistributable.
3. Reusing them gives upstream-sylph regressions a higher chance of
   showing up here.

For multi-genome abundance tests (Phase 5+ real-world tests), we'll graduate
to a separate, larger gated corpus.
