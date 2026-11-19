# HMMER profile fixtures

Inputs for the `read_hmm` tests: profile HMMs from HMMER's and Easel's own source distributions, copied
from the vendored fork (`ext/hmmer`, the-miint/hmmer at `3a9f83bd`, which is HMMER 3.4 with Easel 0.49).
Each test's expected header values are as written in these files. Where HMMER computes a value as it reads
a model (the consensus of a model older than HMMER3/e, and every statistic of a HMMER2 model), the expected
value is what the fork's `hmmconvert` prints for the same file.

## Licence and attribution

The files are distributed as HMMER and Easel distribute them:

- from HMMER's `tutorial/` and `testsuite/`: HMMER's 3-clause BSD licence, `ext/hmmer/LICENSE`
  (Copyright (C) 1992-2023 Sean R. Eddy, President and Fellows of Harvard College, Howard Hughes Medical
  Institute, Washington University School of Medicine, and MRC Laboratory of Molecular Biology);
- from Easel's `demotic/examples/`: Easel's BSD licence, `ext/hmmer/easel/LICENSE`, with the same holders.

Several models are built from public family databases, as their accessions show: `PF…` models from Pfam,
`MADE1` (`DF0000629.2`) from Dfam.

## Files

The larger files are gzipped with `gzip -n` (no name or timestamp in the header); each decompresses to its
source byte for byte.

| File | Source (under `ext/hmmer/`) | Models |
|---|---|---|
| `globins4.hmm.gz` | `tutorial/globins4.hmm` | `globins4`: HMMER3/f, amino; no accession, description or cutoffs |
| `fn3.hmm.gz` | `tutorial/fn3.hmm` | `fn3` (PF00041.13): HMMER3/f, amino; GA, TC and NC |
| `Pkinase.hmm.gz` | `tutorial/Pkinase.hmm` | `Pkinase` (PF00069.17): HMMER3/f, amino; GA, TC and NC |
| `MADE1.hmm.gz` | `tutorial/MADE1.hmm` | `MADE1` (DF0000629.2): HMMER3/f, DNA; MAXL |
| `example.hmm.gz` | `easel/demotic/examples/example.hmm` | `US2`, `Glyco_transf_28`: HMMER3/b, amino; no CONS column, so HMMER computes the consensus |
| `example.hmm2.gz` | `easel/demotic/examples/example.hmm2` | `Orn_DAP_Arg_deC`, `NIF3`: HMMER2.0, amino; HMMER calibrates their statistics as it reads them |
| `XYPPX.hmm` | `testsuite/XYPPX.hmm` | `XYPPX` (PF02162.8): HMMER3/e, amino; GA, TC and NC. Small and plain, so tests derive malformed models from it with `read_text` |
| `XYPPX.hmm.gz` | the same | The same bytes, gzipped |
| `ecori.hmm2` | `testsuite/i15-hmmconvert.pl`, the `h2_nt_file` here-document | `ecori`: HMMER2.0, nucleic (DNA to HMMER 3). Plain, for the same reason as `XYPPX.hmm` |

`ecori.hmm2` was extracted with

```sh
awk '/^sub h2_nt_file/{w=1} /^sub h2_aa_file/{w=0} /<<"EOF"/{on=1; next} /^EOF$/{on=0} on && w' \
    ext/hmmer/testsuite/i15-hmmconvert.pl > data/hmmer/ecori.hmm2
```

One file is not a model, and is made here rather than copied: `truncated_fasta.fa.gz`, a gzipped FASTA
record of 1.2 MB with the last 200 bytes of the gzip stream cut off, so that reading it to the end fails. It
was written with

```sh
python3 -c 'import gzip; body = ">seq1 not a model\n" + ("ACGT" * 15 + "\n") * 20000; data = gzip.compress(body.encode(), mtime=0); open("data/hmmer/truncated_fasta.fa.gz", "wb").write(data[:-200])'
```
