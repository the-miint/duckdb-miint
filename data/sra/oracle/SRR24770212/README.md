# fasterq-dump oracles for SRR24770212

Input: `../../SRR24770212.sra` (394073 bytes, md5 `afe84d2cf5a3e69f7db68d891dff7ed3`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR24770212/SRR24770212`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR24770212.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR24770212.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR24770212.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR24770212.sra
```

fasterq-dump reported:

- `--split-3`: spots read 1,618, reads read 3,153, reads written 3,153
- `--split-spot`: spots read 1,618, reads read 3,153, reads written 3,153
- `--split-spot --include-technical`: spots read 1,618, reads read 3,153, reads written 3,153

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR24770212.fastq.gz` | `c08cacbe6b47450d2bf180dc15c8acfe` | 34424 | 83 |
| split_3 | `split_3/SRR24770212_1.fastq.gz` | `556194e5981c008f078588fa7469f092` | 636450 | 1535 |
| split_3 | `split_3/SRR24770212_2.fastq.gz` | `d8384500e351cfcb76056d270f8a7209` | 636450 | 1535 |
| split_spot | `split_spot/SRR24770212.fastq.gz` | `6dc5b3434e58ea00e92ca2704cb03bf6` | 1319936 | 3153 |
| split_spot_tech | `split_spot/SRR24770212.fastq.gz` | `6dc5b3434e58ea00e92ca2704cb03bf6` | 1319936 | 3153 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
