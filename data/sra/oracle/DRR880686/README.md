# fasterq-dump oracles for DRR880686

Input: `../../DRR880686.sra` (296173 bytes, md5 `80c38c406aaeb7716cf7204fe71092c8`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/DRR880686/DRR880686`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`DRR880686.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./DRR880686.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./DRR880686.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./DRR880686.sra
```

fasterq-dump reported:

- `--split-3`: spots read 261, reads read 261, reads written 261
- `--split-spot`: spots read 261, reads read 261, reads written 261
- `--split-spot --include-technical`: spots read 261, reads read 261, reads written 261

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/DRR880686.fastq.gz` | `6ddf28c16c271bf8f8e0272af402cf87` | 518484 | 261 |
| split_spot | `split_spot/DRR880686.fastq.gz` | `8a1e7891b8801e97468bb84c43927ff9` | 519528 | 261 |
| split_spot_tech | `split_spot/DRR880686.fastq.gz` | `8a1e7891b8801e97468bb84c43927ff9` | 519528 | 261 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
