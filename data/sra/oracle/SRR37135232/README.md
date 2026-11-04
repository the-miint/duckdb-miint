# fasterq-dump oracles for SRR37135232

Input: `../../SRR37135232.sra` (81112 bytes, md5 `2c4383f2f8304be08ba9d17ab504c1ce`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR37135232/SRR37135232`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR37135232.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR37135232.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR37135232.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR37135232.sra
```

fasterq-dump reported:

- `--split-3`: spots read 23, reads read 23, reads written 23
- `--split-spot`: spots read 23, reads read 23, reads written 23
- `--split-spot --include-technical`: spots read 23, reads read 23, reads written 23

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR37135232.fastq.gz` | `4392c15950eabaa0be4365dbf62d2b1b` | 55224 | 23 |
| split_spot | `split_spot/SRR37135232.fastq.gz` | `0d9798bf1ad84d226f554b59a5ea8300` | 55316 | 23 |
| split_spot_tech | `split_spot/SRR37135232.fastq.gz` | `0d9798bf1ad84d226f554b59a5ea8300` | 55316 | 23 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
