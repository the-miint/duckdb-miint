# fasterq-dump oracles for SRR4424282

Input: `../../SRR4424282.sra` (80940 bytes, md5 `b7359ba7b6f04ded74dc50b6bfd02ee3`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4424282/SRR4424282`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR4424282.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR4424282.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR4424282.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR4424282.sra
```

fasterq-dump reported:

- `--split-3`: spots read 3, reads read 6, reads written 6
- `--split-spot`: spots read 3, reads read 6, reads written 6
- `--split-spot --include-technical`: spots read 3, reads read 6, reads written 6

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR4424282_1.fastq.gz` | `242182aeff5c051b3d57e85423590474` | 1668 | 3 |
| split_3 | `split_3/SRR4424282_2.fastq.gz` | `9911ee75cfe782e13afd475d51cb7d64` | 1668 | 3 |
| split_spot | `split_spot/SRR4424282.fastq.gz` | `e6d7988f40790c56c0fac57ddf5c7745` | 3360 | 6 |
| split_spot_tech | `split_spot/SRR4424282.fastq.gz` | `e6d7988f40790c56c0fac57ddf5c7745` | 3360 | 6 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
