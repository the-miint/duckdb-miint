# fasterq-dump oracles for SRR34459993

Input: `../../SRR34459993.sra` (69189 bytes, md5 `ffa6b5b1b1e07271c9ca9c73213572b9`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR34459993/SRR34459993`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR34459993.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR34459993.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR34459993.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR34459993.sra
```

fasterq-dump reported:

- `--split-3`: spots read 39, reads read 78, reads written 78
- `--split-spot`: spots read 39, reads read 78, reads written 78
- `--split-spot --include-technical`: spots read 39, reads read 78, reads written 78

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR34459993_1.fastq.gz` | `8e39a9cd121918a05d788bd5b9d54f73` | 8348 | 39 |
| split_3 | `split_3/SRR34459993_2.fastq.gz` | `e33ef18d6cd000d545d3e176092db7ed` | 8348 | 39 |
| split_spot | `split_spot/SRR34459993.fastq.gz` | `b8119dd0a3dcf024fc3ae29fdf8e7e5b` | 17008 | 78 |
| split_spot_tech | `split_spot/SRR34459993.fastq.gz` | `b8119dd0a3dcf024fc3ae29fdf8e7e5b` | 17008 | 78 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
