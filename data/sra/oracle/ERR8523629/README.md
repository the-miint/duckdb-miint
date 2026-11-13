# fasterq-dump oracles for ERR8523629

Input: `../../ERR8523629.sra` (92582 bytes, md5 `a0ed2a7f70dfcb2bb257445882ca7380`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR8523629/ERR8523629`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`ERR8523629.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./ERR8523629.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./ERR8523629.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./ERR8523629.sra
```

fasterq-dump reported:

- `--split-3`: spots read 312, reads read 312, reads written 312
- `--split-spot`: spots read 312, reads read 312, reads written 312
- `--split-spot --include-technical`: spots read 312, reads read 312, reads written 312

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/ERR8523629.fastq.gz` | `9ec8ad495643bad22f5cc8daaa9dc99f` | 135518 | 312 |
| split_spot | `split_spot/ERR8523629.fastq.gz` | `d8c88e89f70a21af19a18fc688936c19` | 136766 | 312 |
| split_spot_tech | `split_spot/ERR8523629.fastq.gz` | `d8c88e89f70a21af19a18fc688936c19` | 136766 | 312 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
