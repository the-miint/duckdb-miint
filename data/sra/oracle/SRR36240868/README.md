# fasterq-dump oracles for SRR36240868

Input: `../../SRR36240868.sra` (269306 bytes, md5 `5d5d42933872b62366dbd16b99ff5a5a`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR36240868/SRR36240868`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR36240868.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR36240868.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR36240868.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR36240868.sra
```

fasterq-dump reported:

- `--split-3`: spots read 766, reads read 766, reads written 766
- `--split-spot`: spots read 766, reads read 766, reads written 766
- `--split-spot --include-technical`: spots read 766, reads read 766, reads written 766

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR36240868.fastq.gz` | `2c1a565bd17ad7dfcfb286f73b1ae288` | 984550 | 766 |
| split_spot | `split_spot/SRR36240868.fastq.gz` | `55bf4fbecf108b910bf65a64be60cde1` | 987614 | 766 |
| split_spot_tech | `split_spot/SRR36240868.fastq.gz` | `55bf4fbecf108b910bf65a64be60cde1` | 987614 | 766 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
