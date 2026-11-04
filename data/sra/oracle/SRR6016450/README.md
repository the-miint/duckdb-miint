# fasterq-dump oracles for SRR6016450

Input: `../../SRR6016450.sra` (89264 bytes, md5 `b9e97bbb848aa8a4f3265daaa4db740e`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR6016450/SRR6016450`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR6016450.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR6016450.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR6016450.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR6016450.sra
```

fasterq-dump reported:

- `--split-3`: spots read 415, reads read 830, reads written 415, zero-length reads 415
- `--split-spot`: spots read 415, reads read 830, reads written 415, zero-length reads 415
- `--split-spot --include-technical`: spots read 415, reads read 830, reads written 415, zero-length reads 415

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR6016450.fastq.gz` | `f5c8b4d9f9b4e588b6c63388da6908db` | 63478 | 415 |
| split_spot | `split_spot/SRR6016450.fastq.gz` | `ff89f6ddab72fa955bf9d75be00f2d57` | 65138 | 415 |
| split_spot_tech | `split_spot/SRR6016450.fastq.gz` | `ff89f6ddab72fa955bf9d75be00f2d57` | 65138 | 415 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
