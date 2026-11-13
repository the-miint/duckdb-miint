# fasterq-dump oracles for DRR501055

Input: `../../DRR501055.sra` (83181 bytes, md5 `96b7afb6b1ef186326965f90517ab5d5`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/DRR501055/DRR501055`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`DRR501055.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./DRR501055.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./DRR501055.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./DRR501055.sra
```

fasterq-dump reported:

- `--split-3`: spots read 5, reads read 10, reads written 5, zero-length reads 5
- `--split-spot`: spots read 5, reads read 10, reads written 5, zero-length reads 5
- `--split-spot --include-technical`: spots read 5, reads read 10, reads written 5, zero-length reads 5

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/DRR501055.fastq.gz` | `080f15c7d780691a8387e9047fd9245c` | 3178 | 5 |
| split_spot | `split_spot/DRR501055.fastq.gz` | `447742a24069178d57ae9891b4940080` | 3198 | 5 |
| split_spot_tech | `split_spot/DRR501055.fastq.gz` | `447742a24069178d57ae9891b4940080` | 3198 | 5 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
