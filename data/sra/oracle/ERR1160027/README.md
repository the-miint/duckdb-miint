# fasterq-dump oracles for ERR1160027

Input: `../../ERR1160027.sra` (37470 bytes, md5 `47c5df1ac97221102ecc53820d323f46`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR1160027/ERR1160027`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`ERR1160027.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./ERR1160027.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./ERR1160027.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./ERR1160027.sra
```

fasterq-dump reported:

- `--split-3`: spots read 71, reads read 71, reads written 71
- `--split-spot`: spots read 71, reads read 71, reads written 71
- `--split-spot --include-technical`: spots read 71, reads read 71, reads written 71

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/ERR1160027.fastq.gz` | `b751c63ef4b1779b26da9001a4e13d28` | 21594 | 71 |
| split_spot | `split_spot/ERR1160027.fastq.gz` | `00d510dac5ee9809ffbed260c141051a` | 21878 | 71 |
| split_spot_tech | `split_spot/ERR1160027.fastq.gz` | `00d510dac5ee9809ffbed260c141051a` | 21878 | 71 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
