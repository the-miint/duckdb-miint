# fasterq-dump oracles for ERR932865

Input: `../../ERR932865.sra` (40398 bytes, md5 `b44ed5496db86e19cbf0d0aed33de88c`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR932865/ERR932865`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`ERR932865.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./ERR932865.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./ERR932865.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./ERR932865.sra
```

fasterq-dump reported:

- `--split-3`: spots read 148, reads read 148, reads written 148
- `--split-spot`: spots read 148, reads read 148, reads written 148
- `--split-spot --include-technical`: spots read 148, reads read 148, reads written 148

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/ERR932865.fastq.gz` | `9c70c3eaca391954f68ca81ade1024bb` | 31016 | 148 |
| split_spot | `split_spot/ERR932865.fastq.gz` | `3094d9795b61ea25af7272bbecf1187a` | 31608 | 148 |
| split_spot_tech | `split_spot/ERR932865.fastq.gz` | `3094d9795b61ea25af7272bbecf1187a` | 31608 | 148 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
