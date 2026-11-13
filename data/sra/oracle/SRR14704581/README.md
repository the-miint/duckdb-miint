# fasterq-dump oracles for SRR14704581

Input: `../../SRR14704581.sra` (58870 bytes, md5 `e1611efe0cba216262aa5c142fb759a9`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR14704581/SRR14704581`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR14704581.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR14704581.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR14704581.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR14704581.sra
```

fasterq-dump reported:

- `--split-3`: spots read 1, reads read 2, reads written 2
- `--split-spot`: spots read 1, reads read 2, reads written 2
- `--split-spot --include-technical`: spots read 1, reads read 2, reads written 2

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR14704581_1.fastq.gz` | `d858219df278550a25547f823e462ea6` | 928 | 1 |
| split_3 | `split_3/SRR14704581_2.fastq.gz` | `bfb80c717a2e76546148a5302110ad75` | 932 | 1 |
| split_spot | `split_spot/SRR14704581.fastq.gz` | `8322cecf5ffec5e7a7e7697e4a205b4d` | 1868 | 2 |
| split_spot_tech | `split_spot/SRR14704581.fastq.gz` | `8322cecf5ffec5e7a7e7697e4a205b4d` | 1868 | 2 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
