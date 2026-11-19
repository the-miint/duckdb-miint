# fasterq-dump oracles for ERR5009897

Input: `../../ERR5009897.sra` (120647 bytes, md5 `72899eca460d1a7e457ea373428b3378`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR5009897/ERR5009897`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`ERR5009897.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./ERR5009897.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./ERR5009897.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./ERR5009897.sra
```

fasterq-dump reported:

- `--split-3`: spots read 248, reads read 496, reads written 248, zero-length reads 248
- `--split-spot`: spots read 248, reads read 496, reads written 248, zero-length reads 248
- `--split-spot --include-technical`: spots read 248, reads read 496, reads written 248, zero-length reads 248

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/ERR5009897.fastq.gz` | `27bcb2cd3c89d1037a9ee7712481ac44` | 111572 | 248 |
| split_spot | `split_spot/ERR5009897.fastq.gz` | `829498c179116653faefed20d61ae80d` | 112564 | 248 |
| split_spot_tech | `split_spot/ERR5009897.fastq.gz` | `829498c179116653faefed20d61ae80d` | 112564 | 248 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
