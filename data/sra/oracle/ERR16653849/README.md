# fasterq-dump oracles for ERR16653849

Input: `../../ERR16653849.sra` (176893 bytes, md5 `1c7e643476e94d0eb27fb27ced2c1ca9`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/ERR16653849/ERR16653849`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`ERR16653849.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./ERR16653849.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./ERR16653849.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./ERR16653849.sra
```

fasterq-dump reported:

- `--split-3`: spots read 617, reads read 1,234, reads written 1,234
- `--split-spot`: spots read 617, reads read 1,234, reads written 1,234
- `--split-spot --include-technical`: spots read 617, reads read 1,234, reads written 1,234

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/ERR16653849_1.fastq.gz` | `0685fef31074fef6826c398d94f55d96` | 394034 | 617 |
| split_3 | `split_3/ERR16653849_2.fastq.gz` | `f315e8f506fda2eed35609292929d688` | 394034 | 617 |
| split_spot | `split_spot/ERR16653849.fastq.gz` | `90b2bb9a785320f8684bb24ca931c7e3` | 793004 | 1234 |
| split_spot_tech | `split_spot/ERR16653849.fastq.gz` | `90b2bb9a785320f8684bb24ca931c7e3` | 793004 | 1234 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
