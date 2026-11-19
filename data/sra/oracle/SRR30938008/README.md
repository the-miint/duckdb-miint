# fasterq-dump oracles for SRR30938008

Input: `../../SRR30938008.sra` (56910 bytes, md5 `632e0fee49934f8fd03a27c8dd025d57`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30938008/SRR30938008`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR30938008.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR30938008.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR30938008.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR30938008.sra
```

fasterq-dump reported:

- `--split-3`: spots read 1, reads read 1, reads written 1
- `--split-spot`: spots read 1, reads read 1, reads written 1
- `--split-spot --include-technical`: spots read 1, reads read 1, reads written 1

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR30938008.fastq.gz` | `4d24665bf97dd2253088c620a1106002` | 1210 | 1 |
| split_spot | `split_spot/SRR30938008.fastq.gz` | `9dcc54916f2c7ea545d362beffb61c8f` | 1214 | 1 |
| split_spot_tech | `split_spot/SRR30938008.fastq.gz` | `9dcc54916f2c7ea545d362beffb61c8f` | 1214 | 1 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
