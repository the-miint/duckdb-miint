# fasterq-dump oracles for SRR32151205

Input: `../../SRR32151205.sra` (628828 bytes, md5 `3b66210bb87d238c1c357c0ccac0c3fa`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR32151205/SRR32151205`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR32151205.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR32151205.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR32151205.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR32151205.sra
```

fasterq-dump reported:

- `--split-3`: spots read 1,690, reads read 3,380, reads written 3,380
- `--split-spot`: spots read 1,690, reads read 3,380, reads written 3,380
- `--split-spot --include-technical`: spots read 1,690, reads read 3,380, reads written 3,380

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR32151205_1.fastq.gz` | `3b1662aecd058c22295c2e55e177ac07` | 1172326 | 1690 |
| split_3 | `split_3/SRR32151205_2.fastq.gz` | `be12b9d845c7a8abafe2b58d226d3c50` | 1159242 | 1690 |
| split_spot | `split_spot/SRR32151205.fastq.gz` | `645ebfa34939dcf2f82eb9b25f966e03` | 2345088 | 3380 |
| split_spot_tech | `split_spot/SRR32151205.fastq.gz` | `645ebfa34939dcf2f82eb9b25f966e03` | 2345088 | 3380 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
