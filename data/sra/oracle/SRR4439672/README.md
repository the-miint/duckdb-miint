# fasterq-dump oracles for SRR4439672

Input: `../../SRR4439672.sra` (224214 bytes, md5 `36be3fc900122a227ce97a03fd479518`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR4439672/SRR4439672`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR4439672.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR4439672.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR4439672.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR4439672.sra
```

fasterq-dump reported:

- `--split-3`: spots read 2,620, reads read 5,240, reads written 2,620, zero-length reads 2,620
- `--split-spot`: spots read 2,620, reads read 5,240, reads written 2,620, zero-length reads 2,620
- `--split-spot --include-technical`: spots read 2,620, reads read 5,240, reads written 2,620, zero-length reads 2,620

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR4439672.fastq.gz` | `213e9462cd653175c4060288675bf8ee` | 878804 | 2620 |
| split_spot | `split_spot/SRR4439672.fastq.gz` | `6f2e17ad8ce01113a1af4465e96bb87f` | 889284 | 2620 |
| split_spot_tech | `split_spot/SRR4439672.fastq.gz` | `6f2e17ad8ce01113a1af4465e96bb87f` | 889284 | 2620 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
