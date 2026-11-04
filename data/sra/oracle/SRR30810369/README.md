# fasterq-dump oracles for SRR30810369

Input: `../../SRR30810369.sra` (66332 bytes, md5 `c8b04d15dfca753636347ba1741313d7`), downloaded 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/SRR30810369/SRR30810369`. Size and md5 match NCBI's SRA Data Locator.

Generated 2026-11-04 on macOS 15.1 (arm64) with fasterq-dump from bioconda `sra-tools` 3.4.1 (build `h5fa12a8_1`,
`ncbi-vdb` 3.4.1), in a throwaway conda environment outside this repository, with `HOME` set to an empty directory.
The environment and the commands are in `../../README.md`. Each command ran in the directory holding
`SRR30810369.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./SRR30810369.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./SRR30810369.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./SRR30810369.sra
```

fasterq-dump reported:

- `--split-3`: spots read 27, reads read 54, reads written 54
- `--split-spot`: spots read 27, reads read 54, reads written 54
- `--split-spot --include-technical`: spots read 27, reads read 54, reads written 54

The files here are that output compressed with `gzip -9 -n`. The md5, size and record count are those of the
uncompressed FASTQ, and `../oracle_md5.csv` holds the same values.

| Mode | File | md5 (uncompressed) | Bytes | Records |
|---|---|---|---:|---:|
| split_3 | `split_3/SRR30810369_1.fastq.gz` | `021097e50a13ec0d10c547e8f0fec09b` | 11754 | 27 |
| split_3 | `split_3/SRR30810369_2.fastq.gz` | `07fa330843eb099de01d2453f90980b9` | 11754 | 27 |
| split_spot | `split_spot/SRR30810369.fastq.gz` | `65332ef1fb364aa5c2192e646d46220f` | 23724 | 54 |
| split_spot_tech | `split_spot/SRR30810369.fastq.gz` | `65332ef1fb364aa5c2192e646d46220f` | 23724 | 54 |

`split_spot_tech` is fasterq-dump's `--include-technical` output. It was byte-identical to `--split-spot`, so it
points at the same file. These oracles are fixed: a parity test that disagrees with them is a bug in miint.
