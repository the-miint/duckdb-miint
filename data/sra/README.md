# SRA fixtures and fasterq-dump oracles

Fixtures for `read_sra` and `read_ncbi_sra_sequences`. The output of both functions must be byte for byte what
NCBI's `fasterq-dump` writes for the same run, reconstructed spot names included. The oracles in `oracle/` are that
output. They were generated once, outside this repository, and are never regenerated: a change that makes a parity
test fail is a bug in miint, not a stale oracle. miint does not package fasterq-dump, and nothing in this
repository runs it.

## Files

| Path | Contents |
|---|---|
| `<RUN>.sra` | The run's object from NCBI's public AWS bucket, byte for byte |
| `SRR4424282.sralite` | The SRA Lite form of `SRR4424282`, from NCBI's own servers |
| `oracle/<RUN>/split_3/*.fastq.gz` | `fasterq-dump --split-3` output (`<RUN>.fastq` unpaired, `<RUN>_1/_2.fastq` mates), gzipped |
| `oracle/<RUN>/split_spot/<RUN>.fastq.gz` | `fasterq-dump --split-spot` output, gzipped |
| `oracle/<RUN>/README.md` | Provenance for that run's oracles: tool version, exact commands, md5s |
| `oracle/oracle_md5.csv` | Per oracle file: run, mode, file, md5, size and record count of the **uncompressed** FASTQ |
| `real_oracles.csv` | The same for DRR274142 and SRR31286609, whose `.sra` and oracles are not committed (see below) |
| `SHA256SUMS` | sha256 of every other file in this directory; `test/sql/sra_fixtures.test` asserts it |

The oracles are gzipped (`gzip -9 -n`) to keep the repository small: about 3 MB rather than 14 MB. The gzip stream
is not part of the contract. The bytes it decompresses to are, and their md5 is in `oracle_md5.csv`.

`split_spot_tech` rows in `oracle_md5.csv` record `fasterq-dump --split-spot --include-technical`. For every
fixture its output was byte-identical to `--split-spot`, so the rows point at the `split_spot` file and repeat its
md5 rather than commit a second copy. The reason is in "Technical reads" below.

## Accepted fixtures

Every `.sra` was downloaded anonymously on 2026-11-04 from
`https://sra-pub-run-odp.s3.amazonaws.com/sra/<RUN>/<RUN>`. Its size and md5 equal the `size` and `md5` that NCBI's
SRA Data Locator (`https://locate.ncbi.nlm.nih.gov/sdl/2/retrieve?acc=<RUN>`) reported that day for the run's `sra`
file. All are public, and none is human: the organism and BioSample are in the table.

| Run | Bytes | md5 | Spots | Loader | SEQUENCE table type | Platform / layout | Organism | Exercises |
|---|---:|---|---:|---|---|---|---|---|
| SRR32151205 | 628828 | `3b66210bb87d238c1c357c0ccac0c3fa` | 1690 | sharq 3.0.11 | `NCBI:SRA:Illumina:tbl:phred:v2#2` | LS454 (metadata) / paired | mouse gut metagenome | spot names from the `idx/skey` index plus X/Y, a physical `NAME_FMT`, irzip two-series Y, physical ALTREAD (Ns), mates of unequal length |
| SRR34459993 | 69189 | `ffa6b5b1b1e07271c9ca9c73213572b9` | 39 | sharq 3.0.11 | `NCBI:SRA:Illumina:tbl:phred:v2#2` | Illumina / paired | Mus musculus | the same writer, tiny; v1 big-endian READ blobs, two-series X |
| ERR16653849 | 176893 | `1c7e643476e94d0eb27fb27ced2c1ca9` | 617 | sharq 3.4.1 | `NCBI:SRA:Illumina:tbl:phred:v2#2` | Illumina / paired | unidentified | sharq with CRC32 blob checksums |
| ERR5009897 | 120647 | `72899eca460d1a7e457ea373428b3378` | 248 | latf-load 2.9.1 | `NCBI:align:tbl:seq#1.1` | DNBSEQ / paired | soil metagenome | `$R` names, CRC32, physical READ_TYPE, every spot's second read empty, so split-3 sends every spot to the unpaired file |
| SRR4424282 | 80940 | `b7359ba7b6f04ded74dc50b6bfd02ee3` | 3 | latf-load 2.7 | `NCBI:align:tbl:seq#1.1` | Illumina / paired | plant metagenome | ordinary pairs, static READ_LEN, v1 big-endian READ blob |
| SRR6016450 | 89264 | `b9e97bbb848aa8a4f3265daaa4db740e` | 415 | latf-load 2.8.2 | `NCBI:align:tbl:seq#1.1` | Illumina / single | Plasmodium falciparum NF54 | READ_TYPE technical-then-biological; the technical read is empty, so split-spot deflines read `/2` |
| DRR501055 | 83181 | `96b7afb6b1ef186326965f90517ab5d5` | 5 | latf-load 3.0.2 | `NCBI:align:tbl:seq#2` | Illumina / single | metagenome | names from a physical `RAW_NAME`, `ORIGINAL_QUALITY` |
| SRR4439672 | 224214 | `36be3fc900122a227ce97a03fd479518` | 2620 | latf-load 2.7 | `NCBI:align:tbl:seq#1.1` | Illumina / single | Drosophila melanogaster | ALTREAD blobs whose header frame has version 0 |
| SRR24770212 | 394073 | `afe84d2cf5a3e69f7db68d891dff7ed3` | 1618 | sharq 3.0.4 | `NCBI:SRA:GenericFastq:sequence#2` | DNBSEQ / paired | Solanum virginianum | skey names, physical READ_TYPE; split-3 writes paired **and** unpaired spots |
| ERR8523629 | 92582 | `a0ed2a7f70dfcb2bb257445882ca7380` | 312 | fastq-load.py 1.1.2 | `NCBI:SRA:GenericFastq:sequence#1` | Illumina / single | synthetic metagenome | GenericFastq version 1 |
| SRR30810369 | 66332 | `c8b04d15dfca753636347ba1741313d7` | 27 | sharq 3.0.11 | `NCBI:SRA:GenericFastq:sequence#2` | Element / paired | freshwater sediment metagenome | a non-Illumina platform, physical ALTREAD |
| SRR30938008 | 56910 | `632e0fee49934f8fd03a27c8dd025d57` | 1 | fastq-load.py 1.3.5 | `NCBI:SRA:GenericFastq:sequence#2` | BGISEQ / single | Gallus gallus | every column static: READ and quality live in table metadata, with no physical column |
| ERR932865 | 40398 | `b44ed5496db86e19cbf0d0aed33de88c` | 148 | (not recorded) | `NCBI:SRA:Illumina:tbl:phred:v2#1.0.4` | Illumina / single | synthetic construct | a **flat table** (no database), skey names |
| ERR1160027 | 37470 | `47c5df1ac97221102ecc53820d323f46` | 71 | (not recorded) | `NCBI:SRA:Illumina:tbl:phred:v2#1.0.4` | Illumina / single | sponge metagenome | a flat table, names from a physical `SPOT_NAME` |
| SRR14704581 | 58870 | `e1611efe0cba216262aa5c142fb759a9` | 1 | fastq-load.py 1.3 | `NCBI:SRA:GenericFastq:sequence#2` | capillary / paired | Capra hircus | IUPAC ambiguity codes other than N |

### PacBio and Nanopore runs read as sequence tables

fasterq-dump decodes these through the same path as the runs above, because their databases hold no CONSENSUS,
ZMW_METRICS or PASSES table. miint accepts them for the same reason: a run is classified by the tables present, never
by platform (decided 2026-11-04).

| Run | Bytes | md5 | Spots | Loader | SEQUENCE table type | Platform | Organism | Exercises |
|---|---:|---|---:|---|---|---|---|---|
| SRR36240868 | 269306 | `5d5d42933872b62366dbd16b99ff5a5a` | 766 | sharq 3.0.11 | `NCBI:SRA:GenericFastq:sequence#2` | PacBio HiFi | metagenome | qualities of 93 and above, which fasterq-dump prints as `~` |
| DRR880686 | 296173 | `80c38c406aaeb7716cf7204fe71092c8` | 261 | fastq-load.py 1 | `NCBI:SRA:GenericFastq:sequence#2` | Oxford Nanopore | Macaca mulatta | ONT loaded as generic FASTQ |
| SRR37135232 | 81112 | `2c4383f2f8304be08ba9d17ab504c1ce` | 23 | sharq 3.0.11 | `NCBI:SRA:GenericFastq:sequence_nanopore#2` | Oxford Nanopore | indoor metagenome | the table type current ONT loads use |

## Rejection fixtures

miint must refuse these with an error that names what the run is. fasterq-dump reads all of them, so they have no
oracle. What it wrote is recorded so that the refusal is a decision, not an accident.

| File | Bytes | md5 | What it is | fasterq-dump 3.4.1 `--split-3` |
|---|---:|---|---|---|
| SRR26779001.sra | 212860 | `631d6aebf5ded903029d757bf0a7075d` | cSRA: bam-load 3.0.8, tables PRIMARY_ALIGNMENT, REFERENCE, SEQUENCE (wastewater metagenome) | 337 reads. With fewer than 10,000 alignments it treats the run as unaligned and restores the reads from the alignments |
| SRR2513321.sra | 562385 | `ed285e33ed614c616b58c1b15272a446` | ONT with a CONSENSUS table: pore-load.py, `NCBI:SRA:Nanopore:sequence#1` (Saccharomyces cerevisiae W303) | 2 reads, taken from CONSENSUS rather than SEQUENCE (142 spots), and no spot name |
| SRR882123.sra | 302597 | `f55ddde79a98ec4730940e0bee288195` | legacy 454 flat table: sff-load 2.3.2, `NCBI:SRA:_454_:tbl:v2#1.0.7` (marine metagenome) | 238 reads |
| SRR4424282.sralite | 78907 | `0599bb7936ad57889bb597e5e13d99bb` | SRA Lite: `SOFTWARE/delite` 1.1.1 in the **table** metadata, no QUALITY or ORIGINAL_QUALITY column. From `https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos12/sra-pub-zq-1214/SRR004/424/SRR4424282.sralite.1`; size and md5 equal SDL's answer for `filetype=sralite` | 6 reads with every quality `?` (Q30, synthesized) |

## Real runs (downloaded, not committed)

DRR274142 (latf-load, 42,008 spots, 233 of them unpaired under split-3) and SRR31286609 (sharq, 76,368 spots) are
the runs the format was first worked out on. At 15 MB and 30 MB, with 100–200 MB of FASTQ each, they are too large to
commit. `real_oracles.csv` records each object's md5, sha256 and size and the md5, size and record count of every
oracle file, generated the same way as the committed oracles. The tests that use them download the objects into
`data/cache` and are guarded on that download (`run_tests.sh`).

## How the oracles were generated

- Tool: `sra-tools` **3.4.1**, conda package build `h5fa12a8_1` (with `ncbi-vdb` 3.4.1 `haa61c75_0`), from
  bioconda.
- Environment: a throwaway conda environment, created with
  `conda create -p <scratch>/sra-oracle-env --override-channels -c conda-forge -c bioconda "sra-tools=3.4.1=h5fa12a8_1"`.
  It was outside this repository and was discarded.
- Machine: macOS 15.1 (Darwin 24.1.0), arm64.
- Date: 2026-11-04.
- `HOME` pointed at an empty scratch directory, so no user `vdb-config` setting applied. fasterq-dump needed no
  configuration step, and it read only the local file: given a path, it never resolves the accession over the
  network.

Each command ran in the directory holding the fixture, with the fixture named `<RUN>.sra`:

```bash
fasterq-dump --split-3 -O <out> -t <tmp> ./<RUN>.sra
fasterq-dump --split-spot -O <out> -t <tmp> ./<RUN>.sra
fasterq-dump --split-spot --include-technical -O <out> -t <tmp> ./<RUN>.sra
```

fasterq-dump takes the accession it prints (`$ac` in its defline template) from the file name with `.sra` removed,
so every defline begins `@<RUN>.`. `-O` and `-t` only choose where it writes. Every other option was left at its
default: 6 threads, default deflines, technical reads skipped.

fasterq-dump uses one thread for a run of fewer than 10,000 spots, so the committed fixtures say nothing about thread
order. The two real runs are above that limit. Their `--split-3` output with `-e 1` was byte-identical to the default
6-thread output, which confirms that output order is row order and does not depend on the thread count.

## Survey behind the choice

The fixtures were chosen from a structural survey of about 670 small public runs (`mbytes <= 1` in NCBI's AWS
metadata Parquet, sampled across platforms and release years). The survey read each archive's TOC, metadata,
column encodings and blob headers. Findings that shape what miint must support:

- **SEQUENCE table types**, over all 668 runs (370 sampled across platforms, then 298 more Illumina runs from 2020–2024):

  | Table type | Runs |
  |---|---:|
  | `NCBI:align:tbl:seq` (latf-load and bam-load) | 383 |
  | `NCBI:SRA:Illumina:tbl:phred:v2` in a database | 162 |
  | `NCBI:SRA:GenericFastq:sequence` | 64 |
  | `NCBI:SRA:Illumina:tbl:phred:v2` as a **flat table** (no database) | 44 |
  | `_454_`, `IonTorrent`, `Helicos`, `Nanopore`, `sequence_nanopore` and trace together | 15 |

  Flat tables are not only very old runs: SRR14696336 (fastq-load 2.9.1) was released in 2022.
- **Column encodings.** In the accepted table types they are raw, `zip_encoding#1` and `izip_encoding#1`, with one
  exception: ERR263256, a flat Illumina table, stores QUALITY with `delta_average_zip_encoding#1`. Some 454
  columns write the expression with spaces (`< INSDC:quality:phred > izip_encoding#1`).
- **Blob framing.** All 3,027 physical columns have a version 3 header, and none has a compressed idx2 block. Blobs
  are VDB v2, or v1 big-endian for READ. Page maps are versions 1 and 2. Header-frame versions 0–3 all occur;
  version 0 occurs only on `zip_encoding` columns, and no `izip` column used legacy izip.
- **Blob checksums** are CRC32 (1,997 columns) or none (1,030). None is MD5.
- **Non-empty idx0.** 9 of the first 370 runs (about 2.4%) keep blobs in idx0, with nothing in idx1/idx2 for those
  columns. Every one of them was a SARS-CoV-2 or influenza surveillance run from a clinical sample. The 298 non-viral
  runs had none. No idx0 fixture is committed.
- **Not seen in any of the 668 runs:** big-endian archives, chunked or linked TOC entries, and `RNA_FLAG='1'`.
- **Technical reads.** No small run of an accepted table type has a technical read with any bases: every technical
  read found was empty (latf-load's placeholder second read). Technical reads with content occur in the legacy 454
  and IonTorrent tables (key and adapter reads), which are rejected. That is why `--include-technical` changed no
  fixture's output.
- **PacBio and ONT.** Most PacBio runs and nearly all current ONT runs are loaded as `GenericFastq` or `align:tbl:seq`
  tables, with no CONSENSUS, ZMW_METRICS or PASSES table. fasterq-dump reads them through the same path as Illumina
  runs.

## Licence and citation

The `.sra` files are SRA data, distributed by NCBI without restriction under NCBI's data policies; the submitters
may hold rights in the data (https://www.ncbi.nlm.nih.gov/home/about/policies/). The AWS copies are part of the AWS
Open Data registry entry for SRA (NCBI and NIH Genomic Data Sharing policies). sra-tools is a U.S. Government work in
the public domain; see `THIRD_PARTY_LICENSES.md`.

Leinonen R, Sugawara H, Shumway M; International Nucleotide Sequence Database Collaboration. The Sequence Read
Archive. Nucleic Acids Res. 2011;39(Database issue):D19-D21. doi:10.1093/nar/gkq1019
