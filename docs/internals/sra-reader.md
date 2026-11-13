# The SRA reader

`read_sra` (`src/read_sra.cpp`) reads NCBI's SRA archives with its own decoder, written with NCBI's
public-domain sources read as its reference. Most of it is written from the format; the parts that
are ports of NCBI code or follow it closely are marked in the code and listed in
`THIRD_PARTY_LICENSES.md`. Nothing of NCBI's is linked or vendored.
Its contract is `fasterq-dump`'s output, byte for byte, and the tests hold it to committed oracles
that are never edited (`data/sra/README.md`).

## Layers

Each layer is free of DuckDB and has a Catch2 test file of its own (`test/cpp/test_Sra*.cpp`), so
the format can be tested against fixtures without a database.

| File | What it does |
|---|---|
| `sra_source` | `SraByteSource`: positional reads from a file, memory, or (in `read_sra.cpp`) a DuckDB `FileHandle`. |
| `sra_archive` | The KAR container (header and table of contents), a KDB database or bare table inside it, the `md/cur` metadata tree, and static columns, whose single value is kept in metadata. |
| `sra_column` | A physical column: its index (`idx`, the `idx0` journal, `idx1`/`idx2` block types), blob location, CRC32 blob checksums (always verified; MD5 blob checksums, which no public run in the survey used, are refused), and the VDB blob framing: headers and page maps that give each row's length. |
| `sra_codecs` | Column encodings: `zip` (deflate) and `izip`/`irzip` (integer planes with delta and two-series modes). |
| `sra_name` | Spot names: the `skey` text index (a persisted trie with a projection from rows to keys) and `format_spot_name`. |
| `sra_run` | Classification (which runs are read, and the message for each kind refused), `SpotReader` (the logical columns as the run's schema defines them, spot by spot, with `fasterq-dump`'s checks; this is where 2na bases are unpacked and ORed with 4na `ALTREAD`), and routing (`RouteSpot`, `RowsForSpot`: which reads `fasterq-dump` writes and where). |

`src/read_sra.cpp` is the DuckDB layer: binding, the work-unit queue (`GlobalState::Claim`), and
filling the output vectors.

## Things that are easy to get wrong

- **Classification goes by the tables and columns present, never by platform.** PacBio and Nanopore
  runs loaded from FASTQ are ordinary sequence tables and are read; their consensus databases are
  refused. Anything outside the table types and encodings the survey of public runs found is
  refused with `SraUnsupportedError` rather than guessed at. Refusals are `InvalidInputException`;
  damage (`SraFormatError`) is an `IOException`.
- **A logical column is chosen as the schema chooses it**: the first alternative whose inputs exist
  (`ORIGINAL_QUALITY` before `QUALITY`, `READ` as 2na bases ORed with `ALTREAD`, a column from its
  directory or else from a static value). Spot names follow the same rule; the comment on
  `SpotReader::Impl::ResolveNaming` lists the alternatives.
- **Order.** Output is in run and spot order with parallel decoding: each work unit (a range of
  spots, `miint_sra_spots_per_unit`, default 131072) gets a batch index when claimed, and
  `get_partition_data` hands it to DuckDB. Claims happen in order, so batch indexes do too.
- **The name index is shared.** Each decoding thread has its own `SpotReader`, but a run's `skey`
  index (tens of bytes per spot when a loader stores whole names) is read once and shared through
  `share_from`; `SkeyIndex` is immutable after construction for that reason.

## Measurements

On 2026-11-05, `read_sra` read SRR7761321 (1.35 GB `.sra`, 20.25 million spots) from a local file
into a table in 6-8 s with 8 threads; `fasterq-dump` took about 131 s on the same file.

## Tests

- C++ (`[sra]`): each layer against the fixtures, with expected values read independently, and
  every fixture's spots, routing and names against the `fasterq-dump` oracles.
- `read_sra.test`: `read_fastx` over the oracles equals `read_sra`, in both modes.
- `read_sra_fastq_render.test`: the rows rendered as FASTQ match the md5 of every oracle file.
- `read_sra_real_data.test` (`SRA_REAL_DATA`): two whole real runs against `fasterq-dump`'s md5s.
