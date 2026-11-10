# Working with NCBI

MIINT provides integration with NCBI's GenBank and RefSeq.

## Table of contents

- [Accessions are literals, not table names](#accessions-are-literals-not-table-names) - Applies to all three readers below; includes the idiom for driving them from accessions held in a table.
- [Requests: waits and where they go](#requests-waits-and-where-they-go) - Retries, BLAST's waits, what stops them, and the settings that point the requests elsewhere.
- [Read metadata](#read-metadata) - Fetch record metadata from GenBank and RefSeq by accession.
- [Read sequences](#read-sequences) - Fetch sequence data from GenBank and RefSeq by accession.
- [Read annotations](#read-annotations) - Fetch annotation data from GenBank and RefSeq by accession.
- [Taxonomy](#taxonomy) - Obtain and use the NCBI taxonomy: the full tree (bulk/offline), online lineage lookups, offline lineage collapse, and parsing prefixed lineage strings.
  - [Read the taxonomy tree (`read_ncbi_taxdump`)](#read-the-taxonomy-tree-read_ncbi_taxdump) - Bulk/offline: download + cache the taxonomy tree.
  - [Remap retired taxids (`read_ncbi_taxdump_merged`)](#remap-retired-taxids-read_ncbi_taxdump_merged) - The retired→current taxid map.
  - [Read every name class (`read_ncbi_taxdump_names`)](#read-every-name-class-read_ncbi_taxdump_names) - All of `names.dmp`, including `genbank common name`.
  - [List deleted taxids (`read_ncbi_taxdump_deleted`)](#list-deleted-taxids-read_ncbi_taxdump_deleted) - The taxids NCBI deleted outright.
  - [Look up lineages online (`read_ncbi_lineage`)](#look-up-lineages-online-read_ncbi_lineage) - Online/rate-limited: resolve a handful of taxids via E-utilities.
  - [Collapse the tree to lineages (`taxonomy_lineage`)](#collapse-the-tree-to-lineages-taxonomy_lineage) - Offline: rank-collapse a taxdump tree, identical schema to `read_ncbi_lineage`.
  - [Parse lineage strings (`parse_lineage`)](#parse-lineage-strings-parse_lineage) - Turn a `d__…; p__…` lineage string (GTDB, Greengenes2, tax2tree) into the shared schema's rank columns.
  - [Map an NCBI rank to its column (`lineage_rank`)](#map-an-ncbi-rank-to-its-column-lineage_rank) - The shared schema's column an NCBI rank collapses into, so `superkingdom` and `domain` look up alike.
  - [Split GTDB labels (`gtdb_label`, `gtdb_species_label`)](#split-gtdb-labels-gtdb_label-gtdb_species_label) - A GTDB name's base and the suffix GTDB gives the parts of a split taxon (`Bacillota_A`), to find split and polyphyletic taxa.
  - [Read GTDB's genome metadata (`read_gtdb_metadata`, `read_gtdb_qc_failed`)](#read-gtdbs-genome-metadata-read_gtdb_metadata-read_gtdb_qc_failed) - Download and cache a GTDB release's genome metadata and qc_failed list, or read the files; link an NCBI species to its type strain's genome.

### Accessions are literals, not table names

`read_ncbi`, `read_ncbi_fasta` and `read_ncbi_annotation` take **literal** accessions — a VARCHAR or a
VARCHAR[]. Passing the *name of a table* that holds accessions is rejected at bind time:

```sql
CREATE TABLE my_accessions AS SELECT 'NC_001416.1' AS accession;

SELECT * FROM read_ncbi_annotation('my_accessions');
-- Invalid Input Error: read_ncbi_annotation: 'my_accessions' is a table in the catalog,
-- not an NCBI accession. ...
```

To drive one of these functions from accessions held in a table, hoist them into a list first. A
subquery cannot be used as a table-function argument at all (`Binder Error: Table function cannot
contain subqueries`), so use `SET VARIABLE` + `getvariable()`:

```sql
SET VARIABLE accs = (SELECT list(accession) FROM my_accessions);
SELECT * FROM read_ncbi_annotation(getvariable('accs'));
```

The check keys off whether the string resolves to a table or view in the catalog, not off the string's
shape — accession formats are open-ended (bare GenBank ids such as `J02459` carry no distinguishing
prefix), so a pattern test would reject valid input. A name that is *not* in the catalog is still sent
upstream as an accession.

Both bare and qualified names are covered: a bare name is resolved against the session `search_path`,
and a qualified one (`my_schema.accs`, `my_db.my_schema.accs`) is resolved as written. Versioned
accessions are unaffected — `'NC_001416.1'` is *looked up* as schema `NC_001416` / relation `1`, finds
nothing, and is sent upstream unchanged.

### Requests: waits and where they go

The functions here that go online, and `blast`, which searches NCBI BLAST, wait between attempts
and for NCBI:

- `read_ncbi`, `read_ncbi_fasta`, `read_ncbi_annotation` and `read_ncbi_lineage` retry a transient
  failure after 1, 2, 4, 8, 16 and 32 s. `read_ncbi_taxdump` retries its download after 1, 2 and 4 s.
- `blast` backs off between attempts as the readers do. It also waits for NCBI's estimate of how
  long a search takes, at most an hour, then polls every 15 s for up to another hour, after which
  it fails.

Ctrl-C and `max_execution_time` end any of these waits within about a tenth of a second. Before
miint sees the status of a GET answered 429 or 503, DuckDB retries the GET itself, as its
`http_retries`, `http_retry_wait_ms` and `http_retry_backoff` settings say. With the defaults that
takes about a minute, and Ctrl-C does not end it
([#414](https://github.com/the-miint/duckdb-miint/issues/414)). A session `SET` of those three
does not reach miint's requests; use `SET GLOBAL`.

Two settings send the requests elsewhere, such as to a mirror, and a session `SET` changes them.
A trailing slash is dropped, and a query that would use an empty URL fails.

- `miint_ncbi_eutils_url`, the base URL of E-utilities (default
  `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`). The Datasets API, which `read_ncbi_fasta` uses
  for assemblies, keeps its own URL.
- `miint_blast_url`, the URL of BLAST's URL API (default
  `https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi`).

```sql
SET miint_ncbi_eutils_url = 'https://eutils.example.org/entrez/eutils';
RESET miint_ncbi_eutils_url;  -- back to NCBI's
```

### Read metadata

Fetch GenBank and RefSeq metadata from NCBI by accession number. This function queries NCBI's E-utilities API to retrieve sequence metadata without downloading the full sequence.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1', 'NM_001101.5')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits (10 req/s vs 3 req/s)
- `batch_size` (BIGINT, optional, default 500): Number of accessions per epost+efetch round-trip. Larger batches reduce per-accession HTTP overhead; the default matches NCBI's practical guidance for POST.

**Output schema:**
- `accession` (VARCHAR): Accession with version (e.g., 'NC_001416.1')
- `version` (INTEGER): Version number extracted from accession
- `description` (VARCHAR): Sequence description/definition
- `organism` (VARCHAR): Source organism name
- `taxonomy_id` (BIGINT): NCBI taxonomy ID
- `length` (BIGINT): Sequence length in base pairs
- `molecule_type` (VARCHAR): Molecule type (e.g., 'DNA', 'RNA', 'protein')
- `update_date` (DATE): Last modification date

**Behavior:**
- Fetches GenBank XML format from E-utilities using the batched `epost` + `efetch` handshake (one POST + one GET per `batch_size` accessions)
- Automatically detects accession type (RefSeq NC_/NM_/NP_, GenBank, etc.)
- Rate-limited to respect NCBI guidelines (3 req/s without key, 10 req/s with key)
- Retry logic for transient failures (400, 408, 429, 500, 502, 503, 504) with exponential backoff up to 6 retries
- Accessions silently omitted by NCBI (e.g. invalid IDs) are reported via `miint_warnings()` and a stderr `WARNING`; the function does **not** throw, so one bad accession in a batch never aborts the rest of the query

**Examples:**
```sql
-- Get metadata for a single accession
SELECT * FROM read_ncbi('NC_001416.1');

-- Get metadata for multiple accessions
SELECT * FROM read_ncbi(['NC_001416.1', 'NC_001422.1']);

-- Use API key for higher rate limits
SELECT * FROM read_ncbi('NC_001416.1', api_key='your_ncbi_api_key');

-- Get organism and length for a set of sequences
SELECT accession, organism, length
FROM read_ncbi(['NC_001416.1', 'NC_001422.1', 'NC_000913.3']);

-- Filter by molecule type
SELECT accession, description
FROM read_ncbi(['NC_001416.1', 'NM_001101.5'])
WHERE molecule_type = 'DNA';
```

**Notes:**
- Empty accessions raise `InvalidInputException` at bind time
- A table or view name passed where an accession belongs is rejected at bind time — see [Accessions are literals, not table names](#accessions-are-literals-not-table-names)
- Invalid accessions that NCBI silently drops appear in `miint_warnings()` (filter with `SELECT * FROM miint_warnings() WHERE message LIKE '%read_ncbi%'`) — they do not raise
- For bulk downloads, set `api_key` for the higher 10 req/s rate; with the default `batch_size=500` a 10,000-accession query is bound by ~20 round-trips, not by the per-accession rate cap

## Read sequences 

Fetch FASTA sequences from NCBI by accession number. Returns data in the same schema as `read_fastx`, making it easy to combine NCBI sequences with local files.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column showing NCBI URL
- `batch_size` (BIGINT, optional, default 500): Number of sequence accessions per epost+efetch round-trip. Assembly accessions (GCF_/GCA_) are always fetched one at a time via the Datasets API, so this parameter only affects the sequence (NC_/NM_/...) path.

**Output schema (matches `read_fastx`):**
- `sequence_index` (BIGINT): 0-based sequential index
- `read_id` (VARCHAR): Accession number extracted from FASTA header
- `comment` (VARCHAR, nullable): Description from FASTA header
- `sequence1` (VARCHAR): DNA/RNA sequence
- `sequence2` (VARCHAR): Always NULL (single sequences only)
- `qual1` (UINT8[]): Always NULL (FASTA has no quality scores)
- `qual2` (UINT8[]): Always NULL
- `filepath` (VARCHAR, optional): the URL the records were fetched from, when include_filepath=true

**Behavior:**
- Fetches FASTA via the batched `epost` + `efetch` handshake (sequence accessions); assemblies route through the Datasets API one at a time
- Parses pipe-delimited FASTA headers (e.g., `gi|123|ref|NC_001416.1|description`)
- Extracts accession as `read_id`, remainder as `comment`
- Rate-limited and retried on transient failures (400, 408, 429, 5xx) with exponential backoff up to 6 retries
- Accessions silently omitted by NCBI are reported via `miint_warnings()` and stderr `WARNING`; the function does not throw on missing accessions, so one bad ID in a 10,000-accession query never aborts the rest

**Examples:**
```sql
-- Fetch a single sequence
SELECT * FROM read_ncbi_fasta('NC_001416.1');

-- Fetch multiple sequences
SELECT read_id, LENGTH(sequence1) AS seq_length
FROM read_ncbi_fasta(['NC_001416.1', 'NC_001422.1']);

-- Combine NCBI sequences with local files
SELECT read_id, sequence1 FROM read_ncbi_fasta('NC_001416.1')
UNION ALL
SELECT read_id, sequence1 FROM read_fastx('local_sequences.fasta');

-- Export NCBI sequence to local FASTA file
COPY (SELECT * FROM read_ncbi_fasta('NC_001416.1'))
TO 'lambda_phage.fasta' (FORMAT FASTA);

-- Get sequence with source URL tracking
SELECT read_id, LENGTH(sequence1) AS length, filepath
FROM read_ncbi_fasta('NC_001416.1', include_filepath=true);

-- Use as reference for alignment
CREATE TABLE reference AS SELECT * FROM read_ncbi_fasta('NC_001416.1');
CREATE TABLE reads AS SELECT * FROM read_fastx('reads.fastq');
SELECT * FROM align_minimap2('reads', 'reference');
```

**Notes:**
- Large sequences (e.g., complete chromosomes) may take time to download
- Output is compatible with all functions expecting `read_fastx` schema
- Empty accessions raise `InvalidInputException` at bind time
- A table or view name passed where an accession belongs is rejected at bind time — see [Accessions are literals, not table names](#accessions-are-literals-not-table-names)
- Missing/invalid accessions appear in `miint_warnings()` rather than raising — query that table to see what NCBI dropped after a large batched fetch

### Read annotations

Fetch feature annotations from NCBI by accession number. Returns data in the same schema as `read_gff`, making it easy to analyze NCBI annotations with SQL.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI servers (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `accession` (VARCHAR or VARCHAR[]): NCBI accession number(s) (e.g., 'NC_001416.1')
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits
- `include_filepath` (BOOLEAN, optional, default false): Add filepath column showing NCBI URL

**Output schema (matches `read_gff`):**
- `seqid` (VARCHAR): Sequence/chromosome identifier (the accession)
- `source` (VARCHAR): Annotation source ('RefSeq', 'GenBank', or 'NCBI')
- `type` (VARCHAR): Feature type (e.g., 'gene', 'CDS', 'mRNA')
- `position` (INTEGER): 1-based start position
- `stop_position` (INTEGER): 1-based **half-open** end position — i.e. NCBI's annotation `end` **+ 1**, normalized on read so it matches every other miint function. A feature's length is `stop_position - position`, with no `+ 1`. The raw NCBI `end` is `stop_position - 1`. This changed: it previously emitted the closed `end` directly. See [Coordinate conventions](reading.md#coordinate-conventions).
- `score` (DOUBLE, nullable): Feature score (typically NULL for NCBI data)
- `strand` (VARCHAR): Strand ('+' or '-')
- `phase` (INTEGER, nullable): CDS phase (0, 1, 2) derived from codon_start qualifier
- `attributes` (MAP(VARCHAR, VARCHAR)): Key-value attributes (gene, product, locus_tag, etc.)
- `filepath` (VARCHAR, optional): the URL the records were fetched from, when include_filepath=true

**Behavior:**
- Fetches INSDC feature table format from E-utilities
- Converts to GFF3-compatible schema
- Detects source automatically from accession prefix (RefSeq for NC_/NM_, etc.)
- Handles complement strand (reversed positions)
- Parses codon_start qualifier to set correct CDS phase
- Warns on complex locations (join, complement) - uses outer bounds
- A table or view name passed where an accession belongs is rejected at bind time — see [Accessions are literals, not table names](#accessions-are-literals-not-table-names)

**Examples:**
```sql
-- Get all annotations for Lambda phage
SELECT * FROM read_ncbi_annotation('NC_001416.1');

-- Get only genes
SELECT seqid, position, stop_position, attributes['gene'] AS gene_name
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'gene';

-- Get CDS features with their products
SELECT position, stop_position, strand, phase,
       attributes['gene'] AS gene,
       attributes['product'] AS product
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'CDS';

-- Count feature types
SELECT type, COUNT(*) AS count
FROM read_ncbi_annotation('NC_001416.1')
GROUP BY type
ORDER BY count DESC;

-- Combine annotations from multiple genomes
SELECT seqid, type, COUNT(*) AS feature_count
FROM read_ncbi_annotation(['NC_001416.1', 'NC_001422.1'])
GROUP BY seqid, type
ORDER BY seqid, feature_count DESC;

-- Find genes on the minus strand
SELECT attributes['gene'] AS gene_name, position, stop_position
FROM read_ncbi_annotation('NC_001416.1')
WHERE type = 'gene' AND strand = '-';

-- Export to local GFF-like format
COPY (
    SELECT seqid, source, type, position, stop_position,
           NULL AS score, strand, phase
    FROM read_ncbi_annotation('NC_001416.1')
) TO 'lambda_features.tsv' (DELIMITER '\t', HEADER);
```

**Notes:**
- Output is compatible with queries designed for `read_gff`
- Complex feature locations (join, complement) emit a warning and use outer bounds only
- The `phase` column is set from the `codon_start` qualifier for CDS features
- Source is detected from accession prefix: NC_/NM_/NP_ -> 'RefSeq', others -> 'GenBank' or 'NCBI'

## Taxonomy

`read_ncbi(...)` returns a `taxonomy_id` per accession, but a bare taxid is not very
useful on its own. These functions turn taxids into ranks, scientific names, and full
lineages, and let you navigate the NCBI taxonomy as a tree.

There are **two acquisition paths**, chosen by how much of the taxonomy you need:

| | Function | Use when |
|---|---|---|
| **Bulk / offline** | `read_ncbi_taxdump` (+ `taxonomy_lineage`) | You need the whole taxonomy, join against many taxids, or want to work without repeated network calls. Downloads the full dump once (~450 MB extracted) and caches it. |
| **Online / rate-limited** | `read_ncbi_lineage` | You have a handful of taxids and just want their lineages. One E-utilities round-trip per batch; subject to NCBI rate limits. |

`taxonomy_lineage` (offline) and `read_ncbi_lineage` (online) emit the **identical
schema**, so you can develop against the API and switch to the cached dump at scale
without changing downstream SQL.

The NCBI taxonomy is a tree in the exact shape `read_newick` emits, so `read_ncbi_taxdump`
output (`node_index`, `parent_index`, `name`, `is_tip`, plus `rank`) drops straight into
the existing tree tooling.

**Malformed input fails loudly.** All four taxdump readers validate what they parse rather
than coercing it, because a taxonomy that silently reads *slightly* wrong is worse than one
that refuses to read:

- **Field widths are checked exactly** for `names.dmp` (4), `merged.dmp` (2) and
  `delnodes.dmp` (1) — the widths NCBI's `taxdump_readme.txt` documents and whose every
  field is consumed. A row of a different width means the layout moved, and reading it
  positionally would return a `name_class` that is really a `unique_name`: well-formed to
  every downstream check. Errors name the member and line. `nodes.dmp` is checked as a
  lower bound instead, since its documented width differs between `taxdump` (13 fields) and
  `new_taxdump` (18) and only the first three are read.
- **Taxids must be integers outright.** A non-numeric, partially-numeric (`123abc`) or
  out-of-range value is rejected rather than becoming `0`, being truncated, or clamping to
  `INT64_MAX` — the latter two land inside the valid taxid range where no caller could tell
  them from real data.
- **A reader requires the member it reads**, and only that member (see
  [`read_ncbi_taxdump_deleted`](#list-deleted-taxids-read_ncbi_taxdump_deleted)).

### Read the taxonomy tree (`read_ncbi_taxdump`)

Reads NCBI's taxonomy dump into the tree-table schema. With no argument it downloads the
canonical dump from NCBI and caches the extracted files locally, so subsequent calls are
fast and offline.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded) only when downloading or reading a remote source
- Network access to NCBI (ftp.ncbi.nlm.nih.gov) for the default download

**Parameters:**
- `source` (VARCHAR, positional, optional): where to read the dump from. Omit (or pass NULL) to auto-download the default dump. May be:
  - a **directory** of already-extracted `.dmp` files (read directly; no download/untar)
  - a local or remote (`http://`/`https://`) `taxdump.tar.gz` archive (streamed, gunzipped, and untarred)
- `refresh` (BOOLEAN, optional, default false): re-download the default dump even if a cached copy exists

**Output schema (tree-table + `rank`, consistent with `read_newick`):**
- `node_index` (BIGINT): the tax_id
- `parent_index` (BIGINT, nullable): the parent tax_id; **NULL for the root** (NCBI self-parents taxid 1)
- `name` (VARCHAR): scientific name (`name_class = 'scientific name'`). Every other name class — synonyms, `genbank common name`, authorities — is reached via [`read_ncbi_taxdump_names`](#read-every-name-class-read_ncbi_taxdump_names), not here
- `rank` (VARCHAR): raw NCBI rank string (e.g. `superkingdom`, `domain`, `phylum`, ..., `strain`, `no rank`)
- `is_tip` (BOOLEAN): true iff no node lists this node as its parent

**Behavior:**
- Default source is `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
- The **extracted** `.dmp` files are cached under `${XDG_CACHE_HOME:-$HOME/.cache}/miint/taxonomy/`, so download + gunzip + untar happen once; later reads hit the fast local-directory path. Set `MIINT_TAXONOMY_CACHE_DIR` to override the cache location.
- The tree contains **live nodes only**. Merged-away taxids live in `read_ncbi_taxdump_merged`, deleted ones in `read_ncbi_taxdump_deleted`; neither appears here.
- Emitted faithfully — `unclassified`/`environmental`/`no rank` nodes are included, not filtered.

**Examples:**
```sql
-- Download (once) and read the whole NCBI taxonomy
SELECT COUNT(*) FROM read_ncbi_taxdump();

-- Read a previously extracted dump directory (fully offline)
SELECT node_index, parent_index, name, rank
FROM read_ncbi_taxdump('/data/taxonomy/')
WHERE node_index = 562;

-- Read a specific archive (local or remote)
SELECT * FROM read_ncbi_taxdump('/downloads/taxdump.tar.gz');

-- Force a fresh download of the default dump
SELECT COUNT(*) FROM read_ncbi_taxdump(refresh := true);

-- Join sequence metadata's taxonomy_id against the tree
SELECT m.accession, t.name AS organism, t.rank
FROM read_ncbi(['NC_001416.1', 'NC_000913.3']) m
JOIN read_ncbi_taxdump() t ON t.node_index = m.taxonomy_id;
```

### Remap retired taxids (`read_ncbi_taxdump_merged`)

NCBI periodically merges taxids; the merged (retired) ids never appear in the live tree.
This exposes the retired→current map from `merged.dmp` so you can remap taxids in your own
data **before** joining against the tree — a naive join would silently drop merged taxids.

**Parameters:** same `source` / `refresh` as `read_ncbi_taxdump`.

**Output schema:**
- `old_taxid` (BIGINT): the retired taxid
- `new_taxid` (BIGINT): the current taxid it was merged into

**Example:**
```sql
-- Remap possibly-retired taxids to current ones, then join the tree
WITH remap AS (SELECT * FROM read_ncbi_taxdump_merged())
SELECT d.my_taxid, COALESCE(r.new_taxid, d.my_taxid) AS current_taxid, t.name
FROM my_data d
LEFT JOIN remap r ON r.old_taxid = d.my_taxid
LEFT JOIN read_ncbi_taxdump() t ON t.node_index = COALESCE(r.new_taxid, d.my_taxid);
```

### Read every name class (`read_ncbi_taxdump_names`)

`read_ncbi_taxdump` collapses each taxon to its scientific name, which leaves the rest of
`names.dmp` unreachable. This is the **verbatim** view of that member: one row per name,
no class filtered out. Use it when you need a second name for a taxon — most often
`genbank common name`, which is what someone searching for a taxon is likely to type.

**Parameters:** same `source` / `refresh` as `read_ncbi_taxdump`.

**Output schema:**
- `taxid` (BIGINT): the taxon the name belongs to
- `name` (VARCHAR): the name itself (`name_txt`)
- `unique_name` (VARCHAR): NCBI's disambiguator, supplied when a name string is shared across taxa; **empty string** (not NULL) on most rows
- `name_class` (VARCHAR): which kind of name this is

**Behavior:**
- Rows are emitted in file order, so a taxon's names are grouped but not sorted by class
- A taxon has exactly one `scientific name` but may carry several rows of another class, so aggregate rather than assuming one row per (taxid, class) pair
- Reads **only** `names.dmp` — not the whole dump — so it can be joined against the other readers without each one paying for members it does not use
- The classes present in the current dump, by frequency: `scientific name`, `authority`, `type material`, `synonym`, `includes`, `equivalent name`, `genbank common name`, `common name`, `acronym`, `in-part`, `blast name`. Only ~1% of taxa have a `genbank common name`, so expect NULL for most.

**Examples:**
```sql
-- Scientific name alongside the genbank common name, one row per taxon
SELECT taxid,
       MAX(name) FILTER (name_class = 'scientific name')     AS scientific_name,
       MAX(name) FILTER (name_class = 'genbank common name') AS common_name
FROM read_ncbi_taxdump_names()
WHERE taxid IN (9606, 10090, 562)
GROUP BY taxid;
-- 9606  | Homo sapiens     | human
-- 10090 | Mus musculus     | house mouse
-- 562   | Escherichia coli | NULL      <- no genbank common name

-- Every name NCBI records for one taxon
SELECT name_class, name FROM read_ncbi_taxdump_names() WHERE taxid = 9606;

-- Resolve a free-text name to a taxid across all classes
SELECT DISTINCT taxid FROM read_ncbi_taxdump_names() WHERE name ILIKE 'house mouse';
```

### List deleted taxids (`read_ncbi_taxdump_deleted`)

NCBI deletes taxids outright as well as merging them. A deleted taxid appears in neither
the live tree nor `read_ncbi_taxdump_merged`, so without this member it is
indistinguishable from a taxid that never existed — and a stored taxid that stops
resolving cannot be told apart from bad input. This exposes `delnodes.dmp`.

**Parameters:** same `source` / `refresh` as `read_ncbi_taxdump`.

**Output schema:**
- `taxid` (BIGINT): the deleted taxid

There is no replacement column: unlike a merge, a deletion names no successor.

**Behavior:**
- A dump carrying no `delnodes.dmp` is an **error**, not zero rows. "Nothing was deleted" and "the extraction was incomplete" must not look alike to a caller that treats the result as authoritative — an empty answer read as truth would un-retire every taxid NCBI has deleted. Each reader requires only the member it reads, so a dump without `delnodes.dmp` is still perfectly usable via `read_ncbi_taxdump` and `read_ncbi_taxdump_names`.
- Deletions substantially outnumber merges (~778k vs ~100k in the current dump)
- The live, merged, and deleted id sets are mutually disjoint — a taxon is exactly one of the three

**Examples:**
```sql
-- Classify stored taxids as live, merged away, or deleted
SELECT d.my_taxid,
       CASE WHEN t.node_index IS NOT NULL THEN 'live'
            WHEN m.old_taxid  IS NOT NULL THEN 'merged'
            WHEN x.taxid      IS NOT NULL THEN 'deleted'
            ELSE 'unknown' END AS status,
       m.new_taxid AS replaced_by
FROM my_data d
LEFT JOIN read_ncbi_taxdump()         t ON t.node_index = d.my_taxid
LEFT JOIN read_ncbi_taxdump_merged()  m ON m.old_taxid  = d.my_taxid
LEFT JOIN read_ncbi_taxdump_deleted() x ON x.taxid      = d.my_taxid;

-- How many taxids has NCBI deleted?
SELECT COUNT(*) FROM read_ncbi_taxdump_deleted();
```

### Look up lineages online (`read_ncbi_lineage`)

Resolves a small set of taxids to full lineages via NCBI E-utilities (`efetch`,
`db=taxonomy`). This is the **online** path — for large taxid sets or repeated use, prefer
the cached `read_ncbi_taxdump` + `taxonomy_lineage`, which emit the identical schema.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded)
- Network access to NCBI (eutils.ncbi.nlm.nih.gov)

**Parameters:**
- `taxid` (BIGINT or BIGINT[]): the taxon(s) to resolve
- `api_key` (VARCHAR, optional): NCBI API key for higher rate limits (10 req/s vs 3 req/s)
- `batch_size` (BIGINT, optional, default 500, max 10000): taxids per `epost`+`efetch` round-trip

**Output schema (the shared lineage schema, identical to `taxonomy_lineage`):**
- `taxid` (BIGINT): the query taxon
- `name` (VARCHAR): its scientific name
- `rank` (VARCHAR): its rank
- `domain`, `phylum`, `class`, `"order"`, `family`, `genus`, `species`, `strain` (VARCHAR): the rank-collapsed lineage; NULL where that rank is absent. `class` and `order` are SQL reserved words — quote `"order"` (and, when needed, `"class"`) in queries.
- `lineage` (VARCHAR): formatted `d__…;p__…;c__…;o__…;f__…;g__…;s__…`, with `;t__<strain>` appended when a strain is present

**Behavior:**
- Both the legacy `superkingdom` and the newer `domain` rank collapse to the `domain` column (NCBI is mid-transition)
- Taxids silently omitted by NCBI are reported via `miint_warnings()`, not raised; a taxid that NCBI has since merged is recognised as resolved (via the response's `AkaTaxIds`) rather than reported as missing
- Rate-limited and retried on transient failures with exponential backoff

**Examples:**
```sql
-- One taxid
SELECT * FROM read_ncbi_lineage(562);

-- Several taxids
SELECT taxid, name, genus, species, lineage
FROM read_ncbi_lineage([562, 9606]);

-- Higher rate limit with an API key
SELECT taxid, lineage FROM read_ncbi_lineage([562, 9606], api_key := 'your_ncbi_api_key');
```

### Collapse the tree to lineages (`taxonomy_lineage`)

The **offline** counterpart to `read_ncbi_lineage`: walks a taxdump tree (read internally
via `read_ncbi_taxdump`) from each taxon up to the root and rank-collapses the path into
the **same schema** `read_ncbi_lineage` emits. Because the schema is identical, offline and
online are drop-in interchangeable.

**Parameters:**
- `taxids` (BIGINT[], optional, default NULL): the taxa to resolve. NULL collapses **every** live node (see the memory note below).
- `source` (optional, default NULL): passed straight to `read_ncbi_taxdump` — NULL auto-downloads/uses the cache; a directory or `.tar.gz` reads that instead
- `refresh` (BOOLEAN, optional, default false): passed to `read_ncbi_taxdump`

**Output schema:** identical to [`read_ncbi_lineage`](#look-up-lineages-online-read_ncbi_lineage).

**Behavior:**
- Rank→column mapping (the one [`lineage_rank`](#map-an-ncbi-rank-to-its-column-lineage_rank) gives, shared with `read_ncbi_lineage`): `superkingdom` **or** `domain` → `domain`; `phylum`/`class`/`order`/`family`/`genus`/`species`/`strain` → their like-named columns; all other ranks (`clade`, `subspecies`, `no rank`, …) are dropped. Absent ranks are NULL in their column but rendered empty (e.g. `p__`) in the formatted `lineage`.
- A requested taxid that is not a live node produces no row (same omission semantics as `read_ncbi_lineage`; remap first with `read_ncbi_taxdump_merged` if a taxid may be retired).
- The taxdump is parsed once per call, regardless of how many taxids are requested.

**Examples:**
```sql
-- Same call shape and schema as read_ncbi_lineage, but offline over the cached dump
SELECT * FROM taxonomy_lineage(taxids := [562, 9606]);

-- Read from an extracted dump directory instead of the default cache
SELECT taxid, lineage
FROM taxonomy_lineage(taxids := [562], source := '/data/taxonomy/');

-- Attach lineages to sequence metadata by taxonomy_id
SELECT m.accession, l.lineage
FROM read_ncbi(['NC_001416.1', 'NC_000913.3']) m
JOIN taxonomy_lineage(source := '/data/taxonomy/') l ON l.taxid = m.taxonomy_id;
```

**Notes:**
- **Memory:** calling `taxonomy_lineage()` with `taxids := NULL` over the full NCBI tree materializes a lineage for every one of ~2.9M nodes and is memory-heavy (tens of GB peak). For typical use, always pass `taxids := [...]` — resolving thousands of taxids is fast and light. The unfiltered form is intended for whole-taxonomy exports on machines with ample RAM.

### Parse lineage strings (`parse_lineage`)

Turns a prefixed lineage string, as written by GTDB, Greengenes2, tax2tree and
`taxonomy_lineage` itself, into the eight rank columns of the
[shared lineage schema](#look-up-lineages-online-read_ncbi_lineage). Use it to turn a
`(feature_id, lineage)` taxonomy into those columns. It works on any prefixed lineage, not
only NCBI ones.

**Parameters:**
- `lineage` (VARCHAR): for example `d__Bacteria; p__Bacillota; c__Bacilli` or `d__Bacteria;p__Bacillota;c__Bacilli`
- `domain_letter` (VARCHAR, named only, default `'d'`): the prefix letter of the domain, `'d'` or `'k'`. Greengenes 13_8 writes `k__Bacteria`; pass `domain_letter := 'k'` to read it. It must be a constant, so one call cannot mix `d__` and `k__` rows.
- `allow_subset` (BOOLEAN, named only, default `false`): accept a lineage that starts below the domain, such as `c__Bacilli; o__Lactobacillales` or `p__Bacillota`. The internal node labels of a tree decorated by tax2tree are like this: each starts at the rank where decoration began. The ranks above the first field are NULL, and the result type is the same. It must be a constant.

**Returns:** `STRUCT(domain, phylum, class, "order", family, genus, species, strain)`, all
VARCHAR: the shared schema's rank columns, in its order. A NULL `lineage` returns NULL.
DuckDB 2.0 does not accept `(expr).*`, so use `unnest(parse_lineage(lineage))` to turn
the struct into columns.

**Rules:**
- Fields are separated by `;`, and whitespace around a field is ignored, so `'; '` and `';'` both work.
- Each field is `<letter>__<name>`, and its letter must match its position: `d p c o f g s t` (domain through strain). A skipped rank therefore errors rather than shifting the names below it into the wrong columns.
- With `allow_subset := true`, the first field's letter picks the rank the lineage starts at, and each later field must carry the next rank's letter. `c__A; o__B` gives class `A` and order `B`; `c__A; f__B` (order skipped) and `o__B; c__A` (reversed) are errors.
- An empty name (`p__`) is NULL, and a lineage with fewer than eight fields leaves the deeper ranks NULL. `d__A;p__;c__;o__;f__;g__;s__` gives domain `A` and seven NULLs.
- Lineages without rank prefixes are not supported: raw SILVA, and the QIIME releases of SILVA up to 132, whose `D_0__`, `D_1__`, ... prefixes name a depth rather than a rank.
- A `lineage` written by `taxonomy_lineage` or `read_ncbi_lineage` parses back into the columns it was written from, provided no name contains `;`, starts or ends with whitespace, or itself starts with a rank prefix such as `p__`. NCBI has a few names containing `;`, such as the species `Rhodococcus sp. SJ-2;3`; each splits into an extra field that fails its prefix check, so `parse_lineage` errors on those lineages.
- `parse_lineage` does not reject gaps: `d__A; p__; c__C` returns phylum NULL and class `C`. Consumers that need gap-free lineages check that themselves.

**Errors** quote the lineage and name the field:
- an empty lineage, or an empty field (a trailing or doubled `;`)
- a prefix letter out of position, or `k__` without `domain_letter := 'k'` (and `d__` with it)
- a lineage that starts below the domain without `allow_subset := true`
- a field after strain (in a lineage from the domain, a ninth field)
- a name with whitespace after its prefix (`d__ Bacteria`) or a second prefix (`p__p__Bacillota`, usually a prefix that should have been stripped)

**Examples:**
```sql
-- A QIIME 2 taxonomy.tsv (header: Feature ID, Taxon, Confidence). The classifier writes
-- 'Unassigned', which has no prefix, for features it could not classify.
SELECT "Feature ID" AS feature_id, unnest(parse_lineage(Taxon))
FROM read_csv('taxonomy.tsv', delim = '\t', header = true)
WHERE Taxon <> 'Unassigned';

-- Greengenes 13_8, which writes k__ for the domain
SELECT feature_id, unnest(parse_lineage(lineage, domain_letter := 'k'))
FROM gg_13_8_taxonomy;

-- The internal node labels of a tree decorated by tax2tree, such as
-- 'f__Lachnospiraceae; g__Blautia'. (When tax2tree keeps support values it writes
-- '95:f__Lachnospiraceae; g__Blautia'; strip the support first.) read_newick gives an
-- unlabeled node the name ''.
SELECT node_index, unnest(parse_lineage(name, allow_subset := true))
FROM read_newick('decorated.nwk')
WHERE NOT is_tip AND name <> '';

-- An error stops at the first bad row. TRY lists them all.
SELECT feature_id, lineage
FROM raw_taxonomy
WHERE lineage IS NOT NULL AND TRY(parse_lineage(lineage)) IS NULL;
```

### Map an NCBI rank to its column (`lineage_rank`)

Names the [shared lineage schema](#look-up-lineages-online-read_ncbi_lineage)'s column that
an NCBI rank collapses into, the mapping `read_ncbi_lineage` and `taxonomy_lineage` use.
Older taxdumps write the rank `superkingdom` and newer ones `domain`, so a lookup that
filters on the raw rank finds nothing in a taxdump written with the other one;
`lineage_rank` gives `domain` for both.

**Parameters:**
- `ncbi_rank` (VARCHAR): an NCBI rank, as in `read_ncbi_taxdump`'s `rank` column.

**Returns:** VARCHAR, one of `domain`, `phylum`, `class`, `order`, `family`, `genus`,
`species` or `strain`:
- `superkingdom` and `domain` give `domain`; `phylum`, `class`, `order`, `family`, `genus`, `species` and `strain` give their own name.
- Every other rank (`clade`, `kingdom`, `subspecies`, `no rank`, ...) gives NULL: the schema drops it, as the lineage functions do. So does a NULL `ncbi_rank`.
- Ranks are matched exactly as NCBI writes them, in lower case: `Phylum` and `' phylum'` give NULL.

**Examples:**
```sql
-- The domains of a taxdump, whether it ranks them superkingdom or domain.
SELECT node_index, name, rank
FROM read_ncbi_taxdump()
WHERE lineage_rank(rank) = 'domain';

-- How many nodes of each NCBI rank fall into each schema column.
SELECT lineage_rank(rank) AS column_name, rank, count(*) AS nodes
FROM read_ncbi_taxdump()
GROUP BY ALL
ORDER BY column_name NULLS LAST, rank;
```

### Split GTDB labels (`gtdb_label`, `gtdb_species_label`)

When GTDB splits a taxon, the part that holds its type keeps the name and the other parts
carry a suffix: `Bacillota_A`, `Pseudomonas_E`, and at species `Pseudomonas_E fluorescens_B`.
These functions split a name into its base name, which the parts share, and its suffix, so a
split taxon, polyphyletic ones included, can be found by grouping on the base. [`parse_lineage`](#parse-lineage-strings-parse_lineage) keeps names
whole, as the shared schema does; split them afterwards. The functions only split: any
grouping key is the caller's to build from the parts.

**Parameters:**
- `name` (VARCHAR): one name, without its rank prefix, as `parse_lineage` returns it (`Bacillota_A`, not `p__Bacillota_A`).

**Returns:** all fields VARCHAR.
- `gtdb_label(name)`: `STRUCT(base, suffix)`, for the ranks above species.
- `gtdb_species_label(name)`: `STRUCT(genus_base, genus_suffix, epithet, epithet_suffix)`. `epithet` is everything after the genus label and its space, without the epithet's suffix: the specific epithet, or for a name such as `Foo bar subsp. baz`, all of `bar subsp. baz`.

A suffix is returned without its underscore (`Bacillota_A` gives `Bacillota` and `A`), and a
name without one has a NULL suffix. A name that does not fit the rules below, or a NULL
name, gives NULL rather than an error, and so do its fields.

**Rules:**
- A base name is ASCII letters, digits and `-`: `CAG-41`, `UBA1234`.
- A suffix is `_` followed by runs of capital letters and digits joined by single underscores: `_A`, `_A1`, or `_Y_123`, which is one suffix, `Y_123`. `Foo_A_`, `Foo__A` and `Foo_a` do not fit.
- A species is the genus label, one space, then the epithet, each with its own optional suffix. The epithet may hold single spaces and dots, as in `bar subsp. baz`, and its suffix may only sit at its very end: `bar subsp. baz_A` fits, `bar_A subsp. baz` does not.

**Examples:**
```sql
-- A GTDB taxonomy file (bac120_taxonomy.tsv: accession, then the lineage), in the
-- shared schema's columns.
CREATE TABLE gtdb AS
SELECT column0 AS accession, unnest(parse_lineage(column1))
FROM read_csv('bac120_taxonomy.tsv', delim = '\t', header = false);

-- Phyla that GTDB splits: a base name some of whose labels carry a suffix.
SELECT gtdb_label(phylum).base AS base, list(DISTINCT phylum ORDER BY phylum) AS labels
FROM gtdb
GROUP BY base
HAVING bool_or(gtdb_label(phylum).suffix IS NOT NULL);

-- Species that share a genus base and an epithet but differ by suffix. The key, here
-- the genus base and the epithet, is the caller's choice.
SELECT s.genus_base, s.epithet, list(DISTINCT species ORDER BY species) AS labels
FROM (SELECT species, gtdb_species_label(species) AS s FROM gtdb)
WHERE s IS NOT NULL
GROUP BY ALL
HAVING count(DISTINCT species) > 1;

-- Names that do not fit the rules.
SELECT DISTINCT genus FROM gtdb WHERE genus IS NOT NULL AND gtdb_label(genus) IS NULL;
```

### Read GTDB's genome metadata (`read_gtdb_metadata`, `read_gtdb_qc_failed`)

With each release, GTDB (the Genome Taxonomy Database) publishes a table of every genome
it places: its GTDB and NCBI taxonomies, its quality estimates, and whether it is type
material. `read_gtdb_metadata` reads it, from GTDB's `bac120` (bacterial) and `ar53`
(archaeal) files. `read_gtdb_qc_failed` reads the list of genomes that failed GTDB's
quality checks and were left out.

**Requirements:**
- Requires the `httpfs` extension (automatically loaded) only to download a release or read a URL
- Network access to GTDB (data.gtdb.ecogenomic.org) for a release not yet cached

**Parameters** (named only; pass exactly one of `release` and `source`):
- `release` (VARCHAR): a GTDB release, such as `'232'` or `'214.1'`. A bare number is the release's first, `.0`. Its files are downloaded into a cache when the query runs, and read from there after. Releases from 214.1 serve the files these functions read. Earlier ones ship `.tar.gz` archives: extract the `.tsv` and pass it as `source`.
- `source` (VARCHAR): what to read, as it is, never cached:
  - a `.tsv` or `.tsv.gz` file, local or `http(s)://`. Gzip is told by the file's bytes, not its name.
  - a directory holding one release's files, such as a release's cache directory. For `read_gtdb_metadata`, at most one `bac120_metadata_r<N>.tsv[.gz]` and one `ar53_metadata_r<N>.tsv[.gz]` (`ar122_` before release 207), at least one of the two, of the same release. For `read_gtdb_qc_failed`, exactly one `qc_failed_r<N>.tsv[.gz]`. Other files are ignored.
- `refresh` (BOOLEAN, default false): with `release`, download its files again though they are cached. With `source` it is an error.

**Output:**
- `read_gtdb_metadata`: `domain` (VARCHAR), then every column of the file, named and ordered as in its header, all VARCHAR, with GTDB's missing-value marker `none` read as NULL.
  - `domain` is `'Bacteria'` for a `bac120_` file and `'Archaea'` for an `ar53_` or `ar122_` one, matching `parse_lineage`'s domain and GTDB's `d__` names, so the two join.
  - The columns change between releases: release 232 has 113, and 214.1 has 110, without CheckM2's.
  - Cast what you need, for example `TRY_CAST(checkm2_completeness AS DOUBLE)`. GTDB writes booleans as `t` and `f`.
- `read_gtdb_qc_failed`: every column of the file, named as GTDB names them (`"Accession"`, `"NCBI species"`, `"CM1 completeness (%)"`; quote them in SQL), all VARCHAR, values as written (`-` and `s__` included).
- In both, an empty field is NULL, and the order of the rows is not promised.

**Linking an NCBI species to its type strain's genome:** `gtdb_type_designation_ncbi_taxa`
says whether a genome is type material in NCBI's taxonomy (`type strain of species`,
`type strain of heterotypic synonym`, `not type material`, ...), and `ncbi_species_taxid` is
the NCBI species it belongs to. Together they take a species' NCBI taxid to the genome of
its type strain, and so to GTDB's placement of that strain, without matching names. The
examples below show the join.

**Behavior:**
- **The cache.** A release's files are kept, compressed as served, at `<cache>/<release>/<file name>`, such as `~/.cache/miint/gtdb/232.0/bac120_metadata_r232.tsv.gz`. `<cache>` is the `miint_gtdb_cache_dir` setting, or `${XDG_CACHE_HOME:-~/.cache}/miint/gtdb` when that is empty (the default). Every file sits in the release's directory, `qc_failed` included, so the directory is a `source` for either function. Only the files a call reads are downloaded: `read_gtdb_qc_failed` does not fetch the metadata.
- **Downloads.** A file is downloaded to a file of its own beside its final name, which is renamed into place once complete. So a file at its final name is complete, and it is trusted: nothing checks it against GTDB again until `refresh := true`. A failed or interrupted download removes its file. A killed process leaves `<name>.<id>.part`, which is never read and can be deleted.
- **Binding reads only header lines**, which name the columns: from the cache, or for a file not yet cached, from GTDB by a request stopped once the line has arrived. So `DESCRIBE` or `EXPLAIN` of a release not yet cached makes one small request per file and downloads nothing. The download happens when the query runs.
- **Each file's header is read again when the query runs.** A file whose header changed after binding, rewritten under a prepared statement, say, fails the query: prepare it again.
- **Rows are read by DuckDB's CSV reader** in GTDB's dialect: tab-separated, with a value holding `"` quoted and the quote doubled. A row with more or fewer fields than the header is an error.
- **Retries.** A download that fails with HTTP 5xx, 429 or 408 is tried again after 1, 2 and 4 s, and DuckDB's own HTTP retries (`http_retries`) run within each attempt. Ctrl-C and `max_execution_time` end a download or a wait ([Requests: waits and where they go](#requests-waits-and-where-they-go)).

**Settings:**
- `miint_gtdb_url` (default `https://data.gtdb.ecogenomic.org/releases`): where releases are downloaded from, as `<url>/release232/232.0/bac120_metadata_r232.tsv.gz` and `<url>/release232/232.0/auxillary_files/qc_failed_r232.tsv` (GTDB's spelling). A mirror with GTDB's layout works.
- `miint_gtdb_cache_dir` (default empty): where the cache is kept; see above.

**Errors:**
- both or neither of `release` and `source`, a `release` that is not digits with at most one dot, an empty `source`, or `refresh` with `source`
- a `source` that does not exist, or a directory that does not hold one release's files
- a metadata file named for neither domain
- `bac120` and `ar53` files whose headers differ (one release's never do)
- a header with an empty column name, a name that repeats another ignoring case, a name holding `"`, or, for `read_gtdb_metadata`, a column named `domain`
- a file with no header line
- a release GTDB does not serve these files for, which fails with HTTP 404 and the hint to pass the extracted `.tsv` as `source`

**Examples:**
```sql
-- Release 232's genomes, downloaded once (about 300 MB) and read from the cache after
SELECT domain, count(*) FROM read_gtdb_metadata(release := '232') GROUP BY domain;

-- GTDB's representative genomes and their GTDB taxonomy, as the shared lineage schema
SELECT accession, unnest(parse_lineage(gtdb_taxonomy))
FROM read_gtdb_metadata(release := '232')
WHERE gtdb_representative = 't';

-- The genome of each NCBI species' type strain, for records that carry a species taxid
SELECT r.record_id, m.accession, m.gtdb_taxonomy
FROM records r
JOIN read_gtdb_metadata(release := '232') m
  ON m.ncbi_species_taxid = r.species_taxid::VARCHAR
 AND m.gtdb_type_designation_ncbi_taxa = 'type strain of species';

-- Genomes that failed GTDB's CheckM2 completeness check, and why NCBI excludes them
SELECT "Accession", "GTDB taxonomy", "NCBI exclude from RefSeq"
FROM read_gtdb_qc_failed(release := '232')
WHERE "Failed CM2 completeness" = '1';

-- A release before 214.1: extract the file from GTDB's archive and read it
SELECT count(*) FROM read_gtdb_metadata(source := 'bac120_metadata_r207.tsv');

-- The cache directory is a source too, which reads it without going to GTDB
SELECT count(*) FROM read_gtdb_qc_failed(source := '~/.cache/miint/gtdb/232.0');
```
