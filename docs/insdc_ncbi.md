# Working with NCBI

MIINT provides integration with NCBI's GenBank and RefSeq.

## Table of contents

- [Accessions are literals, not table names](#accessions-are-literals-not-table-names) - Applies to all three readers below; includes the idiom for driving them from accessions held in a table.
- [Read metadata](#read-metadata) - Fetch record metadata from GenBank and RefSeq by accession.
- [Read sequences](#read-sequences) - Fetch sequence data from GenBank and RefSeq by accession.
- [Read annotations](#read-annotations) - Fetch annotation data from GenBank and RefSeq by accession.
- [Taxonomy](#taxonomy) - Obtain and use the NCBI taxonomy: the full tree (bulk/offline), online lineage lookups, and offline lineage collapse.
  - [Read the taxonomy tree (`read_ncbi_taxdump`)](#read-the-taxonomy-tree-read_ncbi_taxdump) - Bulk/offline: download + cache the taxonomy tree.
  - [Remap retired taxids (`read_ncbi_taxdump_merged`)](#remap-retired-taxids-read_ncbi_taxdump_merged) - The retired→current taxid map.
  - [Look up lineages online (`read_ncbi_lineage`)](#look-up-lineages-online-read_ncbi_lineage) - Online/rate-limited: resolve a handful of taxids via E-utilities.
  - [Collapse the tree to lineages (`taxonomy_lineage`)](#collapse-the-tree-to-lineages-taxonomy_lineage) - Offline: rank-collapse a taxdump tree, identical schema to `read_ncbi_lineage`.

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
- `filepath` (VARCHAR, optional): NCBI E-utilities URL when include_filepath=true

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
- `filepath` (VARCHAR, optional): NCBI E-utilities URL when include_filepath=true

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
- `name` (VARCHAR): scientific name (`name_class = 'scientific name'`; synonyms/common names are not emitted here)
- `rank` (VARCHAR): raw NCBI rank string (e.g. `superkingdom`, `domain`, `phylum`, ..., `strain`, `no rank`)
- `is_tip` (BOOLEAN): true iff no node lists this node as its parent

**Behavior:**
- Default source is `https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz`
- The **extracted** `.dmp` files are cached under `${XDG_CACHE_HOME:-$HOME/.cache}/miint/taxonomy/`, so download + gunzip + untar happen once; later reads hit the fast local-directory path. Set `MIINT_TAXONOMY_CACHE_DIR` to override the cache location.
- The tree contains **live nodes only**. Retired taxids live in `read_ncbi_taxdump_merged`; deleted taxids are simply absent.
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
- Rank→column mapping (the single source of truth, mirrored from `read_ncbi_lineage`): `superkingdom` **or** `domain` → `domain`; `phylum`/`class`/`order`/`family`/`genus`/`species`/`strain` → their like-named columns; all other ranks (`clade`, `subspecies`, `no rank`, …) are dropped. Absent ranks are NULL in their column but rendered empty (e.g. `p__`) in the formatted `lineage`.
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
