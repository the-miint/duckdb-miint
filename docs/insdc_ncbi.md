# Working with NCBI

MIINT provides integration with NCBI's GenBank and RefSeq.

## Table of contents

- [Read metadata](#read-metadata) - Fetch record metadata from GenBank and RefSeq by accession.
- [Read sequences](#read-sequences) - Fetch sequence data from GenBank and RefSeq by accession.
- [Read annotations](#read-annotations) - Fetch annotation data from GenBank and RefSeq by accession.

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
- `stop_position` (INTEGER): 1-based end position
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
