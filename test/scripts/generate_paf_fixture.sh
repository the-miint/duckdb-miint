#!/usr/bin/env bash
# Regenerate data/paf/: a reference slice, reads cut from it at known coordinates, and the
# PAF minimap2 writes for them. Read names carry the expected alignment, so tests check
# read_paf against how the reads were built rather than against minimap2's own output.
#
# Usage: bash test/scripts/generate_paf_fixture.sh
# Requires: ./build/release/duckdb (miint built) and minimap2 on PATH (or MINIMAP2=...).
#
# The reference is E. coli K-12 W3110 [100000, 120000), named k12_100000_120000. Read
# coordinates are 0-based half-open on that reference record:
#   fwd_<start>_<stop>     the reference slice, forward strand
#   rev_<start>_<stop>     the reverse complement of the slice
#   ovhfwd_<start>_<stop>  the slice followed by 1500 bp from elsewhere in the genome;
#                          aligns over query [0, stop - start)
#   ovhrev_<start>_<stop>  the reverse complement of (slice followed by 1500 bp from
#                          elsewhere); aligns over query [1500, 1500 + stop - start)
#   unmapped               3000 bp from outside the reference
# Both overhang slices end at the reference's end, so an alignment cannot extend into the
# foreign sequence by chance matches.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
DUCKDB=${DUCKDB:-$REPO/build/release/duckdb}
MINIMAP2=${MINIMAP2:-minimap2}
DEST=$REPO/data/paf

# data/paf/README.md describes what this version writes; another may align differently.
FOUND=$("$MINIMAP2" --version)
if [ "$FOUND" != "2.31-r1302" ]; then
    echo "ERROR: expected minimap2 2.31-r1302, found $FOUND" >&2
    exit 1
fi

mkdir -p "$DEST"
cd "$REPO"
"$DUCKDB" -batch <<SQL
CREATE TEMP TABLE k12 AS
SELECT sequence1 AS seq FROM read_fastx('data/sylph/tiny_refs/e.coli-K12.fasta.gz');
CREATE TEMP MACRO ref_slice(s, e) AS (SELECT substr(seq, 100000 + s + 1, e - s) FROM k12);
CREATE TEMP MACRO genome_slice(s, e) AS (SELECT substr(seq, s + 1, e - s) FROM k12);

COPY (SELECT 'k12_100000_120000' AS read_id, ref_slice(0, 20000) AS sequence1)
TO 'data/paf/k12_ref.fa' (FORMAT FASTA);

COPY (
    SELECT * FROM (VALUES
        ('fwd_1000_3000', ref_slice(1000, 3000)),
        ('fwd_5000_9000', ref_slice(5000, 9000)),
        ('rev_10000_13000', sequence_dna_reverse_complement(ref_slice(10000, 13000))),
        ('rev_14000_18000', sequence_dna_reverse_complement(ref_slice(14000, 18000))),
        ('ovhfwd_16000_20000', ref_slice(16000, 20000) || genome_slice(3500000, 3501500)),
        ('ovhrev_17000_20000',
         sequence_dna_reverse_complement(ref_slice(17000, 20000) || genome_slice(3510000, 3511500))),
        ('unmapped', genome_slice(3600000, 3603000))
    ) t(read_id, sequence1)
) TO 'data/paf/k12_reads.fa' (FORMAT FASTA);
SQL

"$MINIMAP2" -x map-hifi -c --paf-no-hit -t 1 data/paf/k12_ref.fa data/paf/k12_reads.fa \
    > data/paf/k12_minimap2.paf
gzip -n -c data/paf/k12_minimap2.paf > data/paf/k12_minimap2.paf.gz
