#!/usr/bin/env bash
# Regenerate data/myloasm/ecoli_sim/: real myloasm output over simulated reads whose
# names encode where each read came from, so tests can check the myloasm readers
# against construction ground truth rather than against myloasm's own claims.
#
# Usage: bash test/scripts/generate_myloasm_fixture.sh WORK_DIR
#   WORK_DIR must not exist; the reads and both full myloasm runs are left there.
#
# Requires: ./build/release/duckdb (miint built), and myloasm in the conda env named
# by MYLOASM_ENV (default myloasm-070; bioconda myloasm=0.7.0).
#
# Every read is cut from E. coli K-12 W3110 (data/sylph/tiny_refs/e.coli-K12.fasta.gz).
# Coordinates in read names are 0-based half-open on that genome. Each name is
#   s<index>_<kind>_<fields...>_<f|r>
# where the trailing f/r says whether the read is the forward strand or the reverse
# complement of the sequence described by the fields:
#   lina  <start>_<stop>                 strain A of the linear region [1000000, 1200000)
#   linb  <start>_<stop>                 strain B: the same region with a substitution at
#                                        every offset where (offset mod 150) = 75
#   circ  <start>_<length>               a read starting at <start> in the circular
#                                        replicon [2000000, 2040000), wrapping past its end
#   low   <start>_<stop>                 the low-coverage region [3000000, 3060000)
#   circb <start>_<length>               like circ, but from a variant of the replicon with
#                                        a substitution at every offset where
#                                        (offset mod 1000) = 7
#   iso   <start>_<stop>                 three reads from [3300000, 3400000) far enough
#                                        apart that none overlaps another
#   chim  <a_start>_<a_stop>_<c_start>_<c_length>
#                                        strain-A sequence [a_start, a_stop) joined to the
#                                        circular-replicon read <c_start>_<c_length>
#   ovh   <start>_<stop>_<tail_start>    strain-A sequence [start, stop) followed by 3000 bp
#                                        from [tail_start, tail_start + 3000), a region no
#                                        other read covers
# A reverse (r) read is the reverse complement of the whole forward construction, so for
# chim and ovh the junction sits at (length of the second part) from the read's start.
set -euo pipefail

REPO=$(cd "$(dirname "$0")/../.." && pwd)
DUCKDB=${DUCKDB:-$REPO/build/release/duckdb}
MYLOASM_ENV=${MYLOASM_ENV:-myloasm-070}
WORK=${1:?usage: generate_myloasm_fixture.sh WORK_DIR}
DEST=$REPO/data/myloasm/ecoli_sim

# data/myloasm/README.md describes what this version writes; another may assemble differently.
FOUND=$(conda run -n "$MYLOASM_ENV" myloasm --version)
if [ "$FOUND" != "myloasm 0.7.0" ]; then
    echo "ERROR: expected myloasm 0.7.0 in $MYLOASM_ENV, found $FOUND" >&2
    exit 1
fi

if [ -e "$WORK" ]; then
    echo "ERROR: $WORK already exists; pass a fresh directory" >&2
    exit 1
fi
mkdir -p "$WORK"
WORK=$(cd "$WORK" && pwd)

cd "$REPO"
"$DUCKDB" -batch <<SQL
CREATE TEMP TABLE k12 AS
SELECT sequence1 AS seq FROM read_fastx('data/sylph/tiny_refs/e.coli-K12.fasta.gz');

-- Strain B of the linear region: one substitution every 150 bp (99.33% identity), enough
-- for myloasm to assemble it as an alternate of strain A.
CREATE TEMP TABLE linb AS
SELECT string_agg(CASE WHEN i % 150 = 75
                       THEN CASE b WHEN 'A' THEN 'C' WHEN 'C' THEN 'G'
                                   WHEN 'G' THEN 'T' WHEN 'T' THEN 'A' END
                       ELSE b END, '' ORDER BY i) AS seq
FROM (SELECT i, substr(k12.seq, 1000000 + i + 1, 1) AS b FROM k12, range(200000) t(i));

-- The circular replicon, doubled so a read can wrap past its end.
CREATE TEMP TABLE circ AS
SELECT substr(seq, 2000000 + 1, 40000) || substr(seq, 2000000 + 1, 40000) AS seq2 FROM k12;

-- A variant of the circular replicon, one substitution every 1000 bp (99.9% identity):
-- close enough that myloasm reports one of the two circles as a duplicate.
CREATE TEMP TABLE circb AS
SELECT s || s AS seq2
FROM (SELECT string_agg(CASE WHEN i % 1000 = 7
                             THEN CASE b WHEN 'A' THEN 'C' WHEN 'C' THEN 'G'
                                         WHEN 'G' THEN 'T' WHEN 'T' THEN 'A' END
                             ELSE b END, '' ORDER BY i) AS s
      FROM (SELECT i, substr(k12.seq, 2000000 + i + 1, 1) AS b FROM k12, range(40000) t(i)));

-- Stratified start positions with a jitter; lengths 8-16 kb. Both come from Knuth's
-- multiplicative hash of the read index: deterministic integer arithmetic, but without
-- the regular pattern a plain (index * prime) mod n leaves, which myloasm reads as
-- coverage structure and breaks contigs on.
-- (kind, global index offset, read count, region start, region length)
CREATE TEMP TABLE plan AS
SELECT * FROM (VALUES
    ('lina',    0, 500, 1000000, 200000),
    ('linb',  500, 333, 1000000, 200000),
    ('circ',  833,  67, 2000000,  40000),
    ('low',   900,  13, 3000000,  60000),
    ('iso',   953,   3, 3300000, 100000),
    ('circb', 956,  50, 2000000,  40000)
) t(kind, base, n, rstart, rlen);

-- Linear regions keep every read inside the region; circular starts cover the whole
-- circle because the doubled sequence lets a read run past the end.
CREATE TEMP TABLE simple_reads AS
SELECT kind, idx, rstart, rlen, len,
       CASE WHEN kind IN ('circ', 'circb') THEN (stratum + jitter) % rlen
            ELSE least(stratum + jitter, rlen - len) END AS rel
FROM (SELECT p.kind, p.base + j AS idx, p.rstart, p.rlen,
             8000 + (((p.base + j) * 2654435761) % 4294967296) % 8001 AS len,
             CASE WHEN p.kind IN ('circ', 'circb') THEN (j * p.rlen) // p.n
                  ELSE (j * (p.rlen - 16000)) // p.n END AS stratum,
             (((p.base + j) * 2246822519) % 4294967296) % 400 AS jitter
      FROM plan p, range(500) t(j)
      WHERE j < p.n);

CREATE TEMP TABLE reads_fwd AS
SELECT idx, kind,
       CASE
           WHEN kind IN ('circ', 'circb') THEN printf('%d_%d', rstart + rel, len)
           ELSE printf('%d_%d', rstart + rel, rstart + rel + len)
       END AS fields,
       CASE kind
           WHEN 'linb' THEN (SELECT substr(seq, rel + 1, len) FROM linb)
           WHEN 'circ' THEN (SELECT substr(seq2, rel + 1, len) FROM circ)
           WHEN 'circb' THEN (SELECT substr(seq2, rel + 1, len) FROM circb)
           ELSE (SELECT substr(seq, rstart + rel + 1, len) FROM k12)
       END AS seq
FROM simple_reads
UNION ALL
-- Chimeras: strain-A sequence joined to a circular-replicon read.
SELECT 913 + j AS idx, 'chim' AS kind,
       printf('%d_%d_%d_%d', a_start, a_start + a_len, 2000000 + c_rel, c_len) AS fields,
       (SELECT substr(seq, a_start + 1, a_len) FROM k12) ||
       (SELECT substr(seq2, c_rel + 1, c_len) FROM circ) AS seq
FROM (SELECT j,
             1000000 + 5000 + j * 9000 AS a_start, 6000 + (j * 7919) % 2001 AS a_len,
             (j * 1999) % 40000 AS c_rel, 6000 + (j * 104729) % 2001 AS c_len
      FROM range(20) t(j))
UNION ALL
-- Overhang reads: strain-A sequence followed by sequence from an uncovered region.
SELECT 933 + j AS idx, 'ovh' AS kind,
       printf('%d_%d_%d', g_start, g_start + 9000, t_start) AS fields,
       (SELECT substr(seq, g_start + 1, 9000) || substr(seq, t_start + 1, 3000) FROM k12) AS seq
FROM (SELECT j, 1000000 + 3000 + j * 9500 AS g_start, 4200000 + j * 5000 AS t_start
      FROM range(20) t(j));

-- Odd indices are reverse-complemented.
CREATE TEMP TABLE reads AS
SELECT printf('s%05d_%s_%s_%s', idx, kind, fields, CASE WHEN idx % 2 = 1 THEN 'r' ELSE 'f' END)
           AS read_id,
       CASE WHEN idx % 2 = 1 THEN sequence_dna_reverse_complement(seq) ELSE seq END AS sequence1
FROM reads_fwd;

-- Fail loudly rather than emit a fixture built from a bad slice.
SELECT CASE WHEN count(*) <> 1006 THEN error('expected 1006 reads, got ' || count(*)) END
FROM reads;
SELECT CASE WHEN count(*) > 0 THEN error('a read contains a non-ACGT base or is empty') END
FROM reads WHERE sequence1 IS NULL OR NOT regexp_full_match(sequence1, '[ACGT]+');
SELECT CASE WHEN (SELECT length(seq) FROM linb) <> 200000 THEN error('strain B is not 200 kb') END;

COPY (SELECT read_id, sequence1 FROM reads ORDER BY read_id) TO '$WORK/reads.fa' (FORMAT FASTA);
SELECT count(*) AS n_reads, sum(length(sequence1)) AS n_bases FROM reads;
SQL

# One thread: with more, myloasm writes the PAF and breakpoint logs in thread order.
# Relative paths keep this machine's directories out of the committed log's COMMAND line.
for run in run1 run2; do
    (cd "$WORK" && conda run -n "$MYLOASM_ENV" myloasm reads.fa -o "$run" --hifi -t 1)
done

# The fixture is committed, so a second run must reproduce it exactly. Logs carry
# timestamps and are excluded.
CONSUMED_PLAIN="assembly_primary.fa alternate_assemblies/assembly_alternate.fa
alternate_assemblies/duplicated_contigs.fa final_contig_graph.gfa"
CONSUMED_GZ="3-mapping/map_to_unitigs.paf.gz 0-cleaning_and_unitigs/read_coverages.txt.gz
0-cleaning_and_unitigs/remap_temp/read_coverages.txt.gz"
for f in $CONSUMED_PLAIN; do
    cmp "$WORK/run1/$f" "$WORK/run2/$f" || { echo "ERROR: $f differs between runs" >&2; exit 1; }
done
for f in $CONSUMED_GZ; do
    cmp <(zcat "$WORK/run1/$f") <(zcat "$WORK/run2/$f") \
        || { echo "ERROR: $f differs between runs" >&2; exit 1; }
done

for f in $CONSUMED_PLAIN $CONSUMED_GZ; do
    mkdir -p "$DEST/$(dirname "$f")"
    cp "$WORK/run1/$f" "$DEST/$f"
done
# myloasm names its log myloasm_<timestamp>.log; a fixed name means regenerating replaces
# the committed log instead of adding a second one beside it. The host name of the
# machine that ran myloasm is redacted.
sed -E 's/(SYSTEM HOST NAME: ).*/\1redacted/' "$WORK"/run1/myloasm*.log > "$DEST/myloasm.log"

du -sh "$DEST"
