#!/usr/bin/env python3
"""Generate paired-end FASTQ fixtures for trim_adapters_pe fastp-parity tests.

These fixtures exercise fastp's *overlap-based* PE adapter trimming: R1 and R2
sequence the same insert from opposite ends, so when the insert is shorter than
the read length each mate reads through into adapter. fastp infers the insert
boundary purely from the R1 / reverse-complement(R2) overlap — no adapter
sequence required.

The committed .fq files are the test INPUT. Their *.golden.{r1,r2}.fq siblings
are the FROZEN expected output — native fastp's overlap-trimmed reads, captured
once and set in stone. The trim_adapters_pe parity test compares against them
directly; there is intentionally no automated regeneration.

Frozen oracle provenance — golden produced with **fastp 1.3.3** via:

    fastp -i data/qc/pe_overlap.r1.fq -I data/qc/pe_overlap.r2.fq \\
          -o data/qc/pe_overlap.golden.r1.fq -O data/qc/pe_overlap.golden.r2.fq \\
          --disable_quality_filtering --disable_length_filtering \\
          --disable_trim_poly_g --thread 1 --dont_eval_duplication

Those isolation flags make fastp run ONLY overlap-based PE adapter trimming.

The pe_fallback.* golden additionally engages fastp's sequence-based adapter
trimming (step 9), so it adds explicit adapter sequences:

    fastp -i data/qc/pe_fallback.r1.fq -I data/qc/pe_fallback.r2.fq \\
          -o data/qc/pe_fallback.golden.r1.fq -O data/qc/pe_fallback.golden.r2.fq \\
          --disable_quality_filtering --disable_length_filtering \\
          --disable_trim_poly_g --thread 1 --dont_eval_duplication \\
          --adapter_sequence AGATCGGAAGAGC --adapter_sequence_r2 AGATCGGAAGAGC

Note: the overlap path is byte-for-byte fastp parity; the fallback matches fastp
on clean exact-adapter cases (as here) but uses miint's adapter matcher, which
diverges from fastp on ambiguous mismatch/indel cases (see docs/qc.md).

If a fixture is ever changed, rerun this script and re-run the matching command
by hand to refresh the golden (and update the version noted here if fastp changed).

Run from the repo root:  python3 test/scripts/generate_qc_pe_fixtures.py
"""

import os

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(s: str) -> str:
    return s.translate(_COMP)[::-1]


# Illumina TruSeq read-through adapters (only their constant 5' portion is ever
# seen in a read). The exact tail bytes are irrelevant to overlap analysis — it
# never inspects them — but using real adapters keeps the fixture realistic.
R1_ADAPTER = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
R2_ADAPTER = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

READLEN = 50


def quals(n: int) -> str:
    # Position-dependent Q30..Q39 ramp (ASCII-33). Varied so the parity test
    # catches any off-by-one in quality slicing, but always valid Phred.
    return "".join(chr(33 + 30 + (i % 10)) for i in range(n))


def make_pair(insert: str):
    """Return (r1_seq, r2_seq) for a fragment of length len(insert).

    R1 reads the insert 5'->3' then into R1_ADAPTER; R2 reads the
    reverse-complement strand then into R2_ADAPTER. Each is truncated to READLEN
    (simulating a fixed-length sequencer)."""
    r1 = (insert + R1_ADAPTER)[:READLEN]
    r2 = (revcomp(insert) + R2_ADAPTER)[:READLEN]
    return r1, r2


# (read_id, insert) — insert length drives the expected overlap behavior:
#   < READLEN  -> adapter read-through; both mates trim to len(insert)
#   >= READLEN -> full insert; forward overlap (offset>=0); no adapter trim
RECORDS = [
    (
        "pair1_readthrough36",
        "GATTACACCGGTATTACGCATGCAGTCAGTCAGTCA",
    ),  # 36bp insert -> trim to 36
    (
        "pair2_fullinsert60",
        "GATTACACCGGTATTACGCATGCAGTCAGTCAGTCATTGGAACCGGTTACGATCGATCG",
    ),  # 60bp >= 50 -> no trim
    (
        "pair3_readthrough45",
        "GATTACACCGGTATTACGCATGCAGTCAGTCAGTCATTGGAACCG",
    ),  # 45bp insert -> trim to 45
]


# Fallback fixture (data/qc/pe_fallback.*): exercises the 11-arg trim_adapters_pe
# adapter-by-sequence FALLBACK (fastp step 9), engaged when overlap analysis finds
# no adapter. Reads are given explicitly as (read_id, r1, r2):
#   ovwin_*  — a read-through pair: overlap (step 8) still wins; fastp + miint both
#              trim to the insert and never reach the fallback.
#   fb_*     — unrelated mates each carrying a clean 3' adapter: overlap fails, so
#              each mate is trimmed by adapter sequence.
# The common adapter AGATCGGAAGAGC is passed to fastp as both --adapter_sequence
# and --adapter_sequence_r2, and to trim_adapters_pe as the single-element list.
FALLBACK_ADAPTER = "AGATCGGAAGAGC"
FALLBACK_RECORDS = [
    ("ovwin_readthrough36",) + make_pair("GATTACACCGGTATTACGCATGCAGTCAGTCAGTCA"),
    (
        "fb_unrelated",
        "GATTACAGCATTGCATGGAACCTTACGATCG" + FALLBACK_ADAPTER,
        "TTGGCCAATTGGCCAATTACGTACGTACGTA" + FALLBACK_ADAPTER,
    ),
]


def main():
    out_dir = os.path.join(os.path.dirname(__file__), "..", "..", "data", "qc")
    out_dir = os.path.normpath(out_dir)
    os.makedirs(out_dir, exist_ok=True)

    r1_path = os.path.join(out_dir, "pe_overlap.r1.fq")
    r2_path = os.path.join(out_dir, "pe_overlap.r2.fq")

    with open(r1_path, "w") as f1, open(r2_path, "w") as f2:
        for rid, insert in RECORDS:
            assert len(insert) >= 30, f"{rid}: insert must be >= overlap_require (30)"
            r1, r2 = make_pair(insert)
            f1.write(f"@{rid}/1\n{r1}\n+\n{quals(len(r1))}\n")
            f2.write(f"@{rid}/2\n{r2}\n+\n{quals(len(r2))}\n")

    print(f"wrote {r1_path}")
    print(f"wrote {r2_path}")

    fb1_path = os.path.join(out_dir, "pe_fallback.r1.fq")
    fb2_path = os.path.join(out_dir, "pe_fallback.r2.fq")
    with open(fb1_path, "w") as f1, open(fb2_path, "w") as f2:
        for rid, r1, r2 in FALLBACK_RECORDS:
            f1.write(f"@{rid}/1\n{r1}\n+\n{quals(len(r1))}\n")
            f2.write(f"@{rid}/2\n{r2}\n+\n{quals(len(r2))}\n")

    print(f"wrote {fb1_path}")
    print(f"wrote {fb2_path}")


if __name__ == "__main__":
    main()
