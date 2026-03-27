#!/usr/bin/env python3
"""Generate synthetic FASTA test data for UCHIME chimera detection tests.

Creates reference sequences, chimeric queries, and de novo input with
abundance annotations. Sequences are designed so chimeric breakpoints
are unambiguous and vsearch produces deterministic output.

Design: 6 references (~300bp), each built from a common template with
~15-20% divergence between any pair. Chimeras are clean joins at known
positions with no additional mutations.
"""

import random
import os

SEED = 42
SEQ_LEN = 300
DATA_DIR = os.path.join(os.path.dirname(__file__), "..", "..", "data", "uchime")


def random_seq(length, rng):
    return "".join(rng.choice("ACGT") for _ in range(length))


def mutate(seq, positions, rng):
    """Mutate seq at given positions to a different base."""
    bases = list(seq)
    for pos in positions:
        orig = bases[pos]
        choices = [b for b in "ACGT" if b != orig]
        bases[pos] = rng.choice(choices)
    return "".join(bases)


def write_fasta(path, records):
    """Write list of (header, sequence) tuples to FASTA."""
    with open(path, "w") as f:
        for header, seq in records:
            f.write(f">{header}\n{seq}\n")


def main():
    rng = random.Random(SEED)
    os.makedirs(DATA_DIR, exist_ok=True)

    # Generate a common template
    template = random_seq(SEQ_LEN, rng)

    # Create 6 references, each with ~50 mutations from template
    # (this gives ~75-85% identity between pairs of references)
    all_positions = list(range(SEQ_LEN))
    refs = []
    for i in range(6):
        rng_ref = random.Random(SEED + i + 1)
        # Pick 50 unique positions to mutate
        mut_positions = sorted(rng_ref.sample(all_positions, 50))
        ref_seq = mutate(template, mut_positions, rng_ref)
        refs.append((f"ref{i+1}", ref_seq))

    # Write reference DB
    write_fasta(os.path.join(DATA_DIR, "chimera_ref.fasta"), refs)

    # Create query sequences
    queries = []

    # 1. Clean (non-chimeric): identical to ref1
    queries.append(("query_clean1", refs[0][1]))

    # 2. Clean: ref3 with 2 small substitutions
    rng_clean2 = random.Random(SEED + 100)
    clean2_positions = rng_clean2.sample(all_positions, 2)
    queries.append(("query_clean2", mutate(refs[2][1], clean2_positions, rng_clean2)))

    # 3. Clear chimera: first half of ref1 + second half of ref2
    mid = SEQ_LEN // 2  # position 150
    chimera1 = refs[0][1][:mid] + refs[1][1][mid:]
    queries.append(("query_chimera1", chimera1))

    # 4. Asymmetric chimera: first 40% of ref3 + last 60% of ref4
    cut_40 = int(SEQ_LEN * 0.4)  # position 120
    chimera2 = refs[2][1][:cut_40] + refs[3][1][cut_40:]
    queries.append(("query_chimera2", chimera2))

    # 5. Chimera with 1 substitution near crossover: ref5 front + ref6 back
    chimera3_raw = refs[4][1][:mid] + refs[5][1][mid:]
    # Add 1 substitution near the crossover (position mid+5)
    rng_c3 = random.Random(SEED + 200)
    chimera3 = mutate(chimera3_raw, [mid + 5], rng_c3)
    queries.append(("query_chimera3", chimera3))

    # 6. No parent: completely random sequence
    queries.append(("query_noparent", random_seq(SEQ_LEN, random.Random(SEED + 300))))

    # 7. Divergent chimera: chimera from ref1+ref2 with 5% extra mutations
    rng_div = random.Random(SEED + 400)
    div_positions = rng_div.sample(all_positions, 15)  # 5% of 300
    divergent = mutate(chimera1, div_positions, rng_div)
    queries.append(("query_divergent_chimera", divergent))

    # 8. Short query: 60bp from ref1
    queries.append(("query_short", refs[0][1][:60]))

    write_fasta(os.path.join(DATA_DIR, "chimera_queries.fasta"), queries)

    # Create de novo input (abundance-annotated)
    # Order doesn't matter in the file; vsearch sorts by abundance internally.
    # We need: 3 abundant non-chimeras, 2 low-abundance chimeras, 1 low non-chimera.
    denovo = []

    # High abundance non-chimeras (based on ref1, ref2, ref3)
    denovo.append((f"seq1;size=1000", refs[0][1]))
    denovo.append((f"seq2;size=800", refs[1][1]))
    denovo.append((f"seq3;size=500", refs[2][1]))

    # Low abundance chimera: first half of seq1(ref1) + second half of seq2(ref2)
    denovo.append((f"seq4;size=10", refs[0][1][:mid] + refs[1][1][mid:]))

    # Low abundance chimera: first half of seq2(ref2) + second half of seq3(ref3)
    denovo.append((f"seq5;size=5", refs[1][1][:mid] + refs[2][1][mid:]))

    # Low abundance non-chimera (based on ref4)
    denovo.append((f"seq6;size=3", refs[3][1]))

    write_fasta(os.path.join(DATA_DIR, "denovo_input.fasta"), denovo)

    print(f"Generated test data in {DATA_DIR}:")
    print(f"  chimera_ref.fasta     - 6 reference sequences ({SEQ_LEN}bp each)")
    print(f"  chimera_queries.fasta - 8 query sequences")
    print(f"  denovo_input.fasta    - 6 abundance-annotated sequences")
    print()

    # Print pairwise identity matrix for references (rough)
    print("Approximate pairwise identity between references:")
    for i in range(6):
        for j in range(i + 1, 6):
            matches = sum(1 for a, b in zip(refs[i][1], refs[j][1]) if a == b)
            pct = 100.0 * matches / SEQ_LEN
            print(f"  ref{i+1} vs ref{j+1}: {pct:.1f}%")


if __name__ == "__main__":
    main()
