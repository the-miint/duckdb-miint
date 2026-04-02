#!/usr/bin/env python3
"""Generate synthetic paired-end FASTQ reads from FASTA sequences.

Splits each input sequence at a known position to create overlapping
forward/reverse read pairs with uniform Q30 quality. The reverse read
is reverse-complemented (as produced by Illumina sequencers).

Usage:
    python generate_test_pairs.py input.fasta output_R1.fq output_R2.fq [--overlap 50] [--max-seqs 50]
"""
import sys

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def reverse_complement(seq):
    return seq.translate(COMPLEMENT)[::-1]

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("input_fasta")
    parser.add_argument("output_r1")
    parser.add_argument("output_r2")
    parser.add_argument("--overlap", type=int, default=50)
    parser.add_argument("--max-seqs", type=int, default=50)
    args = parser.parse_args()

    # Parse FASTA
    sequences = []
    current_id = None
    current_seq = []
    with open(args.input_fasta) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_id is not None:
                    sequences.append((current_id, "".join(current_seq)))
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            sequences.append((current_id, "".join(current_seq)))

    sequences = sequences[:args.max_seqs]

    with open(args.output_r1, "w") as r1, open(args.output_r2, "w") as r2:
        for seq_id, seq in sequences:
            seq = seq.upper()
            # Skip sequences too short for the overlap
            min_read_len = args.overlap + 20
            if len(seq) < min_read_len * 2 - args.overlap:
                continue

            # Split: fwd gets first half + overlap, rev gets overlap + second half
            mid = len(seq) // 2
            fwd_end = mid + args.overlap // 2
            rev_start = mid - args.overlap // 2

            fwd_seq = seq[:fwd_end]
            rev_seq = seq[rev_start:]  # This is the forward strand; we'll RC it

            qual = "I" * len(fwd_seq)  # Q40 (ASCII 73)
            r1.write(f"@{seq_id}\n{fwd_seq}\n+\n{qual}\n")

            rev_rc = reverse_complement(rev_seq)
            qual_rev = "I" * len(rev_rc)
            r2.write(f"@{seq_id}\n{rev_rc}\n+\n{qual_rev}\n")

    print(f"Generated {len(sequences)} read pairs", file=sys.stderr)

if __name__ == "__main__":
    main()
