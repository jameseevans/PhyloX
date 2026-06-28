#!/usr/bin/env python3

import argparse
import re
from Bio import SeqIO


def clean_sequence(seq_str):
    """Remove alignment gaps, then strip leading/trailing Ns (but not internal Ns)."""
    seq_str = seq_str.replace("-", "").replace(".", "")
    seq_str = seq_str.strip("Nn")
    return seq_str


def find_duplicates(input_fasta, output_file, min_length):
    seq_dict = {}
    total = 0

    print("Parsing sequences...")
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        total += 1
        if total % 50000 == 0:
            print(f"  Processed {total:,} sequences...")

        cleaned = clean_sequence(str(seq_record.seq))

        if len(cleaned) < min_length:
            continue

        if cleaned not in seq_dict:
            seq_dict[cleaned] = []
        seq_dict[cleaned].append(seq_record.id)

    print(f"  Done. {total:,} sequences parsed.")

    duplicate_groups = {seq: ids for seq, ids in seq_dict.items() if len(ids) > 1}

    n_groups = len(duplicate_groups)
    n_duplicate_seqs = sum(len(ids) for ids in duplicate_groups.values())

    print(f"\nFound {n_groups:,} duplicate groups ({n_duplicate_seqs:,} sequences involved).")

    with open(output_file, "w") as out:
        out.write("representative\tduplicates\tn_duplicates\n")
        for ids in sorted(duplicate_groups.values(), key=len, reverse=True):
            representative = ids[0]
            duplicates = ids[1:]
            out.write(f"{representative}\t{','.join(duplicates)}\t{len(ids)}\n")

    print(f"Duplicate list written to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=(
            "Find duplicate sequences in a fasta alignment. "
            "Gaps and terminal Ns are removed before comparison."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input fasta alignment")
    parser.add_argument("-o", "--output", required=True, help="Output TSV file")
    parser.add_argument(
        "--min-length", type=int, default=1,
        help="Minimum cleaned sequence length to consider (default: 1)"
    )

    args = parser.parse_args()
    find_duplicates(args.input, args.output, args.min_length)
