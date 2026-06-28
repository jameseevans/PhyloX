#!/usr/bin/env python3

import argparse
from Bio import SeqIO

mito_prefixes = {'BIOD', 'CCCP', 'GBDL', 'HNSP', 'MIZA', 'QINL', 'SPSO', 'SRAA', 'BGLP', 'ZIPC', 'CDBP', 'SRAB', "BPST", "SRR", "GCA_", "ERR"}

def clean_sequence(seq_str):
    """Remove alignment gaps, then strip leading/trailing Ns (but not internal Ns)."""
    seq_str = seq_str.replace("-", "").replace(".", "")
    seq_str = seq_str.strip("Nn")
    return seq_str

def get_longest_sequence(sequences):
    mito_sequences = [seq for seq in sequences if any(seq.id.startswith(prefix) for prefix in mito_prefixes)]
    if mito_sequences:
        return max(mito_sequences, key=lambda seq: len(seq.seq))
    else:
        return max(sequences, key=lambda seq: len(seq.seq))

def remove_duplicates(input_fasta, output_fasta):
    sequences = SeqIO.parse(input_fasta, "fasta")
    unique_sequences = {}
    total = 0

    for seq_record in sequences:
        total += 1
        if total % 50000 == 0:
            print(f"  Processed {total:,} sequences...")

        cleaned = clean_sequence(str(seq_record.seq))

        if cleaned not in unique_sequences:
            unique_sequences[cleaned] = [seq_record]
        else:
            unique_sequences[cleaned].append(seq_record)

    print(f"  Done. {total:,} sequences parsed.")
    print(f"  {len(unique_sequences):,} unique sequences after deduplication ({total - len(unique_sequences):,} duplicates removed).")

    final_sequences = [get_longest_sequence(seq_records) for seq_records in unique_sequences.values()]

    with open(output_fasta, "w") as out_file:
        SeqIO.write(final_sequences, out_file, "fasta")

    print(f"Filtered alignment has been written to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a fasta alignment, retaining mitogenomes or the longest sequences.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta")
    parser.add_argument("-o", "--output", required=True, help="Output fasta")

    args = parser.parse_args()

    remove_duplicates(args.input, args.output)
