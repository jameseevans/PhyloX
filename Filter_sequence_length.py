#!/usr/bin/env python3
import sys
import csv

def filter_fasta(input_fasta, output_fasta, min_len, stats_csv):
    def read_fasta(fp):
        header = None
        seq_chunks = []
        for line in fp:
            line = line.rstrip()
            if line.startswith(">"):
                if header:
                    yield header, "".join(seq_chunks)
                header = line[1:]  # remove ">"
                seq_chunks = []
            else:
                seq_chunks.append(line)
        if header:
            yield header, "".join(seq_chunks)

    removed_stats = []

    with open(input_fasta) as inp, open(output_fasta, "w") as out:
        for header, seq in read_fasta(inp):
            original_len = len(seq)
            clean_seq = seq.replace("-", "").replace("N", "").replace("n", "")
            clean_len = len(clean_seq)

            if clean_len >= min_len:
                out.write(f">{header}\n{seq}\n")
            else:
                removed_stats.append({
                    "header": header,
                    "original_length": original_len,
                    "ungapped_nonN_length": clean_len,
                    "removed_reason": f"length<{min_len}"
                })

    # Write CSV report
    with open(stats_csv, "w", newline="") as csvfile:
        fieldnames = ["header", "original_length", "ungapped_nonN_length", "removed_reason"]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in removed_stats:
            writer.writerow(row)


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python Filter_sequence_length.py <input.fa> <output.fa> <min_length> <stats.csv>")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_fasta = sys.argv[2]
    min_len = int(sys.argv[3])
    stats_csv = sys.argv[4]

    filter_fasta(input_fasta, output_fasta, min_len, stats_csv)
