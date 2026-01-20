#!/usr/bin/env python3
#Processes a TSV file from BOLDsystems v5, to extract specified loci and remove invalid bases and gaps in a fasta file. Outputs a filtered tsv file.

import argparse
import csv
import sys

VALID_BASES = {
    'A', 'T', 'C', 'G', 'N', 'U',
    'R', 'Y', 'S', 'W', 'K', 'M',
    'B', 'D', 'H', 'V'
}

GAP_CHARS = {'-', '.'}

def clean_sequence(seq):
    """
    Remove gaps and replace non-valid bases with N
    """
    cleaned = []
    for base in seq.upper():
        if base in GAP_CHARS:
            continue
        elif base in VALID_BASES:
            cleaned.append(base)
        else:
            cleaned.append('N')
    return ''.join(cleaned)

def tsv_to_fasta(
    input_tsv,
    output_fasta,
    output_tsv=None,
    locus=None
):
    with open(input_tsv, newline="", encoding="utf-8") as infile:
        reader = csv.DictReader(infile, delimiter="\t")

        required_cols = {"processid", "nuc", "marker_code"}
        missing = required_cols - set(reader.fieldnames)
        if missing:
            sys.exit(f"ERROR: Missing required columns: {', '.join(missing)}")

        rows = []
        fasta_count = 0

        for row in reader:
            # Filter by locus if requested
            if locus and row["marker_code"] != locus:
                continue

            raw_seq = row["nuc"].strip()
            if not raw_seq:
                continue

            cleaned_seq = clean_sequence(raw_seq)
            if not cleaned_seq:
                continue

            row["nuc"] = cleaned_seq
            rows.append(row)

        if not rows:
            sys.exit("ERROR: No records passed filters")

        # Write FASTA
        with open(output_fasta, "w", encoding="utf-8") as fasta:
            for row in rows:
                seq_id = row["processid"].strip()
                seq = row["nuc"]

                fasta.write(f">{seq_id}\n")
                for i in range(0, len(seq), 60):
                    fasta.write(seq[i:i+60] + "\n")

                fasta_count += 1

        # Write filtered TSV if requested
        if output_tsv:
            with open(output_tsv, "w", newline="", encoding="utf-8") as outtsv:
                writer = csv.DictWriter(
                    outtsv,
                    fieldnames=reader.fieldnames,
                    delimiter="\t"
                )
                writer.writeheader()
                writer.writerows(rows)

    print(f"Written {fasta_count} sequences to {output_fasta}")
    if output_tsv:
        print(f"Written filtered TSV to {output_tsv}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert BOLD TSV to FASTA, filter by locus, clean sequences"
    )

    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input BOLD TSV file"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output FASTA file"
    )
    parser.add_argument(
        "-f", "--filtered-tsv",
        help="Write filtered records to this TSV file"
    )
    parser.add_argument(
        "-l", "--locus",
        help="Filter by marker_code (e.g. COI-5P)"
    )

    args = parser.parse_args()

    tsv_to_fasta(
        input_tsv=args.input,
        output_fasta=args.output,
        output_tsv=args.filtered_tsv,
        locus=args.locus
    )

if __name__ == "__main__":
    main()
