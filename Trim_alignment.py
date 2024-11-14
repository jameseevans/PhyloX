import argparse
from Bio import SeqIO
import numpy as np

def trim_alignment(input_alignment, input_names, output_file):
    # Load the list of sequence names from the text file
    with open(input_names, "r") as f:
        target_names = set(line.strip() for line in f)

    # Parse the alignment
    alignment_records = list(SeqIO.parse(input_alignment, "fasta"))

    # Find sequences that match the target names
    target_seqs = [record for record in alignment_records if record.id in target_names]

    # Find the column indices that span the region of interest
    alignment_length = len(alignment_records[0].seq)
    start, end = alignment_length, 0

    # Loop through each target sequence to find the min and max bounds
    for record in target_seqs:
        sequence_array = np.array(record.seq)
        non_gap_positions = np.where(sequence_array != "-")[0]
        if len(non_gap_positions) > 0:
            start = min(start, non_gap_positions[0])
            end = max(end, non_gap_positions[-1])

    # Trim each sequence in the alignment to the specified region
    trimmed_records = []
    for record in alignment_records:
        trimmed_seq = record.seq[start:end+1]  # Get the span of columns within the region
        record.seq = trimmed_seq
        trimmed_records.append(record)

    # Write the trimmed alignment to a new file
    with open(output_file, "w") as output_handle:
        SeqIO.write(trimmed_records, output_handle, "fasta")

    print(f"Trimmed alignment saved to {output_file}")

# Set up command-line argument parsing
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Trim alignment to preserve columns spanning a specified region.")
    parser.add_argument("-a", "--alignment", required=True, help="Path to the input alignment FASTA file")
    parser.add_argument("-n", "--names", required=True, help="Path to the text file with sequence names defining the region")
    parser.add_argument("-o", "--output", required=True, help="Path to save the trimmed alignment")

    args = parser.parse_args()

    # Run the trimming function with the specified arguments
    trim_alignment(args.alignment, args.names, args.output)
