import re
from Bio import SeqIO
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description="Filter FASTA and metadata based on MPTP output.")
parser.add_argument("--mptp_file", required=True, help="Path to the MPTP output file.")
parser.add_argument("--fasta_file", required=True, help="Path to the input FASTA file.")
parser.add_argument("--metadata_file", required=True, help="Path to the input metadata CSV file.")
parser.add_argument("--output_fasta", required=True, help="Path to the output filtered FASTA file.")
parser.add_argument("--output_metadata", required=True, help="Path to the output filtered metadata CSV file.")
args = parser.parse_args()

species_dict = {}

with open(args.mptp_file, "r") as file:
    current_species = None
    for line in file:
        if line.startswith("Species "):
            match = re.search(r"Species (\d+):", line)
            if match:
                current_species = int(match.group(1))
                species_dict[current_species] = []
        elif current_species and line.strip():
            species_dict[current_species].extend(line.strip().split(", "))

selected_sequences = {}

sequences = {rec.id: rec for rec in SeqIO.parse(args.fasta_file, "fasta")}

for species, seq_ids in species_dict.items():
    longest_seq = None
    max_length = 0
    for seq_id in seq_ids:
        if seq_id in sequences:
            non_gap_length = len(str(sequences[seq_id].seq).replace("-", ""))
            if non_gap_length > max_length:
                longest_seq = sequences[seq_id]
                max_length = non_gap_length
    if longest_seq:
        selected_sequences[longest_seq.id] = longest_seq

with open(args.output_fasta, "w") as output:
    SeqIO.write(selected_sequences.values(), output, "fasta")

metadata = pd.read_csv(args.metadata_file)

filtered_metadata = metadata[metadata['processid'].isin(selected_sequences.keys())]

filtered_metadata.to_csv(args.output_metadata, index=False)

print(f"Filtered FASTA written to {args.output_fasta}")
print(f"Filtered metadata written to {args.output_metadata}")
