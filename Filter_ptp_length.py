import re
from Bio import SeqIO
import pandas as pd
import argparse

# Argument parser
parser = argparse.ArgumentParser(description="Filter FASTA and metadata based on MPTP output.")
parser.add_argument("--mptp_file", required=True, help="Path to the MPTP output file.")
parser.add_argument("--fasta_file", required=True, help="Path to the input FASTA file.")
parser.add_argument("--metadata_file", required=True, help="Path to the input metadata CSV file.")
parser.add_argument("--output_fasta", required=True, help="Path to the output filtered FASTA file.")
parser.add_argument("--output_metadata", required=True, help="Path to the output filtered metadata CSV file.")
args = parser.parse_args()

# Read MPTP output and extract species information
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

# Extract sequence IDs from FASTA file, using only the first part before '|'
sequences = {}
seq_id_mapping = {}  # Map full FASTA IDs to processid
for record in SeqIO.parse(args.fasta_file, "fasta"):
    processid = record.id.split("|")[0]  # Extract only the first part before '|'
    sequences[processid] = record
    seq_id_mapping[record.id] = processid  # Map full ID to processid

# Select longest sequence for each species
selected_sequences = {}
for species, seq_ids in species_dict.items():
    longest_seq = None
    max_length = 0
    for seq_id in seq_ids:
        processid = seq_id_mapping.get(seq_id, seq_id)  # Convert to processid if needed
        if processid in sequences:
            non_gap_length = len(str(sequences[processid].seq).replace("-", ""))
            if non_gap_length > max_length:
                longest_seq = sequences[processid]
                max_length = non_gap_length
    if longest_seq:
        selected_sequences[longest_seq.id] = longest_seq  # Store using full ID for output

# Write filtered FASTA file
with open(args.output_fasta, "w") as output:
    SeqIO.write(selected_sequences.values(), output, "fasta")

# Read metadata
metadata = pd.read_csv(args.metadata_file)

# Filter metadata using extracted process IDs
filtered_metadata = metadata[metadata['processid'].isin(set(seq_id_mapping[id] for id in selected_sequences.keys()))]

# Write filtered metadata to file
filtered_metadata.to_csv(args.output_metadata, index=False)

print(f"Filtered FASTA written to {args.output_fasta}")
print(f"Filtered metadata written to {args.output_metadata}")
