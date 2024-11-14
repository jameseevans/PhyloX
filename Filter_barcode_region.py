import pandas as pd
from Bio import SeqIO
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Filter sequences with >=500 non-gap bases from alignment and metadata files, optionally removing duplicates based on non-gap content.")
parser.add_argument("-a", "--alignment_file", required=True, help="Path to the input FASTA alignment file")
parser.add_argument("-m", "--metadata_file", required=True, help="Path to the input CSV metadata file")
parser.add_argument("-g", "--genbank_file", required=True, help="Path to the input GenBank file")
parser.add_argument("-fm", "--filtered_metadata_file", required=True, help="Path to the output filtered CSV metadata file")
parser.add_argument("-fg", "--filtered_genbank_file", required=True, help="Path to the output filtered GenBank file")
parser.add_argument("--deduplicate", action="store_true", help="If set, remove duplicate sequences based on non-gap content, keeping the longest sequence.")

# Parse arguments
args = parser.parse_args()

# Define file paths from arguments
alignment_file = args.alignment_file
metadata_file = args.metadata_file
genbank_file = args.genbank_file
filtered_metadata_file = args.filtered_metadata_file
filtered_genbank_file = args.filtered_genbank_file
deduplicate = args.deduplicate

# Step 1: Identify sequences with >=500 non-gap bases and optionally remove duplicates in the alignment file
min_length = 500
valid_sequences = {}  # Dictionary to store unique non-gap sequences and their IDs

for record in SeqIO.parse(alignment_file, "fasta"):
    # Remove gaps and get the non-gap content
    non_gap_sequence = str(record.seq).replace("-", "")
    
    # Count non-gap characters
    non_gap_count = len(non_gap_sequence)
    
    # Only consider sequences with enough non-gap bases
    if non_gap_count >= min_length:
        if deduplicate:
            # Deduplication mode: keep only the longest version of each unique non-gap sequence
            if non_gap_sequence in valid_sequences:
                # Retain the sequence with the longest total length
                if len(record.seq) > len(valid_sequences[non_gap_sequence]["seq"]):
                    valid_sequences[non_gap_sequence] = {"id": record.id, "seq": record.seq}
            else:
                valid_sequences[non_gap_sequence] = {"id": record.id, "seq": record.seq}
        else:
            # No deduplication: store each sequence ID directly if it meets the criteria
            valid_sequences[record.id] = {"id": record.id, "seq": record.seq}

# Extract the final list of valid sequence IDs
final_valid_ids = [data["id"] for data in valid_sequences.values()]

# Step 2: Filter the metadata CSV to retain only rows with valid sequences
metadata = pd.read_csv(metadata_file)
filtered_metadata = metadata[metadata['db_id'].isin(final_valid_ids)]
filtered_metadata.to_csv(filtered_metadata_file, index=False)
print(f"Filtered metadata saved to {filtered_metadata_file}")

# Step 3: Filter the GenBank file to retain only records with valid sequence names and, if --deduplicate, the longest versions
if deduplicate:
    longest_genbank_records = {}
    for record in SeqIO.parse(genbank_file, "genbank"):
        if record.name in final_valid_ids:
            non_gap_sequence = str(record.seq).replace("-", "")
            if non_gap_sequence in longest_genbank_records:
                if len(record.seq) > len(longest_genbank_records[non_gap_sequence].seq):
                    longest_genbank_records[non_gap_sequence] = record
            else:
                longest_genbank_records[non_gap_sequence] = record
    genbank_records_to_write = longest_genbank_records.values()
else:
    # Without deduplication, just keep the records matching final_valid_ids
    genbank_records_to_write = [
        record for record in SeqIO.parse(genbank_file, "genbank") if record.name in final_valid_ids
    ]

# Save the filtered GenBank records to a new GenBank file
with open(filtered_genbank_file, "w") as output_handle:
    SeqIO.write(genbank_records_to_write, output_handle, "genbank")

print(f"Filtered GenBank file saved to {filtered_genbank_file}")
