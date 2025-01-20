import pandas as pd
from Bio import SeqIO
import argparse

parser = argparse.ArgumentParser(description="Filter sequences with >=500 non-gap bases from alignment and metadata files, optionally removing duplicates based on non-gap content.")
parser.add_argument("-a", "--alignment_file", required=True, help="Path to the input FASTA alignment file")
parser.add_argument("-m", "--metadata_file", required=True, help="Path to the input CSV metadata file")
parser.add_argument("-g", "--genbank_file", required=True, help="Path to the input GenBank file")
parser.add_argument("-fm", "--filtered_metadata_file", required=True, help="Path to the output filtered CSV metadata file")
parser.add_argument("-fg", "--filtered_genbank_file", required=True, help="Path to the output filtered GenBank file")
parser.add_argument("--deduplicate", action="store_true", help="If set, remove duplicate sequences based on non-gap content, keeping the longest sequence.")

args = parser.parse_args()

alignment_file = args.alignment_file
metadata_file = args.metadata_file
genbank_file = args.genbank_file
filtered_metadata_file = args.filtered_metadata_file
filtered_genbank_file = args.filtered_genbank_file
deduplicate = args.deduplicate

min_length = 500
valid_sequences = {}

for record in SeqIO.parse(alignment_file, "fasta"):
    non_gap_sequence = str(record.seq).replace("-", "")
    
    non_gap_count = len(non_gap_sequence)
    
    if non_gap_count >= min_length:
        if deduplicate:
            if non_gap_sequence in valid_sequences:
                if len(record.seq) > len(valid_sequences[non_gap_sequence]["seq"]):
                    valid_sequences[non_gap_sequence] = {"id": record.id, "seq": record.seq}
            else:
                valid_sequences[non_gap_sequence] = {"id": record.id, "seq": record.seq}
        else:
            valid_sequences[record.id] = {"id": record.id, "seq": record.seq}

final_valid_ids = [data["id"] for data in valid_sequences.values()]

metadata = pd.read_csv(metadata_file)
filtered_metadata = metadata[metadata['db_id'].isin(final_valid_ids)]
filtered_metadata.to_csv(filtered_metadata_file, index=False)
print(f"Filtered metadata saved to {filtered_metadata_file}")

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
    genbank_records_to_write = [
        record for record in SeqIO.parse(genbank_file, "genbank") if record.name in final_valid_ids
    ]

with open(filtered_genbank_file, "w") as output_handle:
    SeqIO.write(genbank_records_to_write, output_handle, "genbank")

print(f"Filtered GenBank file saved to {filtered_genbank_file}")
