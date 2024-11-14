import argparse
from Bio import SeqIO
import pandas as pd
from collections import defaultdict

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter duplicate sequences from metadata CSV and GenBank files based on taxonomy information.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file of amino acid sequences for identifying duplicates.")
    parser.add_argument("-c", "--csv", required=True, help="Input metadata CSV file with sequence names in the 'db_id' column.")
    parser.add_argument("-g", "--genbank", required=True, help="Input GenBank file with sequence names under 'LOCUS'.")
    parser.add_argument("--csv_output", required=True, help="Output file name for the filtered metadata CSV.")
    parser.add_argument("--genbank_output", required=True, help="Output file name for the filtered GenBank file.")
    return parser.parse_args()

def load_fasta_sequences(fasta_file):
    """Loads sequences from FASTA, removing gaps and detecting duplicates based on non-gap sequences."""
    sequences = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).replace("-", "")  # Remove gaps for duplicate detection
        if seq not in sequences:
            sequences[seq] = []
        sequences[seq].append(record)
    return sequences

def load_metadata(metadata_file):
    df = pd.read_csv(metadata_file)
    taxonomy_cols = ["subspecies", "species", "subgenus", "genus", "subtribe", "tribe", "subfamily", "family", "superfamily", "infraorder", "suborder", "order", "class"]
    return df, taxonomy_cols

def load_genbank_records(genbank_file):
    genbank_records = {record.name: record for record in SeqIO.parse(genbank_file, "genbank")}
    return genbank_records

def count_taxonomy_data(row, taxonomy_cols):
    """Counts non-empty taxonomy fields in a metadata row to determine completeness of taxonomy information."""
    return sum(not pd.isna(row[col]) and row[col] != "" for col in taxonomy_cols)

def filter_duplicates(sequences, metadata_df, taxonomy_cols):
    """Filters duplicates based on taxonomy richness, retaining the entry with the most complete data for each sequence."""
    filtered_metadata_indices = []
    for seq, records in sequences.items():
        if len(records) == 1:
            # Only one record, so keep it
            filtered_metadata_indices.append(records[0].id)
            continue
        # For duplicate sequences, select the metadata entry with the richest taxonomy information
        record_ids = [rec.id for rec in records]
        metadata_subset = metadata_df[metadata_df["db_id"].isin(record_ids)]
        metadata_subset["taxonomy_score"] = metadata_subset.apply(count_taxonomy_data, axis=1, taxonomy_cols=taxonomy_cols)
        best_record = metadata_subset.sort_values("taxonomy_score", ascending=False).iloc[0]
        filtered_metadata_indices.append(best_record["db_id"])
    return metadata_df[metadata_df["db_id"].isin(filtered_metadata_indices)]

def write_output_files(filtered_metadata_df, genbank_records, csv_output, genbank_output):
    # Write the filtered metadata CSV file
    filtered_metadata_df.to_csv(csv_output, index=False)

    # Write the filtered GenBank file based on filtered metadata
    filtered_genbank_records = [genbank_records[seq_id] for seq_id in filtered_metadata_df["db_id"] if seq_id in genbank_records]
    with open(genbank_output, "w") as genbank_out:
        SeqIO.write(filtered_genbank_records, genbank_out, "genbank")

def main():
    args = parse_arguments()
    sequences = load_fasta_sequences(args.fasta)
    metadata_df, taxonomy_cols = load_metadata(args.csv)
    genbank_records = load_genbank_records(args.genbank)
    
    # Filter duplicates in metadata based on taxonomy richness
    filtered_metadata = filter_duplicates(sequences, metadata_df, taxonomy_cols)
    
    # Write the output filtered files
    write_output_files(filtered_metadata, genbank_records, args.csv_output, args.genbank_output)

if __name__ == "__main__":
    main()
