import argparse
from collections import defaultdict
import pandas as pd
from Bio import SeqIO

def parse_fasta(fasta_file):
    """Parse the FASTA file and return a dictionary of sequences with sequence id as key."""
    seq_dict = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).replace("-", "")
        seq_dict[record.id] = seq
    return seq_dict

def parse_metadata(metadata_file):
    """Parse the metadata CSV file into a pandas DataFrame."""
    metadata_df = pd.read_csv(metadata_file)
    return metadata_df

def get_metadata_priority(row, required_columns):
    """Calculate priority based on the completeness of specified columns."""
    return sum(pd.notna(row[col]) for col in required_columns)

def remove_duplicates(fasta_file, metadata_file, output_fasta_file, output_metadata_file, seq_id_column):
    """Remove duplicate sequences from the FASTA file and its associated metadata file."""
    
    seq_dict = parse_fasta(fasta_file)
    metadata_df = parse_metadata(metadata_file)
    
    seen_sequences = {}
    
    required_columns = ['subfamily_name', 'genus_name', 'species_name', 'subspecies_name']
    
    for _, row in metadata_df.iterrows():
        seq_id = row[seq_id_column]
        sequence = seq_dict.get(seq_id)
        
        if sequence:
            if sequence not in seen_sequences:
                seen_sequences[sequence] = row
            else:
                existing_row = seen_sequences[sequence]
                if get_metadata_priority(row, required_columns) > get_metadata_priority(existing_row, required_columns):
                    seen_sequences[sequence] = row

    filtered_seq_ids = [row[seq_id_column] for row in seen_sequences.values()]

    with open(output_fasta_file, 'w') as out_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in filtered_seq_ids:
                SeqIO.write(record, out_fasta, "fasta")
    
    filtered_metadata_df = metadata_df[metadata_df[seq_id_column].isin(filtered_seq_ids)]
    
    filtered_metadata_df.to_csv(output_metadata_file, index=False)

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from FASTA and metadata files.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file containing sequences.")
    parser.add_argument("-m", "--metadata", required=True, help="Input metadata CSV file.")
    parser.add_argument("-o", "--output_fasta", required=True, help="Output FASTA file for filtered sequences.")
    parser.add_argument("-d", "--output_metadata", required=True, help="Output metadata CSV file for filtered sequences.")
    parser.add_argument("-c", "--seq_id_column", default="processid", help="Column name in the metadata file for sequence IDs. Default is 'processid'.")

    args = parser.parse_args()

    remove_duplicates(args.fasta, args.metadata, args.output_fasta, args.output_metadata, args.seq_id_column)

if __name__ == "__main__":
    main()
