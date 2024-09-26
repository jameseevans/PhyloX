### Script to remove dupliacte Genbank accessions from genbank and metadata files ###

import argparse
from Bio import SeqIO
import pandas as pd

def main():
    parser = argparse.ArgumentParser(description="Filter GenBank and metadata files.")
    parser.add_argument("-g", "--input_genbank", required=True, help="Path to the input GenBank file.")
    parser.add_argument("-m", "--input_metadata", required=True, help="Path to the input metadata CSV file.")
    parser.add_argument("-og", "--output_genbank", required=True, help="Path to the output filtered GenBank file.")
    parser.add_argument("-om", "--output_metadata", required=True, help="Path to the output filtered metadata CSV file.")
    
    args = parser.parse_args()
    genbank_records = list(SeqIO.parse(args.input_genbank, "genbank"))
    metadata_df = pd.read_csv(args.input_metadata)
    blank_genbank_accession_df = metadata_df[metadata_df['genbank_accession'].isna()]
    non_blank_genbank_accession_df = metadata_df[metadata_df['genbank_accession'].notna()]
    deduplicated_non_blank_df = non_blank_genbank_accession_df.drop_duplicates(subset=['genbank_accession'], keep='first')
    filtered_metadata_df = pd.concat([blank_genbank_accession_df, deduplicated_non_blank_df], ignore_index=True)
    filtered_locus_names = filtered_metadata_df['db_id'].tolist()
    filtered_genbank_records = [record for record in genbank_records if record.name in filtered_locus_names]
    filtered_metadata_df.to_csv(args.output_metadata, index=False)
    with open(args.output_genbank, "w") as output_handle:
        SeqIO.write(filtered_genbank_records, output_handle, "genbank")

if __name__ == "__main__":
    main()
