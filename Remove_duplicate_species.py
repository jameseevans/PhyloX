import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser(description="Filter GenBank and metadata files to retain the longest sequence per species, keeping sequences without species ID.")
parser.add_argument("-g", "--genbank", required=True, help="Input GenBank file.")
parser.add_argument("-m", "--metadata", required=True, help="Input metadata CSV file.")
parser.add_argument("-og", "--output_genbank", required=True, help="Output GenBank file.")
parser.add_argument("-om", "--output_metadata", required=True, help="Output filtered metadata CSV file.")

args = parser.parse_args()

metadata_df = pd.read_csv(args.metadata)

sequences = []
print("Parsing GenBank file...")

genbank_records = list(SeqIO.parse(args.genbank, 'genbank'))

for record in genbank_records:
    sequences.append({
        'locus': record.name,
        'sequence': record.seq,
        'length': len(record.seq)
    })

sequences_df = pd.DataFrame(sequences)

print(f"Parsed {len(sequences_df)} sequences from GenBank file.")

merged_df = pd.merge(metadata_df, sequences_df, left_on='db_id', right_on='locus', how='left')

print(f"Merged dataset contains {len(merged_df)} records.")

with_species = merged_df[merged_df['species'].notna()]
without_species = merged_df[merged_df['species'].isna()]

print(f"Records with species: {len(with_species)}, without species: {len(without_species)}")

try:
    filtered_with_species = with_species.loc[with_species.groupby('species')['length_y'].idxmax()]
    print(f"Filtered dataset with species contains {len(filtered_with_species)} unique species records.")
except KeyError as e:
    print(f"KeyError: {e}")
    print("Available columns:", merged_df.columns)
    exit(1)

filtered_df = pd.concat([filtered_with_species, without_species])

print(f"Final filtered dataset contains {len(filtered_df)} records (with and without species).")

filtered_records = []
print("Filtering GenBank records...")

genbank_dict = {record.name: record for record in genbank_records}

for index, row in filtered_df.iterrows():
    if row['locus'] in genbank_dict:
        filtered_records.append(genbank_dict[row['locus']])

print(f"Filtered {len(filtered_records)} sequences.")

with open(args.output_genbank, 'w') as output_handle:
    SeqIO.write(filtered_records, output_handle, 'genbank')

print(f"Filtered sequences written to {args.output_genbank}")

filtered_metadata_df = filtered_df.drop(columns=['sequence'])

filtered_metadata_df.to_csv(args.output_metadata, index=False)

print(f"Filtered metadata saved to {args.output_metadata}")
