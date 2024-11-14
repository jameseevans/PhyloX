import argparse
import pandas as pd
from Bio import SeqIO

# Set up command line argument parsing
parser = argparse.ArgumentParser(description="Filter GenBank and metadata files to retain the longest sequence per species, keeping sequences without species ID.")
parser.add_argument("-g", "--genbank", required=True, help="Input GenBank file.")
parser.add_argument("-m", "--metadata", required=True, help="Input metadata CSV file.")
parser.add_argument("-og", "--output_genbank", required=True, help="Output GenBank file.")
parser.add_argument("-om", "--output_metadata", required=True, help="Output filtered metadata CSV file.")

args = parser.parse_args()

# Load the metadata CSV file
metadata_df = pd.read_csv(args.metadata)

# Parse the GenBank file and extract sequence details (locus, length, and sequence)
sequences = []
print("Parsing GenBank file...")

genbank_records = list(SeqIO.parse(args.genbank, 'genbank'))  # Load all GenBank records into memory

for record in genbank_records:
    sequences.append({
        'locus': record.name,  # LOCUS
        'sequence': record.seq,
        'length': len(record.seq)  # Calculate length directly from the sequence
    })

# Convert sequences to a DataFrame for easier processing
sequences_df = pd.DataFrame(sequences)

# Verify sequences were parsed
print(f"Parsed {len(sequences_df)} sequences from GenBank file.")

# Merge sequences with metadata using db_id from metadata and locus from GenBank file
merged_df = pd.merge(metadata_df, sequences_df, left_on='db_id', right_on='locus', how='left')

# Check the number of records after merging
print(f"Merged dataset contains {len(merged_df)} records.")

# Separate rows with and without a species ID
with_species = merged_df[merged_df['species'].notna()]
without_species = merged_df[merged_df['species'].isna()]

# Print the count of sequences with and without species
print(f"Records with species: {len(with_species)}, without species: {len(without_species)}")

# For rows with species, retain only the longest sequence for each species
try:
    filtered_with_species = with_species.loc[with_species.groupby('species')['length_y'].idxmax()]
    print(f"Filtered dataset with species contains {len(filtered_with_species)} unique species records.")
except KeyError as e:
    print(f"KeyError: {e}")
    print("Available columns:", merged_df.columns)
    exit(1)

# Concatenate sequences with species and those without species
filtered_df = pd.concat([filtered_with_species, without_species])

# Print the final count of filtered records
print(f"Final filtered dataset contains {len(filtered_df)} records (with and without species).")

# Save the filtered sequences to a new GenBank file
filtered_records = []
print("Filtering GenBank records...")

# Create a lookup dictionary for GenBank records by locus
genbank_dict = {record.name: record for record in genbank_records}

for index, row in filtered_df.iterrows():
    # Find the corresponding sequence record in the preloaded GenBank records
    if row['locus'] in genbank_dict:
        filtered_records.append(genbank_dict[row['locus']])

# Verify that records were filtered correctly
print(f"Filtered {len(filtered_records)} sequences.")

# Write the filtered sequences to the output GenBank file
with open(args.output_genbank, 'w') as output_handle:
    SeqIO.write(filtered_records, output_handle, 'genbank')

print(f"Filtered sequences written to {args.output_genbank}")

# Drop the 'sequence' column from the metadata before saving
filtered_metadata_df = filtered_df.drop(columns=['sequence'])

# Save the filtered metadata to the output CSV file
filtered_metadata_df.to_csv(args.output_metadata, index=False)

print(f"Filtered metadata saved to {args.output_metadata}")
