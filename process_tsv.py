import csv
import argparse

def split_bold_tsv(tsv_path, fasta_path, metadata_path):
    with open(tsv_path, newline='', encoding='utf-8') as tsvfile, \
         open(fasta_path, 'w', encoding='utf-8') as fasta_out, \
         open(metadata_path, 'w', newline='', encoding='utf-8') as metadata_out:
        
        reader = csv.DictReader(tsvfile, delimiter='\t')

        metadata_fields = [f for f in reader.fieldnames if f not in ('processid', 'nuc')]
        writer = csv.DictWriter(metadata_out, fieldnames=metadata_fields)
        writer.writeheader()

        for row in reader:
            fasta_out.write(f">{row['record_id']}\n{row['nuc']}\n")
            
            metadata_row = {k: v for k, v in row.items() if k in metadata_fields}
            writer.writerow(metadata_row)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Split BOLD tsv into fasta and metadata csv")
    parser.add_argument('-i', '--input', required=True, help="Input BOLD tsv file")
    parser.add_argument('-f', '--fasta', required=True, help="Output fasta file")
    parser.add_argument('-m', '--metadata', required=True, help="Output metadata csv file")
    
    args = parser.parse_args()
    split_bold_tsv(args.input, args.fasta, args.metadata)
