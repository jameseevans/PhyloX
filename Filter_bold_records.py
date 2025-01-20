import csv
import argparse
from Bio import SeqIO

def main(csv_file, fasta_file, output_csv, output_fasta):
    with open(csv_file, mode='r') as csv_in:
        csv_reader = csv.DictReader(csv_in)
        csv_ids = {row['processid'].strip().lower() for row in csv_reader}

    fasta_ids = set()
    fasta_records = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        seq_id = record.id.split('|')[0].strip().lower()
        fasta_ids.add(seq_id)
        fasta_records.append(record)

    csv_only = csv_ids - fasta_ids
    fasta_only = fasta_ids - csv_ids

    print(f"Records in CSV but not in FASTA: {len(csv_only)}")
    print(f"Records in FASTA but not in CSV: {len(fasta_only)}")
    print(f"Example IDs only in CSV: {list(csv_only)[:5]}")
    print(f"Example IDs only in FASTA: {list(fasta_only)[:5]}")

    written_rows = 0
    with open(csv_file, mode='r') as csv_in, open(output_csv, mode='w', newline='') as csv_out:
        csv_reader = csv.DictReader(csv_in)
        csv_writer = csv.DictWriter(csv_out, fieldnames=csv_reader.fieldnames)
        csv_writer.writeheader()

        for row in csv_reader:
            processid = row['processid'].strip().lower()
            if processid in fasta_ids:
                csv_writer.writerow(row)
                written_rows += 1

    print(f"Number of rows written to the new CSV: {written_rows}")

    matched_fasta_records = [record for record in fasta_records if record.id.split('|')[0].strip().lower() in csv_ids]

    with open(output_fasta, mode='w') as fasta_out:
        SeqIO.write(matched_fasta_records, fasta_out, "fasta")

    print(f"Number of sequences written to the new FASTA: {len(matched_fasta_records)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter CSV and FASTA files based on matching IDs.")
    parser.add_argument("-c", "--csv", required=True, help="Input CSV file")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output-csv", required=True, help="Filtered output CSV file")
    parser.add_argument("-p", "--output-fasta", required=True, help="Filtered output FASTA file")
    args = parser.parse_args()

    main(args.csv, args.fasta, args.output_csv, args.output_fasta)
