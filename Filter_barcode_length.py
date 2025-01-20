import argparse
import csv
from Bio import SeqIO

def filter_metadata(fasta_file, metadata_file, output_file, min_length):
    sequence_lengths = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        sequence = str(record.seq).replace("-", "")
        sequence_lengths[record.id] = len(sequence)
    
    print(f"Parsed {len(sequence_lengths)} sequences from FASTA file.")
    
    filtered_metadata = []
    with open(metadata_file, mode="r") as csvfile:
        reader = csv.DictReader(csvfile)
        fieldnames = reader.fieldnames
        print(f"Metadata columns: {fieldnames}")
        for row in reader:
            processid = row["processid"]
            if processid not in sequence_lengths:
                print(f"processid {processid} not found in FASTA sequences.")
            elif sequence_lengths[processid] < min_length:
                print(f"Sequence {processid} length {sequence_lengths[processid]} is less than {min_length}.")
            else:
                print(f"Adding processid {processid} with length {sequence_lengths[processid]} to filtered metadata.")
                filtered_metadata.append(row)
    
    if filtered_metadata:
        with open(output_file, mode="w", newline="") as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(filtered_metadata)
    else:
        print("No metadata rows passed the filtering criteria. Output file will only contain the header.")


def main():
    parser = argparse.ArgumentParser(description="Filter metadata based on sequence lengths.")
    parser.add_argument("-f", "--fasta", required=True, help="Input FASTA file.")
    parser.add_argument("-m", "--metadata", required=True, help="Input metadata CSV file.")
    parser.add_argument("-o", "--output", required=True, help="Output filtered metadata CSV file.")
    parser.add_argument("-l", "--length", type=int, default=500, help="Minimum sequence length (default: 500).")
    args = parser.parse_args()

    filter_metadata(args.fasta, args.metadata, args.output, args.length)

if __name__ == "__main__":
    main()
