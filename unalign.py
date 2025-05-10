
import argparse
from Bio import SeqIO

def remove_gaps(input_file, output_file):
    with open(output_file, "w") as out_handle:
        for record in SeqIO.parse(input_file, "fasta"):
            record.seq = record.seq.ungap("-")
            SeqIO.write(record, out_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Remove gaps (-) from sequences in a fasta file.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta")
    parser.add_argument("-o", "--output", required=True, help="Output fasta")
    args = parser.parse_args()

    remove_gaps(args.input, args.output)

if __name__ == "__main__":
    main()
