import argparse
from Bio import SeqIO

def filter_sequences(input_fasta, output_fasta, min_length=500):
    with open(input_fasta, "r") as infile:
        with open(output_fasta, "w") as outfile:
            for record in SeqIO.parse(infile, "fasta"):
                non_gap_bases = str(record.seq).replace('-', '')
                if len(non_gap_bases) >= min_length:
                    SeqIO.write(record, outfile, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Filter sequences based on a minimum length.")
    parser.add_argument('-i', '--input', required=True, help="Input FASTA file")
    parser.add_argument('-o', '--output', required=True, help="Output FASTA file")
    parser.add_argument('-l', '--min-length', type=int, default=500, help="Minimum length (default: 500)")

    args = parser.parse_args()
    filter_sequences(args.input, args.output, args.min_length)

if __name__ == "__main__":
    main()
