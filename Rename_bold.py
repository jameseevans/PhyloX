import argparse
from Bio import SeqIO

def rename_fasta_sequences(input_fasta, output_fasta):
    """
    Rename sequences in a FASTA file by keeping only the first part of the sequence ID (before '|').
    
    Args:
        input_fasta (str): Path to the input FASTA file.
        output_fasta (str): Path to the output FASTA file with renamed sequences.
    """
    with open(output_fasta, "w") as output_handle:
        for record in SeqIO.parse(input_fasta, "fasta"):

            new_id = record.id.split("|")[0]
            record.id = new_id
            record.description = ""
            SeqIO.write(record, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="Rename FASTA sequences to keep only the first part of the sequence ID (before '|').")
    parser.add_argument("-i", "--input", required=True, help="Path to the input FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Path to the output FASTA file.")
    args = parser.parse_args()

    rename_fasta_sequences(args.input, args.output)

if __name__ == "__main__":
    main()
