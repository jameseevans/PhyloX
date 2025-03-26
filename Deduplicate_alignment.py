### Script to remove duplicate sequences from a fasta alignment file, retaining mitochondrial genomes from the Vogler lab, or the longest sequence ###
### james.evans@nhm.ac.uk ###

import argparse
from Bio import SeqIO

mito_prefixes = {'BIOD', 'CCCP', 'GBDL', 'HNSP', 'MIZA', 'QINL', 'SPSO', 'SRAA', 'BGLP'}

def get_longest_sequence(sequences):
    mito_sequences = [seq for seq in sequences if any(seq.id.startswith(prefix) for prefix in mito_prefixes)]
    if mito_sequences:
        return max(mito_sequences, key=lambda seq: len(seq.seq))
    else:
        return max(sequences, key=lambda seq: len(seq.seq))

def remove_duplicates(input_fasta, output_fasta):
    sequences = SeqIO.parse(input_fasta, "fasta")
    unique_sequences = {}
    for seq_record in sequences:
        seq_str = str(seq_record.seq)
        if seq_str not in unique_sequences:
            unique_sequences[seq_str] = [seq_record]
        else:
            unique_sequences[seq_str].append(seq_record)
    final_sequences = [get_longest_sequence(seq_records) for seq_records in unique_sequences.values()]
    with open(output_fasta, "w") as out_file:
        SeqIO.write(final_sequences, out_file, "fasta")
    print(f"Filtered alignment has been written to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a fasta alignment, retaining mitogenomes or the longest sequences.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta")
    parser.add_argument("-o", "--output", required=True, help="Output fasta")

    args = parser.parse_args()
    
    remove_duplicates(args.input, args.output)
