import argparse
from Bio import SeqIO

special_prefixes = {'BIOD', 'CCCP', 'GBDL', 'HNSP', 'MIZA', 'QINL', 'SPSO', 'SRAA', 'BGLP'}

def get_longest_sequence(sequences):
    """Return the longest sequence, prioritizing those with special prefixes if present."""
    special_sequences = [seq for seq in sequences if any(seq.id.startswith(prefix) for prefix in special_prefixes)]
    
    if special_sequences:
        return max(special_sequences, key=lambda seq: len(seq.seq))
    else:
        return max(sequences, key=lambda seq: len(seq.seq))

def remove_duplicates(input_fasta, output_fasta):
    """Reads a FASTA file, removes duplicates, and writes the longest appropriate sequences to a new file."""
    sequences = SeqIO.parse(input_fasta, "fasta")
    
    unique_sequences = {}
    
    for seq_record in sequences:
        seq_name = seq_record.id.split()[0] 

        if seq_name not in unique_sequences:
            unique_sequences[seq_name] = [seq_record]
        else:
            unique_sequences[seq_name].append(seq_record)

    final_sequences = [get_longest_sequence(seq_records) for seq_records in unique_sequences.values()]
    
    with open(output_fasta, "w") as out_file:
        SeqIO.write(final_sequences, out_file, "fasta")

    print(f"Filtered sequences have been written to {output_fasta}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Remove duplicate sequences from a FASTA file while retaining the longest.")
    parser.add_argument("-i", "--input", required=True, help="Input FASTA file")
    parser.add_argument("-o", "--output", required=True, help="Output FASTA file")

    args = parser.parse_args()
    
    remove_duplicates(args.input, args.output)
