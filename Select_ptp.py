import argparse
from collections import defaultdict

special_prefixes = ["BIOD", "CCCP", "GBDL", "HNSP", "MIZA", "QINL", "SPSO", "SRAA", "BGLP", "ZIPC"]

def parse_fasta(fasta_file):
    seq_lengths = {}
    seq_ids = {}
    with open(fasta_file, 'r') as f:
        seq_id = None
        seq = ''
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if seq_id:
                    seq_lengths[seq_id] = len(seq)
                    seq_ids[seq_id] = seq
                seq_id = line[1:]
                seq = ''
            else:
                seq += line
        if seq_id:
            seq_lengths[seq_id] = len(seq)
            seq_ids[seq_id] = seq
    return seq_lengths, seq_ids

def parse_mptp_output(mptp_file):
    species_groups = defaultdict(list)
    with open(mptp_file, 'r') as f:
        current_species = None
        for line in f:
            line = line.strip()
            if line.startswith('Species'):
                current_species = int(line.split()[1].strip(':'))
            elif line and not line.startswith('Species'):
                species_groups[current_species].extend(line.split())
    return species_groups

def select_sequences(species_groups, seq_lengths, seq_ids):
    selected_sequences = []
    
    for species, sequences in species_groups.items():
        longest_sequence = None
        longest_length = 0
        special_sequence = None
        longest_special_sequence = None
        longest_special_length = 0
        
        for seq in sequences:
            if seq in seq_lengths:
                seq_len = seq_lengths[seq]
                
                if any(seq.startswith(prefix) for prefix in special_prefixes):
                    if seq_len > longest_special_length:
                        longest_special_sequence = seq
                        longest_special_length = seq_len
                else:
                    if seq_len > longest_length:
                        longest_sequence = seq
                        longest_length = seq_len
        
        if longest_special_sequence:
            selected_sequences.append(longest_special_sequence)
        else:
            selected_sequences.append(longest_sequence)
    
    return selected_sequences

def write_selected_sequences(fasta_file, selected_sequences, seq_ids, output_file):
    with open(output_file, 'w') as f:
        for seq_id in selected_sequences:
            if seq_id in seq_ids:
                f.write(f">{seq_id}\n{seq_ids[seq_id]}\n")

def main():
    parser = argparse.ArgumentParser(description="Select longest sequence per species from a MPTP output and FASTA alignment.")
    
    parser.add_argument('-f', '--fasta', required=True, help="Input FASTA file containing the sequence alignment.")
    parser.add_argument('-m', '--mptp', required=True, help="Input MPTP output file with species groupings.")
    parser.add_argument('-o', '--output', required=True, help="Output FASTA file to save the selected sequences.")
    
    args = parser.parse_args()
    
    seq_lengths, seq_ids = parse_fasta(args.fasta)
    species_groups = parse_mptp_output(args.mptp)
    
    selected_sequences = select_sequences(species_groups, seq_lengths, seq_ids)
    
    write_selected_sequences(args.fasta, selected_sequences, seq_ids, args.output)
    
    print(f"Selected sequences saved to {args.output}")

if __name__ == "__main__":
    main()
