import argparse

parser = argparse.ArgumentParser(description="Filter FASTA sequences by header containing 'COI-5P'")
parser.add_argument('-i', '--input', type=str, required=True, help="Input FASTA file")
parser.add_argument('-o', '--output', type=str, required=True, help="Output FASTA file")

args = parser.parse_args()

removed_count = 0

with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
    write_sequence = False
    for line in infile:
        if line.startswith('>'):
            if 'COI-5P' in line:
                write_sequence = True
                outfile.write(line)
            else:
                write_sequence = False
                removed_count += 1
        elif write_sequence:
            outfile.write(line)

print(f"Number of sequences removed: {removed_count}")
