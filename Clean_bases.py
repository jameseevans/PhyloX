import argparse

parser = argparse.ArgumentParser(description="Clean sequences by removing gaps and replacing non-standard bases.")
parser.add_argument('-i', '--input', type=str, required=True, help="Input FASTA file")
parser.add_argument('-o', '--output', type=str, required=True, help="Output FASTA file")

args = parser.parse_args()

valid_bases = {'A', 'T', 'C', 'G', 'N', 'U', 'R', 'Y', 'S', 'W', 'K', 'M', 'B', 'D', 'H', 'V', '-'}

with open(args.input, 'r') as infile, open(args.output, 'w') as outfile:
    for line in infile:
        if line.startswith('>'):
            outfile.write(line)
        else:
            cleaned_sequence = ''.join(
                char if char.upper() in valid_bases else 'N' for char in line.strip()
            )
            outfile.write(cleaned_sequence + '\n')
