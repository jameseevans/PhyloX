### Remove sequences from a Genbank file that contian an identical COX1 region, based on total sequence length ###

import argparse
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

cox1_synonyms = [
    "COX1", "Cytochrome c oxidase subunit 1", "CYTOCHROME C OXIDASE SUBUNIT 1",
    "CYTOCHROME OXIDASE SUBUNIT I", "CYTOCHROME C OXIDASE SUBUNIT I", "COXI",
    "CO1", "COI", "CYTOCHROME COXIDASE SUBUNIT I", "CYTOCHROME OXIDASE SUBUNIT 1",
    "CYTOCHROME OXYDASE SUBUNIT 1"
]

def extract_cox1_sequence(record):
    for feature in record.features:
        if feature.type == "CDS":
            if "gene" in feature.qualifiers:
                gene_names = feature.qualifiers["gene"]
                if any(name.upper() in cox1_synonyms for name in gene_names):
                    return feature.extract(record.seq)
            if "product" in feature.qualifiers:
                product_names = feature.qualifiers["product"]
                if any(name.upper() in cox1_synonyms for name in product_names):
                    return feature.extract(record.seq)
    return None

def main(input_file, output_file):
    cox1_dict = {}
    removed_records = []
    for record in SeqIO.parse(input_file, "genbank"):
        cox1_seq = extract_cox1_sequence(record)
        if cox1_seq:
            cox1_seq_str = str(cox1_seq)
            if cox1_seq_str not in cox1_dict:
                cox1_dict[cox1_seq_str] = record
            else:
                if len(record.seq) > len(cox1_dict[cox1_seq_str].seq):
                    removed_records.append(cox1_dict[cox1_seq_str].name)
                    cox1_dict[cox1_seq_str] = record
                else:
                    removed_records.append(record.name)
    filtered_records = list(cox1_dict.values())
    SeqIO.write(filtered_records, output_file, "genbank")
    with open("removed.txt", "w") as removed_file:
        for record_name in removed_records:
            removed_file.write(record_name + "\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter GenBank sequences to remove identical COX1 sequences based on total length")
    parser.add_argument("-i", "--input", required=True, help="Input GenBank file")
    parser.add_argument("-o", "--output", required=True, help="Output GenBank file")

    args = parser.parse_args()
    main(args.input, args.output)
