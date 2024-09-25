### Script to remove sequences from genbank files that do not contain a COX1 annotation ###

synonyms = set([
    "COX1", "Cytochrome c oxidase subunit 1", "CYTOCHROME C OXIDASE SUBUNIT 1",
    "CYTOCHROME OXIDASE SUBUNIT I", "COXI", "CO1", "COI",
    "CYTOCHROME COXIDASE SUBUNIT I", "CYTOCHROME OXIDASE SUBUNIT 1", "CYTOCHROME OXYDASE SUBUNIT 1"
])

def parse_arguments():
    parser = argparse.ArgumentParser(description="Filter GenBank sequences and metadata based on COX1 presence.")
    parser.add_argument('-g', '--input-genbank', required=True, help="Path to the input GenBank file.")
    parser.add_argument('-m', '--input-metadata', required=True, help="Path to the input metadata csv file.")
    parser.add_argument('-og', '--output-genbank', required=True, help="Path to the output GenBank file.")
    parser.add_argument('-om', '--output-metadata', required=True, help="Path to the output metadata CSV file.")
    return parser.parse_args()

def filter_genbank_file(input_genbank, output_genbank, synonyms):
    sequences_with_cox1 = set()
    for record in SeqIO.parse(input_genbank, "genbank"):
        locus_id = record.name
        if any(feature.qualifiers.get('product', [''])[0] in synonyms for feature in record.features):
            sequences_with_cox1.add(locus_id)
    with open(input_genbank, "r") as input_handle, open(output_genbank, "w") as output_handle:
        for record in SeqIO.parse(input_handle, "genbank"):
            locus_id = record.name
            if locus_id in sequences_with_cox1:
                SeqIO.write(record, output_handle, "genbank")
    return sequences_with_cox1

def filter_metadata_csv(input_metadata, output_metadata, sequences_with_cox1):
    df = pd.read_csv(input_metadata)
    df['db_id'] = df['db_id'].str.strip().str.upper()
    sequences_with_cox1 = set(seq.upper() for seq in sequences_with_cox1)
    filtered_df = df[df['db_id'].isin(sequences_with_cox1)]
    filtered_df.to_csv(output_metadata, index=False)
  
def main():
    args = parse_arguments()
    sequences_with_cox1 = filter_genbank_file(args.input_genbank, args.output_genbank, synonyms)
    filter_metadata_csv(args.input_metadata, args.output_metadata, sequences_with_cox1)

if __name__ == "__main__":
    main()
