### Script to download barcode sequences and associated metadata ("specimen") from BOLD (https://www.boldsystems.org/) ###

import argparse
import requests
import csv

def download_data(taxon):
    sequence_url = f"http://www.boldsystems.org/index.php/API_Public/sequence?taxon={taxon}"
    specimen_url = f"http://www.boldsystems.org/index.php/API_Public/specimen?taxon={taxon}&format=tsv"
    print(f"Downloading barcodes for taxon: {taxon}")
    sequences = requests.get(sequence_url).text
    with open(f"{taxon}_barcodes.fasta", "w") as file:
        file.write(sequences)
    print(f"Sequences saved to {taxon}_barcodes.fasta")
    print(f"Downloading metadata for taxon: {taxon}")
    metadata_response = requests.get(specimen_url)
    metadata_tsv = metadata_response.text
    temp = f"{taxon}_metadata.tsv"
    with open(temp, "w") as file:
        file.write(metadata_tsv)
    metadata_csv = f"{taxon}_metadata.csv"
    with open(temp, "r") as tsvfile, open(metadata_csv, "w", newline='') as csvfile:
        tsv_reader = csv.reader(tsvfile, delimiter='\t')
        csv_writer = csv.writer(csvfile)
        for row in tsv_reader:
            csv_writer.writerow(row)
    print(f"Metadata saved to {taxon}_metadata.csv")
    print("Download complete")

def main():
    parser = argparse.ArgumentParser(description="Download barcode sequences and associated metadata from BOLD")
    parser.add_argument("-t", "--taxon", required=True, help="The taxon for which to download data")
    args = parser.parse_args()
    download_data(args.taxon)

if __name__ == "__main__":
    main()
