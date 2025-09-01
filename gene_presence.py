#!/usr/bin/env python3
import argparse
import csv
import re
from Bio import SeqIO

GENES = [
    ("ATP6", ["ATP6", "ATP SYNTHASE F0 SUBUNIT 6", "APT6", "ATP SYNTHASE A0 SUBUNIT 6",
              "ATP SYNTHASE SUBUNIT 6", "ATP SYNTHASE FO SUBUNIT 6", "ATPASE6", "ATPASE SUBUNIT 6"]),
    ("ATP8", ["ATP8", "ATP SYNTHASE F0 SUBUNIT 8", "APT8", "ATP SYNTHASE A0 SUBUNIT 8",
              "ATP SYNTHASE SUBUNIT 8", "ATP SYNTHASE FO SUBUNIT 8", "ATPASE8", "ATPASE SUBUNIT 8"]),
    ("COX1", ["COX1", "CYTOCHROME C OXIDASE SUBUNIT 1", "CYTOCHROME OXIDASE SUBUNIT I",
              "CYTOCHROME C OXIDASE SUBUNIT I", "COXI", "CO1", "COI",
              "CYTOCHROME COXIDASE SUBUNIT I", "CYTOCHROME OXIDASE SUBUNIT 1",
              "CYTOCHROME OXYDASE SUBUNIT 1"]),
    ("COX2", ["COX2", "CYTOCHROME C OXIDASE SUBUNIT 2", "CYTOCHROME OXIDASE SUBUNIT II",
              "CYTOCHROME C OXIDASE SUBUNIT II", "COXII", "CO2", "COII",
              "CYTOCHROME COXIDASE SUBUNIT II", "CYTOCHROME OXIDASE SUBUNIT 2",
              "CYTOCHROME OXYDASE SUBUNIT 2"]),
    ("COX3", ["COX3", "CYTOCHROME C OXIDASE SUBUNIT 3", "CYTOCHROME OXIDASE SUBUNIT III",
              "CYTOCHROME C OXIDASE SUBUNIT III", "COXIII", "CO3", "COIII",
              "CYTOCHROME COXIDASE SUBUNIT III", "CYTOCHROME OXIDASE SUBUNIT 3",
              "CYTOCHROME OXYDASE SUBUNIT 3"]),
    ("CYTB", ["CYTB", "CYTOCHROME B", "CYB", "COB", "COB / CYTB"]),
    ("ND1", ["ND1", "NAD1", "NSD1", "NADH1", "NADH DEHYDROGENASE SUBUNIT I",
             "NADH DEHYDROGENASE SUBUNIT 1", "NADH DESHYDROGENASE SUBUNIT 1", "NAD1-0"]),
    ("ND2", ["ND2", "NAD2", "NSD2", "NADH2", "NADH DEHYDROGENASE SUBUNIT II",
             "NADH DEHYDROGENASE SUBUNIT 2", "NADH DESHYDROGENASE SUBUNIT 2", "NAD2-0"]),
    ("ND3", ["ND3", "NAD3", "NSD3", "NADH3", "NADH DEHYDROGENASE SUBUNIT III",
             "NADH DEHYDROGENASE SUBUNIT 3", "NADH DESHYDROGENASE SUBUNIT 3", "NAD3-0"]),
    ("ND4", ["ND4", "NAD4", "NSD4", "NADH4", "NADH DEHYDROGENASE SUBUNIT IV",
             "NADH DEHYDROGENASE SUBUNIT 4", "NADH DESHYDROGENASE SUBUNIT 4", "NAD4-0"]),
    ("ND4L", ["ND4L", "NAD4L", "NSD4L", "NADH4L", "NADH DEHYDROGENASE SUBUNIT IVL",
              "NADH DEHYDROGENASE SUBUNIT 4L", "NADH DESHYDROGENASE SUBUNIT 4L", "NAD4L-0"]),
    ("ND5", ["ND5", "NAD5", "NSD5", "NADH5", "NADH DEHYDROGENASE SUBUNIT V",
             "NADH DEHYDROGENASE SUBUNIT 5", "NADH DESHYDROGENASE SUBUNIT 5", "NAD5-0"]),
    ("ND6", ["ND6", "NAD6", "NSD6", "NADH6", "NADH DEHYDROGENASE SUBUNIT VI",
             "NADH DEHYDROGENASE SUBUNIT 6", "NADH DESHYDROGENASE SUBUNIT 6", "NAD6-0"]),
    ("LSU", ["LSU", "RRNL", "L", "LARGE", "L-RRNA", "16S RIBOSOMAL RNA", "16SRRN", "16S-RRNA",
             "RRN16S", "16S", "16S-RNA", "16S RRNA", "RRNL RNA"]),
    ("SSU", ["SSU", "RRNS", "S", "SMALL", "S-RRNA", "12S RIBOSOMAL RNA", "12SRRN", "12S-RRNA",
             "RRN12S", "12S", "12S-RNA", "12S RRNA", "RRNS RNA"]),
]

PCG_COUNT = 13


def _normalize(name: str) -> str:
    """Normalize gene name string."""
    return re.sub(r"[^A-Z0-9]", "", name.upper())

NORM_GENES = [(gene, {_normalize(s) for s in synonyms}) for gene, synonyms in GENES]


def extract_gene_presence(record, include_rrna=False):
    """Return gene presence string (X/-) for a SeqRecord."""
    found = set()
    use_genes = NORM_GENES if include_rrna else NORM_GENES[:PCG_COUNT]

    for feature in record.features:
        if feature.type not in {"gene", "CDS", "rRNA"}:
            continue
        for vals in feature.qualifiers.values():
            for raw in vals:
                nval = _normalize(raw)
                for gene, nsyns in use_genes:
                    if gene in found:
                        continue
                    if nval in nsyns:
                        found.add(gene)
                        break

    return "".join("X" if gene in found else "-" for gene, _ in use_genes)


def main():
    parser = argparse.ArgumentParser(
        description="Count gene presence/absence from GenBank into CSV"
    )
    parser.add_argument("-g", "--genbank", required=True, help="Input GenBank file (.gb)")
    parser.add_argument("-c", "--csv", required=True, help="Input CSV file (sequence names in first column)")
    parser.add_argument("-o", "--output", required=True, help="Output CSV file")
    parser.add_argument("--include-rrna", action="store_true", help="Include rRNA genes (15 total instead of 13)")
    parser.add_argument("--summary-out", help="Optional summary stats output file")
    args = parser.parse_args()

    gb_records = {rec.name: rec for rec in SeqIO.parse(args.genbank, "genbank")}

    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        rows = list(reader)
        fieldnames = reader.fieldnames or []

    if "gene_presence" not in fieldnames:
        fieldnames.append("gene_presence")

    updated_rows = []
    gene_count = len(GENES) if args.include_rrna else PCG_COUNT

    for row in rows:
        seqname = row[fieldnames[0]]
        if seqname in gb_records:
            presence = extract_gene_presence(gb_records[seqname], args.include_rrna)
            row["gene_presence"] = presence
        else:
            row["gene_presence"] = "-" * gene_count
        updated_rows.append(row)

    with open(args.output, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(updated_rows)

    total = len(updated_rows)
    gene_count = len(GENES) if args.include_rrna else PCG_COUNT
    complete = sum(1 for r in updated_rows if r["gene_presence"].count("X") == gene_count)
    incomplete = total - complete

    summary_lines = [
        "Gene Presence Summary:",
        f"  Total sequences: {total}",
        f"  Complete sequences ({gene_count}/{gene_count} genes): {complete}",
        f"  Incomplete sequences: {incomplete}",
    ]
    summary = "\n".join(summary_lines)
    print("\n" + summary)

    if args.summary_out:
        with open(args.summary_out, "w") as fh:
            fh.write(summary + "\n")


if __name__ == "__main__":
    main()
