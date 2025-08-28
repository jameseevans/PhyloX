#!/usr/bin/env python3
"""
Mitochondrial Gene Presence/Absence Analysis Script
Analyzes GenBank files for the presence of 15 standard mitochondrial genes
and adds gene presence patterns to an existing CSV file for phylogenetic analysis.

Author: Generated for phylogenetic analysis
Usage: python gene_presence_analyzer.py -i input.gb -c existing_data.csv -o output.csv
"""

import argparse
import sys
import pandas as pd
from Bio import SeqIO
import re
from collections import defaultdict

def create_gene_lookup():
    """
    Create a comprehensive lookup dictionary for all gene name variants.
    Returns a dictionary mapping all possible gene names to their standard names.
    """
    gene_variants = {
        'ATP6': ['ATP6', 'ATP SYNTHASE F0 SUBUNIT 6', 'APT6', 'ATP SYNTHASE A0 SUBUNIT 6', 
                'ATP SYNTHASE SUBUNIT 6', 'ATP SYNTHASE FO SUBUNIT 6', 'ATPASE6', 'ATPASE SUBUNIT 6'],
        'ATP8': ['ATP8', 'ATP SYNTHASE F0 SUBUNIT 8', 'APT8', 'ATP SYNTHASE A0 SUBUNIT 8',
                'ATP SYNTHASE SUBUNIT 8', 'ATP SYNTHASE FO SUBUNIT 8', 'ATPASE8', 'ATPASE SUBUNIT 8'],
        'COX1': ['COX1', 'CYTOCHROME C OXIDASE SUBUNIT 1', 'CYTOCHROME OXIDASE SUBUNIT I',
                'CYTOCHROME C OXIDASE SUBUNIT I', 'COXI', 'CO1', 'COI', 'CYTOCHROME COXIDASE SUBUNIT I',
                'CYTOCHROME OXIDASE SUBUNIT 1', 'CYTOCHROME OXYDASE SUBUNIT 1'],
        'COX2': ['COX2', 'CYTOCHROME C OXIDASE SUBUNIT 2', 'CYTOCHROME OXIDASE SUBUNIT II',
                'CYTOCHROME C OXIDASE SUBUNIT II', 'COXII', 'CO2', 'COII', 'CYTOCHROME COXIDASE SUBUNIT II',
                'CYTOCHROME OXIDASE SUBUNIT 2', 'CYTOCHROME OXYDASE SUBUNIT 2'],
        'COX3': ['COX3', 'CYTOCHROME C OXIDASE SUBUNIT 3', 'CYTOCHROME OXIDASE SUBUNIT III',
                'CYTOCHROME C OXIDASE SUBUNIT III', 'COXIII', 'CO3', 'COIII', 'CYTOCHROME COXIDASE SUBUNIT III',
                'CYTOCHROME OXIDASE SUBUNIT 3', 'CYTOCHROME OXYDASE SUBUNIT 3'],
        'CYTB': ['CYTB', 'CYTOCHROME B', 'CYB', 'COB', 'COB / CYTB'],
        'ND1': ['ND1', 'NAD1', 'NSD1', 'NADH1', 'NADH DEHYDROGENASE SUBUNIT I',
               'NADH DEHYDROGENASE SUBUNIT 1', 'NADH DESHYDROGENASE SUBUNIT 1', 'NAD1-0'],
        'ND2': ['ND2', 'NAD2', 'NSD2', 'NADH2', 'NADH DEHYDROGENASE SUBUNIT II',
               'NADH DEHYDROGENASE SUBUNIT 2', 'NADH DESHYDROGENASE SUBUNIT 2', 'NAD2-0'],
        'ND3': ['ND3', 'NAD3', 'NSD3', 'NADH3', 'NADH DEHYDROGENASE SUBUNIT III',
               'NADH DEHYDROGENASE SUBUNIT 3', 'NADH DESHYDROGENASE SUBUNIT 3', 'NAD3-0'],
        'ND4': ['ND4', 'NAD4', 'NSD4', 'NADH4', 'NADH DEHYDROGENASE SUBUNIT IV',
               'NADH DEHYDROGENASE SUBUNIT 4', 'NADH DESHYDROGENASE SUBUNIT 4', 'NAD4-0'],
        'ND4L': ['ND4L', 'NAD4L', 'NSD4L', 'NADH4L', 'NADH DEHYDROGENASE SUBUNIT IVL',
                'NADH DEHYDROGENASE SUBUNIT 4L', 'NADH DESHYDROGENASE SUBUNIT 4L', 'NAD4L-0'],
        'ND5': ['ND5', 'NAD5', 'NSD5', 'NADH5', 'NADH DEHYDROGENASE SUBUNIT V',
               'NADH DEHYDROGENASE SUBUNIT 5', 'NADH DESHYDROGENASE SUBUNIT 5', 'NAD5-0'],
        'ND6': ['ND6', 'NAD6', 'NSD6', 'NADH6', 'NADH DEHYDROGENASE SUBUNIT VI',
               'NADH DEHYDROGENASE SUBUNIT 6', 'NADH DESHYDROGENASE SUBUNIT 6', 'NAD6-0'],
        'LSU': ['LSU', 'RRNL', 'L', 'LARGE', 'L-RRNA', '16S RIBOSOMAL RNA', '16SRRN',
               '16S-RRNA', 'RRN16S', '16S', '16S-RNA', '16S RRNA', 'RRNL RNA'],
        'SSU': ['SSU', 'RRNS', 'S', 'SMALL', 'S-RRNA', '12S RIBOSOMAL RNA', '12SRRN',
               '12S-RRNA', 'RRN12S', '12S', '12S-RNA', '12S RRNA', 'RRNS RNA']
    }
    
    # Create reverse lookup dictionary
    lookup = {}
    for standard_name, variants in gene_variants.items():
        for variant in variants:
            lookup[variant.upper()] = standard_name
    
    return lookup

def extract_gene_name(feature):
    """
    Extract gene name from a GenBank feature, trying multiple annotation fields.
    Returns the standardized gene name or None if not found.
    """
    gene_lookup = create_gene_lookup()
    
    # List of qualifiers to check for gene names
    qualifiers_to_check = ['gene', 'product', 'note', 'label']
    
    for qualifier in qualifiers_to_check:
        if qualifier in feature.qualifiers:
            for value in feature.qualifiers[qualifier]:
                # Clean and normalize the gene name
                clean_name = re.sub(r'[^\w\s/-]', '', value.upper().strip())
                clean_name = re.sub(r'\s+', ' ', clean_name)
                
                # Direct lookup
                if clean_name in gene_lookup:
                    return gene_lookup[clean_name]
                
                # Try partial matches for complex annotations
                for variant, standard in gene_lookup.items():
                    if variant in clean_name or clean_name in variant:
                        return standard
    
    return None

def analyze_genbank_sequence(record):
    """
    Analyze a single GenBank sequence record and return gene presence information.
    """
    # Define the gene order for the presence/absence string
    gene_order = ['ATP6', 'ATP8', 'COX1', 'COX2', 'COX3', 'CYTB', 'ND1', 'ND2', 
                  'ND3', 'ND4', 'ND4L', 'ND5', 'ND6', 'LSU', 'SSU']
    
    # Initialize gene presence dictionary
    genes_present = defaultdict(bool)
    
    # Check all features in the sequence
    for feature in record.features:
        if feature.type in ['gene', 'CDS', 'rRNA', 'tRNA']:
            gene_name = extract_gene_name(feature)
            if gene_name and gene_name in gene_order:
                genes_present[gene_name] = True
    
    # Create presence/absence pattern string
    pattern = ""
    for gene in gene_order:
        if genes_present[gene]:
            pattern += "X"
        else:
            pattern += "-"
    
    # Add hyphens at positions 4, 7, and 10 for readability
    formatted_pattern = (pattern[:4] + "-" + pattern[4:6] + "-" + 
                       pattern[6:8] + "-" + "-" + pattern[8:])
    
    return formatted_pattern

def process_genbank_file(genbank_file):
    """
    Process GenBank file and return dictionary mapping sequence IDs to gene presence patterns.
    """
    gene_patterns = {}
    processed_count = 0
    
    try:
        with open(genbank_file, 'r') as handle:
            for record in SeqIO.parse(handle, "genbank"):
                pattern = analyze_genbank_sequence(record)
                gene_patterns[record.id] = pattern
                processed_count += 1
                
        print(f"Processed {processed_count} sequences from GenBank file")
        return gene_patterns
        
    except Exception as e:
        print(f"Error processing GenBank file '{genbank_file}': {e}", file=sys.stderr)
        sys.exit(1)

def load_csv_file(csv_file):
    """
    Load existing CSV file and determine the ID column.
    """
    try:
        df = pd.read_csv(csv_file)
        
        # Determine which column contains sequence IDs
        id_column = None
        if 'mt_id' in df.columns:
            id_column = 'mt_id'
            print(f"Using 'mt_id' column for sequence identifiers")
        else:
            id_column = df.columns[0]
            print(f"Using first column '{id_column}' for sequence identifiers")
        
        return df, id_column
        
    except Exception as e:
        print(f"Error loading CSV file '{csv_file}': {e}", file=sys.stderr)
        sys.exit(1)

def update_csv_with_gene_presence(df, id_column, gene_patterns):
    """
    Add gene_presence column to the DataFrame.
    """
    # Initialize gene_presence column
    if 'gene_presence' not in df.columns:
        df['gene_presence'] = None
    
    matches_found = 0
    no_matches = []
    
    # Update gene presence for matching sequences
    for idx, row in df.iterrows():
        seq_id = str(row[id_column])
        
        # Try exact match first
        if seq_id in gene_patterns:
            df.at[idx, 'gene_presence'] = gene_patterns[seq_id]
            matches_found += 1
        else:
            # Try to find partial matches (in case of version numbers or prefixes)
            found_match = False
            for gb_id in gene_patterns:
                # Remove version numbers and compare
                clean_seq_id = seq_id.split('.')[0]
                clean_gb_id = gb_id.split('.')[0]
                
                if clean_seq_id == clean_gb_id or clean_seq_id in clean_gb_id or clean_gb_id in clean_seq_id:
                    df.at[idx, 'gene_presence'] = gene_patterns[gb_id]
                    matches_found += 1
                    found_match = True
                    break
            
            if not found_match:
                no_matches.append(seq_id)
                df.at[idx, 'gene_presence'] = "---------------"  # All genes absent for unmatched sequences
    
    print(f"Matched {matches_found} sequences with GenBank data")
    if no_matches:
        print(f"No GenBank matches found for {len(no_matches)} sequences:")
        for seq_id in no_matches[:10]:  # Show first 10
            print(f"  - {seq_id}")
        if len(no_matches) > 10:
            print(f"  ... and {len(no_matches) - 10} more")
    
    return df

def print_gene_presence_summary(df):
    """
    Print summary statistics about gene presence patterns.
    """
    if 'gene_presence' not in df.columns:
        return
    
    # Count different patterns
    pattern_counts = df['gene_presence'].value_counts()
    
    # Calculate completeness statistics
    complete_genomes = len(df[df['gene_presence'] == 'XXXX-XX-XX--XX-'])
    total_sequences = len(df[df['gene_presence'].notna()])
    
    print(f"\nGene Presence Summary:")
    print(f"Total sequences in CSV: {len(df)}")
    print(f"Sequences with gene data: {total_sequences}")
    print(f"Complete mitochondrial genomes (15/15 genes): {complete_genomes}")
    print(f"Incomplete genomes: {total_sequences - complete_genomes}")
    
    if len(pattern_counts) <= 20:  # Show all patterns if not too many
        print(f"\nGene presence patterns:")
        for pattern, count in pattern_counts.head(10).items():
            print(f"  {pattern}: {count} sequences")
        if len(pattern_counts) > 10:
            print(f"  ... and {len(pattern_counts) - 10} other patterns")

def main():
    parser = argparse.ArgumentParser(
        description='Add mitochondrial gene presence/absence patterns to existing CSV file',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
This script processes a GenBank file to extract gene presence information
and adds it to an existing CSV file in a new 'gene_presence' column.

Gene Order and Pattern Format:
ATP6-ATP8-COX1-COX2-COX3-CYTB-ND1-ND2-ND3-ND4-ND4L-ND5-ND6-LSU-SSU

Pattern Example: XXXX-XX-XX--XX-
Where X = gene present, - = gene absent

The script looks for sequence IDs in the 'mt_id' column, or uses the first column if 'mt_id' is not found.
        """
    )
    
    parser.add_argument('-i', '--input', '--genbank', required=True,
                       help='Input GenBank file (.gb or .gbk)')
    parser.add_argument('-c', '--csv', required=True,
                       help='Input CSV file with existing data')
    parser.add_argument('-o', '--output', required=True,
                       help='Output CSV file with added gene_presence column')
    parser.add_argument('--summary', action='store_true',
                       help='Print detailed summary statistics')
    
    args = parser.parse_args()
    
    # Process GenBank file
    print(f"Processing GenBank file: {args.input}")
    gene_patterns = process_genbank_file(args.input)
    
    if not gene_patterns:
        print("No sequences found in GenBank file.", file=sys.stderr)
        sys.exit(1)
    
    # Load existing CSV
    print(f"Loading CSV file: {args.csv}")
    df, id_column = load_csv_file(args.csv)
    
    # Update CSV with gene presence data
    print(f"Updating gene presence information...")
    df_updated = update_csv_with_gene_presence(df, id_column, gene_patterns)
    
    # Save updated CSV
    try:
        df_updated.to_csv(args.output, index=False)
        print(f"Updated CSV saved to: {args.output}")
    except Exception as e:
        print(f"Error saving output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Print summary
    if args.summary or len(gene_patterns) > 0:
        print_gene_presence_summary(df_updated)
    
    print(f"\nAnalysis complete. Updated {len(df)} CSV records with gene presence data.")

if __name__ == "__main__":
    main()
