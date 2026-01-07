#!/usr/bin/env python3
"""
Summarize AMR and Virulence Results
Creates consolidated matrices from pipeline outputs
"""

import argparse
import pandas as pd
import sys
from pathlib import Path

def parse_amr_file(amr_file):
    """Parse AMRFinderPlus output file"""
    try:
        df = pd.read_csv(amr_file, sep='\t')
        if df.empty:
            return set()
        # Extract unique gene symbols
        genes = set(df['Gene symbol'].unique())
        return genes
    except Exception as e:
        print(f"Warning: Could not parse {amr_file}: {e}", file=sys.stderr)
        return set()

def parse_virulence_file(vir_file):
    """Parse virulence BLAST output"""
    try:
        df = pd.read_csv(vir_file, sep='\t', header=None,
                        names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                               'gapopen', 'qstart', 'qend', 'sstart', 'send',
                               'evalue', 'bitscore', 'stitle'])
        if df.empty:
            return set()
        # Filter for high quality hits
        df_filt = df[(df['pident'] >= 70) & (df['length'] >= 50)]
        # Extract virulence factor names
        vfs = set()
        for title in df_filt['stitle'].unique():
            # Clean up the title
            vf_name = title.split('(')[0].strip()
            vfs.add(vf_name)
        return vfs
    except Exception as e:
        print(f"Warning: Could not parse {vir_file}: {e}", file=sys.stderr)
        return set()

def create_presence_absence_matrix(sample_dict, output_file):
    """Create presence/absence matrix from dictionary of samples and features"""
    # Get all unique features
    all_features = sorted(set().union(*sample_dict.values()))
    
    # Create matrix
    matrix = []
    for sample, features in sorted(sample_dict.items()):
        row = [sample]
        for feature in all_features:
            row.append('1' if feature in features else '0')
        matrix.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(matrix, columns=['Sample'] + all_features)
    
    # Save to file
    df.to_csv(output_file, sep='\t', index=False)
    print(f"Created matrix: {output_file}")
    print(f"  Samples: {len(sample_dict)}")
    print(f"  Features: {len(all_features)}")

def main():
    parser = argparse.ArgumentParser(description='Summarize AMR and virulence results')
    parser.add_argument('--amr_dir', type=str, required=True,
                       help='Directory containing AMR results')
    parser.add_argument('--virulence_dir', type=str, required=True,
                       help='Directory containing virulence results')
    parser.add_argument('--output_dir', type=str, required=True,
                       help='Output directory for summary files')
    
    args = parser.parse_args()
    
    # Create output directory
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Process AMR files
    amr_dir = Path(args.amr_dir)
    amr_samples = {}
    for amr_file in amr_dir.glob('*_amr.tsv'):
        sample = amr_file.stem.replace('_amr', '')
        genes = parse_amr_file(amr_file)
        if genes:
            amr_samples[sample] = genes
    
    if amr_samples:
        create_presence_absence_matrix(
            amr_samples,
            output_dir / 'amr_matrix.tsv'
        )
    else:
        print("No AMR results found", file=sys.stderr)
    
    # Process virulence files
    vir_dir = Path(args.virulence_dir)
    vir_samples = {}
    for vir_file in vir_dir.glob('*_virulence.tsv'):
        sample = vir_file.stem.replace('_virulence', '')
        vfs = parse_virulence_file(vir_file)
        if vfs:
            vir_samples[sample] = vfs
    
    if vir_samples:
        create_presence_absence_matrix(
            vir_samples,
            output_dir / 'virulence_matrix.tsv'
        )
    else:
        print("No virulence results found", file=sys.stderr)

if __name__ == '__main__':
    main()
