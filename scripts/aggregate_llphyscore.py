#!/usr/bin/env python3
"""
FASTA headerからprotein_id→gene_id/gene_symbolマッピングを抽出し、
LLPhyScoreスコアと統合して遺伝子レベルに集約

入力: 
  - EuRBPDB FASTA (*.RBP.fa) 
  - llphyscore_summary.tsv
出力:
  - llphyscore_per_gene.tsv (gene-level MAX scores)

使い方:
  python3 aggregate_llphyscore.py /Users/kyotayasuda/rbp_pfam/
"""

import sys
import os
import re

def parse_fasta_headers(fasta_dir):
    """FASTAヘッダーからprotein_id→gene情報を抽出"""
    mapping = []
    
    species_files = {
        'Caenorhabditis_elegans': 'Caenorhabditis_elegans.RBP.fa',
        'Drosophila_melanogaster': 'Drosophila_melanogaster.RBP.fa',
        'Danio_rerio': 'Danio_rerio.RBP.fa',
        'Xenopus_tropicalis': 'Xenopus_tropicalis.RBP.fa',
        'Mus_musculus': 'Mus_musculus.RBP.fa',
        'Homo_sapiens': 'Homo_sapiens.RBP.fa',
    }
    
    for species, filename in species_files.items():
        filepath = os.path.join(fasta_dir, filename)
        if not os.path.exists(filepath):
            print(f"WARNING: {filepath} not found")
            continue
        
        count = 0
        with open(filepath) as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]  # Remove '>'
                    
                    # Extract protein_id (first field, remove version)
                    parts = header.split()
                    protein_id_full = parts[0]
                    
                    # Parse key-value pairs from header
                    gene_id = ''
                    gene_symbol = ''
                    
                    for part in parts:
                        if part.startswith('gene:'):
                            gene_id = part.replace('gene:', '')
                        elif part.startswith('gene_symbol:'):
                            gene_symbol = part.replace('gene_symbol:', '')
                    
                    # If no gene: tag found, try other patterns
                    if not gene_id:
                        # Some FASTA headers might have different formats
                        # Try to find ENSG/WBGene/FBgn pattern
                        for part in parts:
                            if 'ENSG' in part or 'WBGene' in part or 'FBgn' in part:
                                gene_id = part
                                break
                    
                    mapping.append({
                        'protein_id': protein_id_full,
                        'gene_id': gene_id,
                        'gene_symbol': gene_symbol,
                        'species': species,
                    })
                    count += 1
        
        print(f"  {species}: {count} protein headers parsed")
    
    return mapping


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 aggregate_llphyscore.py <rbp_pfam_dir>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    output_dir = os.path.join(base_dir, 'output')
    
    # Step 1: Parse FASTA headers
    print("=== Step 1: Parsing FASTA headers ===")
    mapping = parse_fasta_headers(base_dir)
    print(f"Total mappings: {len(mapping)}")
    
    # Show first few
    print("\nSample mappings:")
    for m in mapping[:5]:
        print(f"  {m['protein_id']} -> {m['gene_id']} ({m['gene_symbol']})")
    
    # Step 2: Load LLPhyScore results
    print("\n=== Step 2: Loading LLPhyScore results ===")
    import pandas as pd
    
    scores = pd.read_csv(os.path.join(output_dir, 'llphyscore_summary.tsv'), sep='\t')
    scores['raw_score'] = pd.to_numeric(scores['raw_score'], errors='coerce')
    scores['percentile'] = pd.to_numeric(scores['percentile'], errors='coerce')
    print(f"Total scores: {len(scores)}")
    
    # Step 3: Create mapping DataFrame
    map_df = pd.DataFrame(mapping)
    
    # Strip version from protein_id in both dataframes for matching
    scores['protein_id_clean'] = scores['protein_id'].str.split('.').str[0]
    map_df['protein_id_clean'] = map_df['protein_id'].str.split('.').str[0]
    
    print(f"\nMapping df: {len(map_df)}")
    print(f"Scores df: {len(scores)}")
    
    # Step 4: Merge
    merged = scores.merge(map_df[['protein_id_clean', 'gene_id', 'gene_symbol', 'species']], 
                          on=['protein_id_clean'], how='left', suffixes=('', '_map'))
    
    # Use species from scores if mapping species is missing
    matched = merged['gene_id'].notna() & (merged['gene_id'] != '')
    print(f"Matched: {matched.sum()}/{len(merged)} ({matched.sum()/len(merged)*100:.1f}%)")
    
    # Show unmatched
    if (~matched).sum() > 0:
        unmatched_species = merged.loc[~matched, 'species'].value_counts()
        print(f"Unmatched by species:\n{unmatched_species}")
    
    # Step 5: Aggregate per gene (MAX score)
    print("\n=== Step 3: Aggregating per gene (MAX raw score) ===")
    
    # For matched entries
    matched_df = merged[matched].copy()
    
    gene_agg = matched_df.groupby(['gene_id', 'gene_symbol', 'species']).agg(
        max_raw_score=('raw_score', 'max'),
        mean_raw_score=('raw_score', 'mean'),
        max_percentile=('percentile', 'max'),
        mean_percentile=('percentile', 'mean'),
        n_isoforms=('protein_id', 'count'),
    ).reset_index()
    
    print(f"Gene-level entries: {len(gene_agg)}")
    
    for sp in ['Caenorhabditis_elegans', 'Drosophila_melanogaster', 'Danio_rerio',
               'Xenopus_tropicalis', 'Mus_musculus', 'Homo_sapiens']:
        n = len(gene_agg[gene_agg['species'] == sp])
        print(f"  {sp}: {n} genes")
    
    # Save
    out_file = os.path.join(output_dir, 'llphyscore_per_gene.tsv')
    gene_agg.to_csv(out_file, sep='\t', index=False)
    print(f"\nSaved: {out_file}")
    
    # Quick correlation check
    from scipy import stats as st
    
    neuron_map = {
        'Caenorhabditis_elegans': 302,
        'Drosophila_melanogaster': 100000,
        'Danio_rerio': 10000000,
        'Xenopus_tropicalis': 160000000,
        'Mus_musculus': 71000000,
        'Homo_sapiens': 86000000000,
    }
    species_order = ['Caenorhabditis_elegans', 'Drosophila_melanogaster', 
                     'Danio_rerio', 'Xenopus_tropicalis', 'Mus_musculus', 'Homo_sapiens']
    
    print("\n=== Species-level summary ===")
    rows = []
    for sp in species_order:
        sub = gene_agg[gene_agg['species'] == sp]
        high5 = (sub['max_raw_score'] > 5).sum() / len(sub) * 100
        high0 = (sub['max_raw_score'] > 0).sum() / len(sub) * 100
        rows.append({
            'species': sp, 'neurons': neuron_map[sp],
            'n_genes': len(sub),
            'mean_max_raw': sub['max_raw_score'].mean(),
            'high_pct_raw5': high5,
            'high_pct_raw0': high0,
        })
        print(f"  {sp}: {len(sub)} genes, mean_max={sub['max_raw_score'].mean():.1f}, "
              f"high(>5)={high5:.1f}%, high(>0)={high0:.1f}%")
    
    res = pd.DataFrame(rows)
    print("\n=== Correlations vs neurons ===")
    for col, label in [
        ('mean_max_raw', 'Mean max raw score'),
        ('high_pct_raw5', 'High-LLPS % (raw>5)'),
        ('high_pct_raw0', 'High-LLPS % (raw>0)'),
    ]:
        rho, p = st.spearmanr(res['neurons'], res[col])
        sig = '★' if p < 0.05 else 'ns'
        print(f"  {label:35s}: ρ={rho:.3f}, p={p:.4f} {sig}")


if __name__ == '__main__':
    main()
