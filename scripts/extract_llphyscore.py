#!/usr/bin/env python3
"""
LLPhyScore出力からsum scoreだけ抽出

入力: llphyscore_*.txt (種ごと)
出力: llphyscore_summary.tsv (1ファイルに統合)

使い方:
  python3 extract_llphyscore.py /Users/kyotayasuda/rbp_pfam/output/
"""

import sys
import os
import re

def extract_scores(filepath, species):
    """LLPhyScore出力ファイルからprotein IDとsum scoreを抽出"""
    results = []
    current_protein = None
    
    with open(filepath, 'r') as f:
        for line in f:
            line = line.rstrip()
            if not line:
                continue
            
            parts = line.split()
            if len(parts) < 2:
                continue
            
            # New protein ID line (starts with non-space)
            if not line.startswith(' ') and not line.startswith('\t'):
                # This line has protein ID as first field
                current_protein = parts[0]
            
            # Check for "8-feature sum" line
            if '8-feature sum' in line:
                # Format: [protein_id?]  8-feature sum  raw  percentile  zscore  modified_zscore
                # Find the numbers after "8-feature sum"
                idx = line.index('8-feature sum')
                after = line[idx + len('8-feature sum'):].split()
                if len(after) >= 4:
                    raw_score = after[0]
                    percentile = after[1]
                    zscore = after[2]
                    mod_zscore = after[3]
                    
                    if current_protein:
                        results.append({
                            'protein_id': current_protein,
                            'species': species,
                            'raw_score': raw_score,
                            'percentile': percentile,
                            'zscore': zscore,
                            'mod_zscore': mod_zscore,
                        })
    
    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 extract_llphyscore.py <output_dir>")
        sys.exit(1)
    
    out_dir = sys.argv[1]
    
    species_files = {
        'Caenorhabditis_elegans': 'llphyscore_Caenorhabditis_elegans.txt',
        'Drosophila_melanogaster': 'llphyscore_Drosophila_melanogaster.txt',
        'Danio_rerio': 'llphyscore_Danio_rerio.txt',
        'Xenopus_tropicalis': 'llphyscore_Xenopus_tropicalis.txt',
        'Mus_musculus': 'llphyscore_Mus_musculus.txt',
        'Homo_sapiens': 'llphyscore_Homo_sapiens.txt',
    }
    
    all_results = []
    
    for species, filename in species_files.items():
        filepath = os.path.join(out_dir, filename)
        if not os.path.exists(filepath):
            print(f"WARNING: {filepath} not found, skipping")
            continue
        
        print(f"Processing {species}...")
        results = extract_scores(filepath, species)
        print(f"  Extracted {len(results)} proteins")
        all_results.extend(results)
    
    # Write output
    out_tsv = os.path.join(out_dir, 'llphyscore_summary.tsv')
    with open(out_tsv, 'w') as f:
        f.write('protein_id\tspecies\traw_score\tpercentile\tzscore\tmod_zscore\n')
        for r in all_results:
            f.write(f"{r['protein_id']}\t{r['species']}\t{r['raw_score']}\t{r['percentile']}\t{r['zscore']}\t{r['mod_zscore']}\n")
    
    print(f"\n{'='*60}")
    print(f"DONE: {len(all_results)} proteins extracted")
    print(f"Output: {out_tsv}")
    
    # Quick stats
    for species in species_files:
        n = sum(1 for r in all_results if r['species'] == species)
        print(f"  {species}: {n}")
    print(f"{'='*60}")

if __name__ == '__main__':
    main()
