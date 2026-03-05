#!/usr/bin/env python3
"""
PhaSePred JSONから新EuRBPDBリストのRBPスコアを抽出（修正版）

JSON構造: { "Q9H553": {"catGRANULE": ..., "PScore": ..., ...}, ... }

使い方:
  python3 extract_phasepred_v2.py /Users/kyotayasuda/rbp_pfam/output/eurbpdb_pfam_mapping.tsv /Users/kyotayasuda/Desktop/python/phasepred_work/
"""

import sys
import os
import json
import csv

def load_accessions(pfam_tsv):
    acc_by_species = {}
    with open(pfam_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            sp = row.get('species', '')
            acc = row.get('accession', '')
            if sp and acc and acc != 'NOT_FOUND' and sp != 'species':
                if sp not in acc_by_species:
                    acc_by_species[sp] = set()
                acc_by_species[sp].add(acc)
    return acc_by_species


def extract_scores(json_path, target_accs, species_name):
    print(f"  Loading {os.path.basename(json_path)}...")
    with open(json_path) as f:
        data = json.load(f)
    
    print(f"  JSON entries: {len(data)}")
    
    results = []
    found = 0
    
    for acc in target_accs:
        if acc in data:
            entry = data[acc]
            found += 1
            
            result = {
                'accession': acc,
                'species': species_name,
                'gene_name': entry.get('Gene names', ''),
            }
            
            # Extract predictor scores
            for key in ['catGRANULE', 'PScore', 'PLAAC', 'DeepPhase']:
                val = entry.get(key, {})
                if isinstance(val, dict):
                    # These might have sub-keys like 'score', 'rank', etc.
                    for subkey, subval in val.items():
                        result[f'{key}_{subkey}'] = subval
                else:
                    result[key] = val
            
            results.append(result)
    
    print(f"  Matched: {found}/{len(target_accs)} RBPs ({found/max(len(target_accs),1)*100:.1f}%)")
    return results


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 extract_phasepred_v2.py <pfam_mapping.tsv> <phasepred_dir>")
        sys.exit(1)
    
    pfam_tsv = sys.argv[1]
    phasepred_dir = sys.argv[2]
    
    print("=== Loading EuRBPDB accessions ===")
    acc_by_species = load_accessions(pfam_tsv)
    for sp, accs in sorted(acc_by_species.items()):
        print(f"  {sp}: {len(accs)}")
    
    species_json = {
        'Homo_sapiens': 'human_reviewed.json',
        'Mus_musculus': 'mouse_reviewed.json',
        'Danio_rerio': 'zebrafish_reviewed.json',
        'Drosophila_melanogaster': 'fruit-fly_reviewed.json',
        'Caenorhabditis_elegans': 'caenorhabditis-elegans_reviewed.json',
    }
    
    print("\n=== Extracting PhaSePred scores ===")
    all_results = []
    
    for species, json_file in species_json.items():
        json_path = os.path.join(phasepred_dir, json_file)
        if not os.path.exists(json_path):
            print(f"  WARNING: {json_path} not found!")
            continue
        
        target_accs = acc_by_species.get(species, set())
        if not target_accs:
            continue
        
        results = extract_scores(json_path, target_accs, species)
        all_results.extend(results)
    
    print(f"\n=== Total: {len(all_results)} RBPs with PhaSePred scores ===")
    
    if all_results:
        all_keys = sorted(set().union(*[r.keys() for r in all_results]))
        print(f"Columns: {all_keys}")
        print(f"\nSample entry:")
        for k, v in sorted(all_results[0].items()):
            val_str = str(v)[:80]
            print(f"  {k}: {val_str}")
        
        out_path = os.path.join(os.path.dirname(pfam_tsv), 'phasepred_eurbpdb_scores.tsv')
        with open(out_path, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=all_keys, delimiter='\t')
            writer.writeheader()
            for r in all_results:
                writer.writerow({k: r.get(k, '') for k in all_keys})
        print(f"\nSaved: {out_path}")


if __name__ == '__main__':
    main()
