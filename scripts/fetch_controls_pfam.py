#!/usr/bin/env python3
"""
Kinase/GPCR の UniProt accession から Pfam domain を取得

入力: integrated_controls_all15species_final.csv (プロジェクトファイル)
出力: controls_pfam_mapping.tsv

使い方:
  python3 fetch_controls_pfam.py <integrated_controls.csv> <output_dir>

所要時間: 30-40分（約3,000件）
レジューム対応
"""

import sys
import os
import csv
import json
import time
import urllib.request
import urllib.parse
import ssl

ssl._create_default_https_context = ssl._create_unverified_context

MAIN6 = ['worm', 'fly', 'zebrafish', 'frog', 'mouse', 'human']


def fetch_pfam_by_accession(accession, retries=3):
    params = {
        'query': f'(accession:{accession})',
        'format': 'json',
        'fields': 'accession,gene_names,organism_name,length,xref_pfam',
        'size': 1,
    }
    url = 'https://rest.uniprot.org/uniprotkb/search?' + urllib.parse.urlencode(params)
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            with urllib.request.urlopen(req, timeout=30) as resp:
                data = json.loads(resp.read().decode('utf-8'))
                results = data.get('results', [])
                if results:
                    return results[0]
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return None


def extract_pfam(entry):
    if not entry:
        return None
    acc = entry.get('primaryAccession', '')
    genes = entry.get('genes', [])
    gene_name = genes[0].get('geneName', {}).get('value', '') if genes else ''
    length = entry.get('sequence', {}).get('length', 0)
    xrefs = entry.get('uniProtKBCrossReferences', [])
    pfam_list = []
    for xref in xrefs:
        if xref.get('database') == 'Pfam':
            xid = xref.get('id', '')
            ename = ''
            for p in xref.get('properties', []):
                if p.get('key') == 'EntryName':
                    ename = p.get('value', '')
            pfam_list.append(f"{xid}:{ename}" if ename else xid)
    return {
        'accession': acc, 'gene_name': gene_name,
        'length': length, 'pfam_ids': ';'.join(pfam_list),
    }


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 fetch_controls_pfam.py <integrated_controls.csv> <output_dir>")
        sys.exit(1)

    import pandas as pd
    controls_csv = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    ctrl = pd.read_csv(controls_csv)
    ctrl = ctrl[(ctrl['species'].isin(MAIN6)) & (ctrl['protein_type'].isin(['Kinase', 'GPCR']))]
    print(f"Total: {len(ctrl)} entries (Kinase + GPCR, main 6 species)")

    fieldnames = ['accession', 'gene_name', 'species', 'protein_type',
                  'old_family', 'length', 'pfam_ids']

    out_tsv = os.path.join(out_dir, 'controls_pfam_mapping.tsv')

    # Resume
    done_keys = set()
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                done_keys.add((row['accession'], row['species']))
        print(f"Resuming: {len(done_keys)} already done")

    write_header = len(done_keys) == 0
    tsv_file = open(out_tsv, 'a', newline='')
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    if write_header:
        writer.writeheader()

    total_queries = 0
    total_found = 0

    for ptype in ['Kinase', 'GPCR']:
        sub = ctrl[ctrl['protein_type'] == ptype]
        sp_found = 0
        sp_skip = 0

        print(f"\n{'='*60}")
        print(f"{ptype}: {len(sub)} entries")
        print(f"{'='*60}")

        for i, (_, row) in enumerate(sub.iterrows()):
            acc = row['accession']
            sp = row['species']
            old_fam = row.get('family', '')

            if (acc, sp) in done_keys:
                sp_skip += 1
                continue

            total_queries += 1
            entry = fetch_pfam_by_accession(acc)
            result = extract_pfam(entry)

            if result:
                writer.writerow({
                    'accession': acc, 'gene_name': result['gene_name'],
                    'species': sp, 'protein_type': ptype,
                    'old_family': old_fam, 'length': result['length'],
                    'pfam_ids': result['pfam_ids'],
                })
                total_found += 1
                sp_found += 1
            else:
                writer.writerow({
                    'accession': acc, 'gene_name': '', 'species': sp,
                    'protein_type': ptype, 'old_family': old_fam,
                    'length': '', 'pfam_ids': '',
                })

            if (i + 1 - sp_skip) % 200 == 0:
                print(f"  {i+1}/{len(sub)} processed, {sp_found} found...")

            time.sleep(0.35)

        print(f"  {ptype}: {sp_found}/{len(sub)-sp_skip} found, {sp_skip} resumed")

    tsv_file.close()
    pct = total_found / max(total_queries, 1) * 100
    print(f"\n{'='*60}")
    print(f"DONE: {total_found}/{total_queries} found ({pct:.1f}%)")
    print(f"Output: {out_tsv}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
