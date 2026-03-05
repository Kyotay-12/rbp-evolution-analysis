#!/usr/bin/env python3
"""
MobiDB APIからIDR情報を取得（修正版）

入力: eurbpdb_pfam_mapping.tsv
出力: idr_mobidb_results_v2.tsv

使い方:
  python3 fetch_idr_mobidb_v2.py <eurbpdb_pfam_mapping.tsv> <output_dir>

レジューム対応。所要時間: 1.5-2時間
"""

import sys
import os
import csv
import json
import time
import urllib.request
import ssl

ssl._create_default_https_context = ssl._create_unverified_context


def fetch_mobidb(accession, retries=3):
    url = f'https://mobidb.bio.unipd.it/api/download?acc={accession}&format=json'
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url)
            req.add_header('Accept', 'application/json')
            with urllib.request.urlopen(req, timeout=30) as resp:
                return json.loads(resp.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            if e.code == 404:
                return None
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
        except Exception:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
    return None


def extract_idr(data, accession):
    """MobiDB JSONからIDR情報を抽出
    
    MobiDB response structure:
      {
        "acc": "P09651",
        "length": 372,
        "prediction-disorder-mobidb_lite": {  <-- consensus IDR prediction
          "regions": [[start, end], ...],
          "content_fraction": 0.xxx,
          "content_count": N
        },
        "curated-disorder-merge": { ... },    <-- curated IDR
        ...
      }
    """
    if not data:
        return None

    result = {
        'accession': accession,
        'length': data.get('length', 0),
        'disorder_content': 0.0,
        'n_disordered_residues': 0,
        'n_disordered_regions': 0,
        'longest_idr': 0,
        'idr_source': '',
    }

    # Priority order for disorder data:
    # 1. prediction-disorder-mobidb_lite (consensus, most widely used)
    # 2. curated-disorder-merge (curated, fewer entries)
    # 3. prediction-disorder-th_50 (threshold-based)
    # 4. homology-disorder-merge
    source_keys = [
        'prediction-disorder-mobidb_lite',
        'curated-disorder-merge',
        'prediction-disorder-th_50',
        'homology-disorder-merge',
    ]

    for key in source_keys:
        if key in data and data[key]:
            entry = data[key]
            regions = entry.get('regions', [])
            content_fraction = entry.get('content_fraction', 0)
            content_count = entry.get('content_count', 0)

            if regions:
                longest = 0
                for r in regions:
                    if isinstance(r, list) and len(r) >= 2:
                        rlen = r[1] - r[0] + 1
                        longest = max(longest, rlen)

                result['disorder_content'] = round(content_fraction, 4) if content_fraction else 0.0
                result['n_disordered_residues'] = content_count if content_count else 0
                result['n_disordered_regions'] = len(regions)
                result['longest_idr'] = longest
                result['idr_source'] = key

                # If content_fraction not provided, calculate from regions
                if not content_fraction and result['length'] > 0 and content_count:
                    result['disorder_content'] = round(content_count / result['length'], 4)
                elif not content_fraction and result['length'] > 0:
                    total_res = sum(r[1] - r[0] + 1 for r in regions if isinstance(r, list) and len(r) >= 2)
                    result['disorder_content'] = round(total_res / result['length'], 4)
                    result['n_disordered_residues'] = total_res

                break  # Use first available source

    return result


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 fetch_idr_mobidb_v2.py <eurbpdb_pfam_mapping.tsv> <output_dir>")
        sys.exit(1)

    import pandas as pd

    input_tsv = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    pfam = pd.read_csv(input_tsv, sep='\t')
    pfam = pfam[(pfam['accession'].notna()) & (pfam['accession'] != 'NOT_FOUND') & (pfam['accession'] != '')]

    unique_acc = pfam[['accession', 'species']].drop_duplicates(subset='accession')
    print(f"Total unique accessions: {len(unique_acc)}")

    fieldnames = ['accession', 'species', 'length', 'disorder_content',
                  'n_disordered_residues', 'n_disordered_regions',
                  'longest_idr', 'idr_source']

    out_tsv = os.path.join(out_dir, 'idr_mobidb_results_v2.tsv')

    # Resume
    done_accs = set()
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                done_accs.add(row['accession'])
        print(f"Resuming: {len(done_accs)} already done")

    write_header = len(done_accs) == 0
    tsv_file = open(out_tsv, 'a', newline='')
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    if write_header:
        writer.writeheader()

    acc_species = dict(zip(unique_acc['accession'], unique_acc['species']))
    accs_to_process = [a for a in unique_acc['accession'] if a not in done_accs]
    print(f"To process: {len(accs_to_process)}")

    found = 0
    with_idr = 0

    for i, acc in enumerate(accs_to_process):
        sp = acc_species.get(acc, '')
        data = fetch_mobidb(acc)
        result = extract_idr(data, acc)

        if result:
            result['species'] = sp
            writer.writerow({k: result.get(k, '') for k in fieldnames})
            found += 1
            if result['disorder_content'] > 0:
                with_idr += 1
        else:
            writer.writerow({
                'accession': acc, 'species': sp, 'length': '',
                'disorder_content': '', 'n_disordered_residues': '',
                'n_disordered_regions': '', 'longest_idr': '', 'idr_source': '',
            })

        if (i + 1) % 200 == 0:
            print(f"  {i+1}/{len(accs_to_process)} | found={found} with_idr={with_idr}")
            tsv_file.flush()

        time.sleep(0.4)

    tsv_file.close()
    total = len(accs_to_process)
    print(f"\n{'='*60}")
    print(f"DONE: {found}/{total} found ({found/max(total,1)*100:.1f}%)")
    print(f"With IDR>0: {with_idr}/{found} ({with_idr/max(found,1)*100:.1f}%)")
    print(f"Output: {out_tsv}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
