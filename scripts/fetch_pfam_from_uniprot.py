#!/usr/bin/env python3
"""
EuRBPDB FASTA + _RBP.txt からUniProt Pfam domainを統一取得

配列はFASTAから取得済み。このスクリプトはgene symbolで
UniProt検索してPfam domainのみ取得する（配列DL不要で高速）。

入力ディレクトリに以下を配置:
  Homo_sapiens_RBP.txt       (EuRBPDB, 5列)
  Mus_musculus_RBP.txt       (EuRBPDB, 6列)
  Danio_rerio_RBP.txt        (同上)
  Drosophila_melanogaster_RBP.txt
  Caenorhabditis_elegans_RBP.txt
  Xenopus_tropicalis_RBP.txt

使い方:
  python3 fetch_pfam_from_uniprot.py <dir_with_txt_files> <output_dir>

出力:
  eurbpdb_pfam_mapping.tsv  (ensembl_id → UniProt accession → Pfam domains)

所要時間: 2-3時間（約13,300件）
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

SPECIES = {
    'Homo_sapiens':            {'taxon': '9606',  'ncols': 5},
    'Mus_musculus':            {'taxon': '10090', 'ncols': 6},
    'Danio_rerio':             {'taxon': '7955',  'ncols': 6},
    'Drosophila_melanogaster': {'taxon': '7227',  'ncols': 6},
    'Caenorhabditis_elegans':  {'taxon': '6239',  'ncols': 6},
    'Xenopus_tropicalis':      {'taxon': '8364',  'ncols': 6},
}


def _uniprot_search(query, retries=3):
    params = {
        'query': query,
        'format': 'json',
        'fields': 'accession,gene_names,length,xref_pfam',
        'size': 5,
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
                    reviewed = [r for r in results
                                if r.get('entryType') == 'UniProtKB reviewed (Swiss-Prot)']
                    return reviewed[0] if reviewed else results[0]
                return None
        except Exception as e:
            if attempt < retries - 1:
                time.sleep(2 * (attempt + 1))
            else:
                return None


def fetch_by_gene(gene_symbol, taxon_id):
    return _uniprot_search(f'(gene_exact:{gene_symbol}) AND (organism_id:{taxon_id})')


def fetch_by_ensembl(ensembl_id, taxon_id):
    return _uniprot_search(f'(xref:ensembl-{ensembl_id}) AND (organism_id:{taxon_id})')


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
    reviewed = 'reviewed' if 'reviewed' in entry.get('entryType', '').lower() else 'unreviewed'
    return {
        'accession': acc, 'gene_name_uniprot': gene_name,
        'length': length, 'pfam_ids': ';'.join(pfam_list),
        'reviewed': reviewed,
    }


def load_eurbpdb_txt(filepath, ncols):
    genes = []
    with open(filepath, encoding='utf-8') as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            ens_id = parts[0] if len(parts) > 0 else ''
            gene_sym = parts[1] if len(parts) > 1 else ''
            if gene_sym == '\\N':
                gene_sym = ''
            rbd = parts[5].strip() if ncols == 6 and len(parts) >= 6 else ''
            genes.append((ens_id, gene_sym, rbd))
    return genes


def main():
    if len(sys.argv) < 3:
        print("Usage: python3 fetch_pfam_from_uniprot.py <dir_with_txt_files> <output_dir>")
        sys.exit(1)

    input_dir = sys.argv[1]
    out_dir = sys.argv[2]
    os.makedirs(out_dir, exist_ok=True)

    fieldnames = ['ensembl_id', 'query_gene', 'species', 'eurbpdb_rbd',
                  'accession', 'gene_name_uniprot', 'length',
                  'pfam_ids', 'reviewed']

    out_tsv = os.path.join(out_dir, 'eurbpdb_pfam_mapping.tsv')

    # Resume support
    done_keys = set()
    if os.path.exists(out_tsv):
        with open(out_tsv) as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                done_keys.add((row['ensembl_id'], row['species']))
        print(f"Resuming: {len(done_keys)} already done")

    write_header = len(done_keys) == 0
    tsv_file = open(out_tsv, 'a', newline='')
    writer = csv.DictWriter(tsv_file, fieldnames=fieldnames, delimiter='\t')
    if write_header:
        writer.writeheader()

    total_queries = 0
    total_found = 0
    total_with_pfam = 0

    for sp, info in SPECIES.items():
        filename = f'{sp}.RBP.txt'
        filepath = os.path.join(input_dir, filename)
        if not os.path.exists(filepath):
            print(f"SKIP: {filepath} not found")
            continue

        genes = load_eurbpdb_txt(filepath, info['ncols'])
        taxon = info['taxon']
        sp_found = 0
        sp_pfam = 0
        sp_skip = 0

        print(f"\n{'='*60}")
        print(f"{sp}: {len(genes)} genes ({filename})")
        print(f"{'='*60}")

        for i, (ens_id, gene_sym, rbd) in enumerate(genes):
            if (ens_id, sp) in done_keys:
                sp_skip += 1
                continue

            total_queries += 1

            # Try gene symbol first, then Ensembl ID fallback
            entry = None
            if gene_sym and gene_sym not in ['-', '']:
                entry = fetch_by_gene(gene_sym, taxon)
            if not entry and ens_id:
                entry = fetch_by_ensembl(ens_id, taxon)

            result = extract_pfam(entry)

            if result:
                writer.writerow({
                    'ensembl_id': ens_id, 'query_gene': gene_sym, 'species': sp,
                    'eurbpdb_rbd': rbd, **result,
                })
                total_found += 1
                sp_found += 1
                if result['pfam_ids']:
                    total_with_pfam += 1
                    sp_pfam += 1
            else:
                writer.writerow({
                    'ensembl_id': ens_id, 'query_gene': gene_sym, 'species': sp,
                    'eurbpdb_rbd': rbd, 'accession': 'NOT_FOUND',
                    'gene_name_uniprot': '', 'length': '',
                    'pfam_ids': '', 'reviewed': '',
                })

            if (i + 1 - sp_skip) % 200 == 0:
                print(f"  {i+1}/{len(genes)} processed, {sp_found} found ({sp_pfam} with Pfam)...")

            time.sleep(0.35)

        print(f"  {sp}: {sp_found}/{len(genes)-sp_skip} found, {sp_pfam} with Pfam, {sp_skip} resumed")

    tsv_file.close()

    pct = total_found / max(total_queries, 1) * 100
    pfam_pct = total_with_pfam / max(total_found, 1) * 100
    print(f"\n{'='*60}")
    print(f"DONE: {total_found}/{total_queries} found ({pct:.1f}%)")
    print(f"With Pfam: {total_with_pfam}/{total_found} ({pfam_pct:.1f}%)")
    print(f"Output: {out_tsv}")
    print(f"{'='*60}")


if __name__ == '__main__':
    main()
