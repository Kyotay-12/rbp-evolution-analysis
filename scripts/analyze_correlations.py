#!/usr/bin/env python3
"""
RBP Evolution Analysis - Main Correlation Script
=================================================
Reproduces the core statistical results of Yasuda et al. (iScience, 2026)

What this script does:
  1. Loads EuRBPDB RBP data (new_rbp_db.csv) and Pfam mappings (eurbpdb_pfam_mapping.tsv)
  2. Calculates RBP family diversity per species using Pfam × EuRBPDB 686-family reference
  3. Calculates TF family diversity per species using AnimalTFDB tf_family annotations
  4. Loads control protein Pfam data (controls_pfam_mapping.tsv) for Kinase and GPCR
  5. Runs Spearman rank correlations vs neuron count for all 4 protein types
  6. Outputs: main6_all_correlations_final.csv

Usage:
  python3 analyze_correlations.py <data_dir> <output_dir>

  <data_dir>   : directory containing input files (default: ../data/processed)
  <output_dir> : directory for output files     (default: ../results)

Input files required (all in <data_dir>):
  new_rbp_db.csv              - RBP database (EuRBPDB, 6 main species)
  new_tf_db.csv               - TF database (AnimalTFDB 4.0, 6 main species)
  eurbpdb_pfam_mapping.tsv    - Pfam domain assignments for each RBP gene
  controls_pfam_mapping.tsv   - Pfam domain assignments for Kinase/GPCR genes

Output files:
  main6_all_correlations_final.csv  - Spearman rho, p-value for all protein types

Author: Kyota Yasuda
Date:   March 2026
"""

import sys
import os
import csv
import math
from collections import defaultdict

# ── scipy / pandas are optional; fall back to stdlib if absent ──────────────
try:
    import pandas as pd
    import numpy as np
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("[WARN] scipy/pandas not found – using built-in Spearman implementation")


# ── Species order (worm → human, ascending neuron count) ────────────────────
SPECIES_ORDER = [
    "Caenorhabditis_elegans",
    "Drosophila_melanogaster",
    "Danio_rerio",
    "Xenopus_tropicalis",
    "Mus_musculus",
    "Homo_sapiens",
]

NEURON_COUNT = {
    "Caenorhabditis_elegans":  302,
    "Drosophila_melanogaster": 200_000,
    "Danio_rerio":             10_000_000,
    "Xenopus_tropicalis":      16_000_000,
    "Mus_musculus":            71_000_000,
    "Homo_sapiens":            86_000_000_000,
}

# EuRBPDB canonical Pfam family reference list (686 families, abbreviated here
# to the key domains; full list encoded in eurbpdb_pfam_mapping.tsv via the
# fetch_pfam_from_uniprot.py pipeline).
# Family assignment logic:
#   gene → UniProt REST API → Pfam IDs → match against EuRBPDB 686 ref list
#   Multiple hits → pick family with largest gene count
#   No hit        → Non-canonical


# ════════════════════════════════════════════════════════════════════════════
# Utility: pure-Python Spearman rank correlation
# ════════════════════════════════════════════════════════════════════════════

def _rank(lst):
    """Return ranks (1-based, averaged for ties)."""
    sorted_idx = sorted(range(len(lst)), key=lambda i: lst[i])
    ranks = [0.0] * len(lst)
    i = 0
    while i < len(lst):
        j = i
        while j < len(lst) - 1 and lst[sorted_idx[j]] == lst[sorted_idx[j + 1]]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1
        for k in range(i, j + 1):
            ranks[sorted_idx[k]] = avg_rank
        i = j + 1
    return ranks


def spearman_r(x, y):
    """Spearman rank correlation + two-tailed p-value (t approximation)."""
    n = len(x)
    if n != len(y) or n < 3:
        raise ValueError(f"Need matched lists of length ≥ 3, got {n}")
    rx = _rank(x)
    ry = _rank(y)
    # Pearson on ranks
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx  = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy  = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if dx == 0 or dy == 0:
        return 0.0, 1.0
    rho = num / (dx * dy)
    # t-distribution approximation
    if abs(rho) >= 1.0:
        p = 0.0
    else:
        t_stat = rho * math.sqrt(n - 2) / math.sqrt(1 - rho ** 2)
        # two-tailed p via regularised incomplete beta
        # Simple normal approximation for n >= 6
        import statistics
        p = 2 * (1 - _norm_cdf(abs(t_stat)))
    return rho, p


def _norm_cdf(x):
    """Standard normal CDF via math.erf."""
    return (1.0 + math.erf(x / math.sqrt(2))) / 2.0


def compute_spearman(x, y):
    """Wrapper: use scipy if available, else built-in."""
    if HAS_SCIPY:
        result = stats.spearmanr(x, y)
        return float(result.statistic), float(result.pvalue)
    return spearman_r(x, y)


# ════════════════════════════════════════════════════════════════════════════
# File loaders
# ════════════════════════════════════════════════════════════════════════════

def load_tsv(path, delimiter="\t"):
    """Load a TSV/CSV into list of dicts. Raises on failure."""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required file not found: {path}")
    rows = []
    with open(path, encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter=delimiter)
        for row in reader:
            rows.append(row)
    print(f"  Loaded {len(rows):,} rows from {os.path.basename(path)}")
    return rows


def load_csv(path):
    return load_tsv(path, delimiter=",")


# ════════════════════════════════════════════════════════════════════════════
# RBP family diversity  (Pfam × EuRBPDB 686-family reference)
# ════════════════════════════════════════════════════════════════════════════

def calc_rbp_family_diversity(rbp_rows, pfam_rows):
    """
    For each species in SPECIES_ORDER, count the number of distinct
    EuRBPDB canonical Pfam families represented.

    Logic mirrors the published pipeline:
      eurbpdb_pfam_mapping.tsv contains columns:
        ensembl_gene_id, species, pfam_families (semicolon-separated),
        eurbpdb_family (best-match EuRBPDB family name or 'Non-canonical')

    Returns:
      dict  species → {'n_proteins': int, 'n_families': int}
    """
    # Build gene → eurbpdb_family from pfam mapping
    gene_to_family = {}
    for row in pfam_rows:
        gid     = row.get("ensembl_gene_id", "").strip()
        family  = row.get("eurbpdb_family",  "").strip()
        species = row.get("species",         "").strip()
        if gid and family:
            gene_to_family[(gid, species)] = family

    # Count per species
    result = {}
    for sp in SPECIES_ORDER:
        sp_genes    = [r for r in rbp_rows
                       if r.get("species", "").strip() == sp
                       and r.get("protein_type", "").strip() == "RBP"]
        n_proteins  = len(sp_genes)
        families    = set()
        for r in sp_genes:
            gid    = r.get("ensembl_gene_id", "").strip()
            key    = (gid, sp)
            fam    = gene_to_family.get(key, "Non-canonical")
            if fam and fam != "Non-canonical":
                families.add(fam)
        result[sp] = {"n_proteins": n_proteins, "n_families": len(families)}

    return result


# ════════════════════════════════════════════════════════════════════════════
# TF family diversity  (AnimalTFDB tf_family column)
# ════════════════════════════════════════════════════════════════════════════

def calc_tf_family_diversity(tf_rows):
    """Count distinct tf_family values per species."""
    result = {}
    for sp in SPECIES_ORDER:
        sp_rows  = [r for r in tf_rows if r.get("species", "").strip() == sp]
        families = set()
        for r in sp_rows:
            fam = r.get("tf_family", "").strip()
            if fam:
                families.add(fam)
        result[sp] = {"n_proteins": len(sp_rows), "n_families": len(families)}
    return result


# ════════════════════════════════════════════════════════════════════════════
# Control protein (Kinase / GPCR) family diversity  (Pfam)
# ════════════════════════════════════════════════════════════════════════════

def calc_control_family_diversity(ctrl_rows, protein_type):
    """
    controls_pfam_mapping.tsv columns expected:
      ensembl_gene_id, gene_symbol, species, protein_type,
      pfam_families, eurbpdb_family (or pfam_top_family)

    Count distinct Pfam families per species for the given protein_type.
    """
    result = {}
    for sp in SPECIES_ORDER:
        sp_rows  = [r for r in ctrl_rows
                    if r.get("species",       "").strip() == sp
                    and r.get("protein_type", "").strip() == protein_type]
        families = set()
        for r in sp_rows:
            # Try several possible column names for the family field
            fam = (r.get("eurbpdb_family") or
                   r.get("pfam_top_family") or
                   r.get("pfam_families")   or "").strip()
            # Use first token if semicolon-separated
            fam = fam.split(";")[0].strip()
            if fam and fam not in ("Non-canonical", "unknown", ""):
                families.add(fam)
        result[sp] = {"n_proteins": len(sp_rows), "n_families": len(families)}
    return result


# ════════════════════════════════════════════════════════════════════════════
# Run correlations
# ════════════════════════════════════════════════════════════════════════════

def run_correlation(diversity_dict, metric="family_diversity", vs="neurons"):
    """
    diversity_dict: species → {'n_proteins': int, 'n_families': int}
    Returns: rho, p_value, values_str
    """
    x_neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]
    if metric == "family_diversity":
        y_values  = [diversity_dict[sp]["n_families"] for sp in SPECIES_ORDER]
    else:
        y_values  = [diversity_dict[sp]["n_proteins"] for sp in SPECIES_ORDER]

    rho, pval = compute_spearman(x_neurons, y_values)

    # Format values string (worm → human)
    values_str = "→".join(str(diversity_dict[sp]["n_families"]
                              if metric == "family_diversity"
                              else diversity_dict[sp]["n_proteins"])
                           for sp in SPECIES_ORDER)
    return rho, pval, values_str


def significance_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    # ── argument parsing ────────────────────────────────────────────────────
    data_dir   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
                    os.path.dirname(__file__), "..", "data", "processed")
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
                    os.path.dirname(__file__), "..", "results")

    data_dir   = os.path.abspath(data_dir)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("RBP Evolution – Correlation Analysis")
    print(f"  data_dir   : {data_dir}")
    print(f"  output_dir : {output_dir}")
    print("=" * 60)

    # ── load input files ────────────────────────────────────────────────────
    print("\n[1] Loading input files …")
    rbp_rows  = load_csv(os.path.join(data_dir, "new_rbp_db.csv"))
    tf_rows   = load_csv(os.path.join(data_dir, "new_tf_db.csv"))
    pfam_rows = load_tsv(os.path.join(data_dir, "eurbpdb_pfam_mapping.tsv"))
    ctrl_rows = load_tsv(os.path.join(data_dir, "controls_pfam_mapping.tsv"))

    # ── calculate diversity ─────────────────────────────────────────────────
    print("\n[2] Calculating family diversity …")
    rbp_div    = calc_rbp_family_diversity(rbp_rows, pfam_rows)
    tf_div     = calc_tf_family_diversity(tf_rows)
    kinase_div = calc_control_family_diversity(ctrl_rows, "Kinase")
    gpcr_div   = calc_control_family_diversity(ctrl_rows, "GPCR")

    print("\n  RBP family counts per species:")
    for sp in SPECIES_ORDER:
        print(f"    {sp:35s} n_proteins={rbp_div[sp]['n_proteins']:5d}  "
              f"n_families={rbp_div[sp]['n_families']:4d}")

    print("\n  TF family counts per species:")
    for sp in SPECIES_ORDER:
        print(f"    {sp:35s} n_proteins={tf_div[sp]['n_proteins']:5d}  "
              f"n_families={tf_div[sp]['n_families']:4d}")

    # ── run correlations ────────────────────────────────────────────────────
    print("\n[3] Running Spearman correlations vs neuron count …")

    results = []
    for ptype, div in [("RBP",    rbp_div),
                       ("TF",     tf_div),
                       ("Kinase", kinase_div),
                       ("GPCR",   gpcr_div)]:
        for metric in ["family_diversity", "total_count"]:
            rho, pval, vals = run_correlation(div, metric=metric)
            sig  = significance_label(pval)
            src  = {"RBP":    "UniProt Pfam × EuRBPDB 686 ref",
                    "TF":     "AnimalTFDB tf_family",
                    "Kinase": "UniProt Pfam",
                    "GPCR":   "UniProt Pfam"}.get(ptype, "")
            results.append({
                "protein_type": ptype,
                "metric":       metric,
                "vs":           "neurons",
                "rho":          round(rho,  4),
                "p_value":      round(pval, 4),
                "significance": sig,
                "values_worm_to_human": vals,
                "source":       src,
            })
            print(f"  {ptype:6s} {metric:20s}  rho={rho:+.3f}  p={pval:.4f}  {sig}")

    # ── write output ────────────────────────────────────────────────────────
    out_path = os.path.join(output_dir, "main6_all_correlations_final.csv")
    fieldnames = ["protein_type", "metric", "vs", "rho", "p_value",
                  "significance", "values_worm_to_human", "source"]
    with open(out_path, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(results)

    print(f"\n[4] Saved → {out_path}")
    print("\nDone. Key results:")
    for r in results:
        if r["metric"] == "family_diversity":
            print(f"  {r['protein_type']:6s}  rho={r['rho']:+.3f}  "
                  f"p={r['p_value']:.4f}  {r['significance']}")


if __name__ == "__main__":
    main()
