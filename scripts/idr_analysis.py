#!/usr/bin/env python3
"""
IDR Analysis Script
====================
Reproduces IDR (Intrinsically Disordered Region) results in Yasuda et al. (iScience, 2026)
Corresponds to Main Figure 4.

What this script does:
  1. Loads per-gene IDR data from MobiDB (idr_mobidb_results_v2.tsv)
  2. Computes species-level summary statistics:
       - Mean IDR content (%)           → Figure 4A  [ns, rho=0.714]
       - Mean longest IDR length (aa)   → Figure 4B  [*, rho=0.829, p=0.042]
       - Mean total disordered residues → Figure 4C  [*, rho=0.829, p=0.042]
  3. Runs Spearman correlations vs neuron count
  4. Outputs: idr_species_summary.csv, idr_correlations.csv

Usage:
  python3 idr_analysis.py <data_dir> <output_dir>

  <data_dir>   : directory containing input files (default: ../data/processed)
  <output_dir> : directory for output files     (default: ../results)

Input files required:
  idr_mobidb_results_v2.tsv   - per-gene MobiDB disorder predictions
  new_rbp_db.csv              - RBP database (for species membership)

Expected columns in idr_mobidb_results_v2.tsv:
  uniprot_acc, gene_symbol, species,
  idr_pct (% disordered residues),
  longest_idr_aa (length of longest contiguous IDR),
  n_disordered_residues (total disordered aa count),
  n_idr_regions (number of distinct IDR regions)

Output files:
  idr_species_summary.csv    - species-level IDR statistics
  idr_correlations.csv       - Spearman correlations vs neurons

Author: Kyota Yasuda
Date:   March 2026
"""

import sys
import os
import csv
import math
from collections import defaultdict

try:
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("[WARN] scipy not found – using built-in Spearman implementation")


# ── Species metadata ─────────────────────────────────────────────────────────

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


# ════════════════════════════════════════════════════════════════════════════
# Spearman utilities
# ════════════════════════════════════════════════════════════════════════════

def _rank(lst):
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


def _norm_cdf(x):
    return (1.0 + math.erf(x / math.sqrt(2))) / 2.0


def spearman_builtin(x, y):
    n = len(x)
    rx, ry = _rank(x), _rank(y)
    mx = sum(rx) / n
    my = sum(ry) / n
    num = sum((rx[i] - mx) * (ry[i] - my) for i in range(n))
    dx  = math.sqrt(sum((rx[i] - mx) ** 2 for i in range(n)))
    dy  = math.sqrt(sum((ry[i] - my) ** 2 for i in range(n)))
    if dx == 0 or dy == 0:
        return 0.0, 1.0
    rho = num / (dx * dy)
    if abs(rho) >= 1.0:
        return rho, 0.0
    t_stat = rho * math.sqrt(n - 2) / math.sqrt(1 - rho ** 2)
    p = 2 * (1 - _norm_cdf(abs(t_stat)))
    return rho, p


def spearman(x, y):
    if HAS_SCIPY:
        r = stats.spearmanr(x, y)
        return float(r.statistic), float(r.pvalue)
    return spearman_builtin(x, y)


def sig_label(p):
    if p < 0.001: return "***"
    if p < 0.01:  return "**"
    if p < 0.05:  return "*"
    return "ns"


# ════════════════════════════════════════════════════════════════════════════
# File loaders
# ════════════════════════════════════════════════════════════════════════════

def load_csv(path, delimiter=","):
    if not os.path.isfile(path):
        raise FileNotFoundError(f"Required file not found: {path}")
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh, delimiter=delimiter))
    print(f"  Loaded {len(rows):,} rows ← {os.path.basename(path)}")
    return rows


def load_tsv(path):
    return load_csv(path, delimiter="\t")


def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved → {path}")


# ════════════════════════════════════════════════════════════════════════════
# IDR summary computation
# ════════════════════════════════════════════════════════════════════════════

def compute_species_summary(idr_rows, rbp_genes_per_species):
    """
    Aggregate per-gene IDR data to species-level means.

    idr_rows: list of dicts from idr_mobidb_results_v2.tsv
    rbp_genes_per_species: dict species → set of gene_symbols (from new_rbp_db)
      Used to restrict to EuRBPDB RBPs only.

    Returns list of dicts (one per species in SPECIES_ORDER).
    """
    # Group IDR rows by species
    by_species = defaultdict(list)
    for row in idr_rows:
        sp = row.get("species", "").strip()
        if sp in SPECIES_ORDER:
            by_species[sp].append(row)

    summary = []
    for sp in SPECIES_ORDER:
        rows = by_species[sp]
        rbp_set = rbp_genes_per_species.get(sp, set())

        # Filter to RBPs present in EuRBPDB if rbp_set is available
        if rbp_set:
            rows = [r for r in rows
                    if r.get("gene_symbol", "").strip() in rbp_set
                    or r.get("uniprot_acc",  "").strip() in rbp_set]

        if not rows:
            raise ValueError(
                f"No IDR data found for {sp}. "
                f"Check that idr_mobidb_results_v2.tsv contains this species."
            )

        def safe_float(v):
            try:
                return float(v)
            except (ValueError, TypeError):
                return None

        # Extract numeric columns
        idr_pcts      = [safe_float(r.get("idr_pct"))              for r in rows]
        longest_idrs  = [safe_float(r.get("longest_idr_aa"))       for r in rows]
        disord_res    = [safe_float(r.get("n_disordered_residues")) for r in rows]
        n_regions     = [safe_float(r.get("n_idr_regions"))        for r in rows]

        def mean_nonnull(lst):
            valid = [v for v in lst if v is not None]
            return sum(valid) / len(valid) if valid else None

        def median_nonnull(lst):
            valid = sorted(v for v in lst if v is not None)
            if not valid:
                return None
            n = len(valid)
            return (valid[n // 2] if n % 2 else
                    (valid[n // 2 - 1] + valid[n // 2]) / 2)

        def sd_nonnull(lst):
            valid = [v for v in lst if v is not None]
            if len(valid) < 2:
                return None
            m = sum(valid) / len(valid)
            return math.sqrt(sum((v - m) ** 2 for v in valid) / (len(valid) - 1))

        n_with_data = sum(1 for v in idr_pcts if v is not None)

        summary.append({
            "species":                  sp,
            "n_rbps_total":             len(rows),
            "n_rbps_with_idr_data":     n_with_data,
            "coverage_pct":             round(100.0 * n_with_data / len(rows), 1)
                                        if rows else 0,
            "neurons":                  NEURON_COUNT[sp],
            "mean_idr_pct":             round(mean_nonnull(idr_pcts),     2)
                                        if mean_nonnull(idr_pcts) is not None else "NA",
            "median_idr_pct":           round(median_nonnull(idr_pcts),   2)
                                        if median_nonnull(idr_pcts) is not None else "NA",
            "sd_idr_pct":              round(sd_nonnull(idr_pcts),       2)
                                        if sd_nonnull(idr_pcts) is not None else "NA",
            "mean_longest_idr_aa":      round(mean_nonnull(longest_idrs), 1)
                                        if mean_nonnull(longest_idrs) is not None else "NA",
            "mean_n_idr_regions":       round(mean_nonnull(n_regions),    2)
                                        if mean_nonnull(n_regions) is not None else "NA",
            "mean_disordered_residues": round(mean_nonnull(disord_res),   1)
                                        if mean_nonnull(disord_res) is not None else "NA",
        })

    return summary


# ════════════════════════════════════════════════════════════════════════════
# Correlation analysis
# ════════════════════════════════════════════════════════════════════════════

def compute_correlations(summary):
    """
    Run Spearman correlations for each IDR metric vs neurons.
    Corresponds to Figure 4A (IDR%), 4B (longest IDR), 4C (disordered residues).
    """
    neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]

    metrics = [
        ("Mean IDR content (%)",         "mean_idr_pct",
         "Figure 4A – expected ns (rho≈0.714)"),
        ("Median IDR content (%)",        "median_idr_pct",
         "Figure 4A alternative"),
        ("Mean longest IDR length (aa)",  "mean_longest_idr_aa",
         "Figure 4B – expected * (rho=0.829, p=0.042)"),
        ("Mean disordered residues (aa)", "mean_disordered_residues",
         "Figure 4C – expected * (rho=0.829, p=0.042)"),
        ("Mean n IDR regions",            "mean_n_idr_regions",
         "Supplementary"),
    ]

    results = []
    for label, col, note in metrics:
        y_vals = []
        valid  = True
        for row in summary:
            v = row.get(col, "NA")
            if v == "NA":
                print(f"  [WARN] Missing {col} for {row['species']} – skipping metric")
                valid = False
                break
            y_vals.append(float(v))

        if not valid:
            continue

        rho, pval = spearman(neurons, y_vals)
        sig       = sig_label(pval)
        vals_str  = "→".join(str(v) for v in y_vals)

        results.append({
            "metric":              label,
            "vs":                  "neurons",
            "rho":                 round(rho,  3),
            "p_value":             round(pval, 4),
            "significance":        sig,
            "values_worm_to_human": vals_str,
            "note":                note,
        })
        print(f"  {label:40s}  rho={rho:+.3f}  p={pval:.4f}  {sig}")

    return results


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    data_dir   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
                    os.path.dirname(__file__), "..", "data", "processed")
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
                    os.path.dirname(__file__), "..", "results")

    data_dir   = os.path.abspath(data_dir)
    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("IDR Analysis (Main Figure 4)")
    print(f"  data_dir   : {data_dir}")
    print(f"  output_dir : {output_dir}")
    print("=" * 60)

    # ── Load data ────────────────────────────────────────────────────────────
    print("\n[1] Loading data …")
    idr_rows = load_tsv(os.path.join(data_dir, "idr_mobidb_results_v2.tsv"))
    rbp_rows = load_csv(os.path.join(data_dir, "new_rbp_db.csv"))

    # Build species → gene symbol set from EuRBPDB
    rbp_genes = defaultdict(set)
    for r in rbp_rows:
        sp  = r.get("species", "").strip()
        sym = r.get("gene_symbol", "").strip()
        if sp and sym:
            rbp_genes[sp].add(sym)

    # ── Compute species summary ──────────────────────────────────────────────
    print("\n[2] Computing species-level IDR statistics …")
    summary = compute_species_summary(idr_rows, rbp_genes)

    print(f"\n  {'Species':35s} {'n_RBPs':>7} {'cov%':>6} "
          f"{'IDR%':>7} {'longest(aa)':>12} {'disord(aa)':>11}")
    for r in summary:
        print(f"  {r['species']:35s} {r['n_rbps_total']:>7} "
              f"{r['coverage_pct']:>6} {str(r['mean_idr_pct']):>7} "
              f"{str(r['mean_longest_idr_aa']):>12} "
              f"{str(r['mean_disordered_residues']):>11}")

    write_csv(
        os.path.join(output_dir, "idr_species_summary.csv"),
        summary,
        ["species", "n_rbps_total", "n_rbps_with_idr_data", "coverage_pct",
         "neurons", "mean_idr_pct", "median_idr_pct", "sd_idr_pct",
         "mean_longest_idr_aa", "mean_n_idr_regions", "mean_disordered_residues"],
    )

    # ── Run correlations ─────────────────────────────────────────────────────
    print("\n[3] Running Spearman correlations vs neuron count …")
    corr_results = compute_correlations(summary)

    write_csv(
        os.path.join(output_dir, "idr_correlations.csv"),
        corr_results,
        ["metric", "vs", "rho", "p_value", "significance",
         "values_worm_to_human", "note"],
    )

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Expected results (from manuscript):")
    print("  Mean IDR content (%)        rho=0.714  ns   (Figure 4A)")
    print("  Mean longest IDR (aa)       rho=0.829  *    (Figure 4B)")
    print("  Mean disordered residues    rho=0.829  *    (Figure 4C)")
    print("  Range: worm 45.9 aa → human 78.5 aa")
    print("\nDone.")


if __name__ == "__main__":
    main()
