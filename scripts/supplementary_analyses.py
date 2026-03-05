#!/usr/bin/env python3
"""
Supplementary Analyses Script
===============================
Reproduces supplementary analyses in Yasuda et al. (iScience, 2026)

Analyses performed:
  1. Non-canonical RBP analysis
       - Count non-canonical RBPs per species (RBPWorld)
       - Spearman correlation vs neurons
       - Combined canonical + non-canonical correlation
  2. RBP–TF overlap
       - Intersect EuRBPDB gene symbols with AnimalTFDB gene symbols
       - Report overlap % per species
  3. Vertebrate domain expansion
       - Compute vertebrate/invertebrate gene count ratio per Pfam domain
       - Identify domains with ≥3-fold enrichment or vertebrate-specific emergence

Usage:
  python3 supplementary_analyses.py <data_dir> <output_dir>

  <data_dir>   : directory containing input files (default: ../data/processed)
  <output_dir> : directory for output files     (default: ../results)

Input files required:
  new_rbp_db.csv               - RBP database (EuRBPDB, canonical + non-canonical flag)
  new_tf_db.csv                - TF database (AnimalTFDB 4.0)
  eurbpdb_pfam_mapping.tsv     - gene-level Pfam assignments

Output files:
  noncanonical_correlations.csv  - non-canonical RBP counts and correlations
  rbp_tf_overlap.csv             - RBP-TF overlap per species
  domain_expansion.csv           - vertebrate expansion ratios per Pfam domain

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

# RBPWorld available for 5 species (Xenopus excluded)
SPECIES_NONCANONICAL = [
    "Caenorhabditis_elegans",
    "Drosophila_melanogaster",
    "Danio_rerio",
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

# Vertebrate vs invertebrate classification
INVERTEBRATES = {"Caenorhabditis_elegans", "Drosophila_melanogaster"}
VERTEBRATES   = {"Danio_rerio", "Xenopus_tropicalis", "Mus_musculus", "Homo_sapiens"}


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
# 1. Non-canonical RBP analysis
# ════════════════════════════════════════════════════════════════════════════

def run_noncanonical(rbp_rows):
    """
    Uses is_canonical column in new_rbp_db.csv to separate canonical / non-canonical.
    Non-canonical = experimentally identified RBPs lacking canonical RBDs (from RBPWorld).

    Expected result (from manuscript):
      Non-canonical count vs neurons: rho=0.900, p=0.037 (5 species, excl. Xenopus)
      Combined (canonical + non-canonical) vs neurons: rho=0.975, p=0.005
    """
    # Count by species and canonical status
    canonical_counts    = defaultdict(int)
    noncanonical_counts = defaultdict(int)

    for r in rbp_rows:
        sp  = r.get("species", "").strip()
        can = r.get("is_canonical", "").strip().lower()
        if sp not in SPECIES_ORDER:
            continue
        if can in ("true", "1", "yes"):
            canonical_counts[sp] += 1
        elif can in ("false", "0", "no"):
            noncanonical_counts[sp] += 1

    # Check non-canonical data exists
    total_noncan = sum(noncanonical_counts.values())
    if total_noncan == 0:
        print("  [WARN] No non-canonical RBPs found in new_rbp_db.csv.")
        print("         Ensure is_canonical=False entries from RBPWorld are present.")

    # Correlation: non-canonical vs neurons (5 species, excl. Xenopus)
    sp5       = SPECIES_NONCANONICAL
    neurons5  = [NEURON_COUNT[sp] for sp in sp5]
    noncan5   = [noncanonical_counts[sp] for sp in sp5]
    rho_nc, p_nc = spearman(neurons5, noncan5)

    # Correlation: combined vs neurons (5 species)
    combined5 = [canonical_counts[sp] + noncanonical_counts[sp] for sp in sp5]
    rho_cb, p_cb = spearman(neurons5, combined5)

    print(f"\n  Non-canonical RBP counts (5 species):")
    for sp in sp5:
        print(f"    {sp:35s}  canonical={canonical_counts[sp]:5d}  "
              f"non-canonical={noncanonical_counts[sp]:5d}")

    print(f"\n  Non-canonical vs neurons:  rho={rho_nc:.3f}  p={p_nc:.4f}  "
          f"{sig_label(p_nc)}  (expected: rho=0.900, p=0.037)")
    print(f"  Combined vs neurons:       rho={rho_cb:.3f}  p={p_cb:.4f}  "
          f"{sig_label(p_cb)}  (expected: rho=0.975, p=0.005)")

    results = []
    for sp in SPECIES_ORDER:
        results.append({
            "species":             sp,
            "neurons":             NEURON_COUNT[sp],
            "canonical_count":     canonical_counts[sp],
            "noncanonical_count":  noncanonical_counts[sp],
            "combined_count":      canonical_counts[sp] + noncanonical_counts[sp],
            "xenopus_excluded":    sp not in SPECIES_NONCANONICAL,
        })

    # Append correlation summary
    results.append({
        "species":            "CORRELATION: non-canonical vs neurons (n=5)",
        "neurons":            "",
        "canonical_count":    "",
        "noncanonical_count": round(rho_nc, 3),
        "combined_count":     round(p_nc, 4),
        "xenopus_excluded":   sig_label(p_nc),
    })
    results.append({
        "species":            "CORRELATION: combined vs neurons (n=5)",
        "neurons":            "",
        "canonical_count":    "",
        "noncanonical_count": round(rho_cb, 3),
        "combined_count":     round(p_cb, 4),
        "xenopus_excluded":   sig_label(p_cb),
    })

    return results


# ════════════════════════════════════════════════════════════════════════════
# 2. RBP–TF overlap
# ════════════════════════════════════════════════════════════════════════════

def run_rbp_tf_overlap(rbp_rows, tf_rows):
    """
    Intersect EuRBPDB gene symbols with AnimalTFDB gene symbols per species.

    Expected result (from manuscript):
      Overlap ranges 2.6–8.3% (excluding zebrafish where zf-C2H2 inflates to 16.0%)
      Human: 194/2961 = 6.6% overlap (86 zf-C2H2)
    """
    # Build gene symbol sets
    rbp_genes = defaultdict(set)
    tf_genes  = defaultdict(set)

    for r in rbp_rows:
        sp  = r.get("species", "").strip()
        sym = r.get("gene_symbol", "").strip().upper()
        if sp and sym:
            rbp_genes[sp].add(sym)

    for r in tf_rows:
        sp  = r.get("species", "").strip()
        sym = r.get("gene_symbol", "").strip().upper()
        if sp and sym:
            tf_genes[sp].add(sym)

    results = []
    print(f"\n  {'Species':35s} {'RBPs':>6} {'TFs':>6} {'Overlap':>8} {'%':>7}")
    for sp in SPECIES_ORDER:
        rbp_set  = rbp_genes[sp]
        tf_set   = tf_genes[sp]
        overlap  = rbp_set & tf_set
        n_rbp    = len(rbp_set)
        n_tf     = len(tf_set)
        n_ov     = len(overlap)
        pct      = 100.0 * n_ov / n_rbp if n_rbp > 0 else 0.0

        note = ""
        if sp == "Danio_rerio":
            note = "inflated by zf-C2H2 multi-annotation"

        print(f"  {sp:35s} {n_rbp:>6} {n_tf:>6} {n_ov:>8} {pct:>6.1f}%"
              + (f"  [{note}]" if note else ""))

        results.append({
            "species":           sp,
            "n_rbps":            n_rbp,
            "n_tfs":             n_tf,
            "n_overlap":         n_ov,
            "overlap_pct":       round(pct, 1),
            "note":              note,
        })

    return results


# ════════════════════════════════════════════════════════════════════════════
# 3. Vertebrate domain expansion
# ════════════════════════════════════════════════════════════════════════════

def run_domain_expansion(pfam_rows):
    """
    Compute mean gene count per Pfam domain family in invertebrates vs vertebrates.
    Identify domains with ≥3-fold enrichment or vertebrate-specific emergence.

    Expected (from manuscript Figure 3A):
      PARP domain:   ~8.9-fold
      YTH domain:    ~8.0-fold
      RAP domain:    ~11.0-fold
      RNase A domain: vertebrate-specific (mammalian)
      OAS domain:    vertebrate-specific
    """
    # Count genes per (family, species)
    family_species_count = defaultdict(lambda: defaultdict(int))

    for row in pfam_rows:
        sp     = row.get("species", "").strip()
        family = row.get("eurbpdb_family", "").strip()
        if sp not in SPECIES_ORDER:
            continue
        if not family or family == "Non-canonical":
            continue
        family_species_count[family][sp] += 1

    results = []
    for family, sp_counts in family_species_count.items():
        invert_counts = [sp_counts.get(sp, 0) for sp in INVERTEBRATES]
        vert_counts   = [sp_counts.get(sp, 0) for sp in VERTEBRATES]

        mean_invert = sum(invert_counts) / len(invert_counts)
        mean_vert   = sum(vert_counts)   / len(vert_counts)

        total_genes = sum(sp_counts.values())
        n_species   = len(sp_counts)

        # Expansion ratio
        if mean_invert == 0 and mean_vert > 0:
            ratio = float("inf")
            expansion_type = "vertebrate-specific"
        elif mean_invert == 0:
            ratio = 1.0
            expansion_type = "absent"
        else:
            ratio = mean_vert / mean_invert
            if ratio >= 3.0:
                expansion_type = "vertebrate-expanded (>=3x)"
            elif ratio >= 2.0:
                expansion_type = "vertebrate-enriched (2-3x)"
            else:
                expansion_type = "conserved"

        results.append({
            "pfam_family":      family,
            "n_species":        n_species,
            "total_genes":      total_genes,
            "mean_invertebrate": round(mean_invert, 2),
            "mean_vertebrate":  round(mean_vert,   2),
            "vert_invert_ratio": round(ratio, 2) if ratio != float("inf") else "inf",
            "expansion_type":   expansion_type,
            "worm":   sp_counts.get("Caenorhabditis_elegans",  0),
            "fly":    sp_counts.get("Drosophila_melanogaster", 0),
            "zfish":  sp_counts.get("Danio_rerio",             0),
            "frog":   sp_counts.get("Xenopus_tropicalis",      0),
            "mouse":  sp_counts.get("Mus_musculus",            0),
            "human":  sp_counts.get("Homo_sapiens",            0),
        })

    # Sort by ratio descending (inf first)
    def sort_key(r):
        v = r["vert_invert_ratio"]
        return -float("inf") if v == "inf" else -float(v)

    results.sort(key=sort_key)

    # Print top expanded
    print(f"\n  Top vertebrate-expanded Pfam domain families (≥3-fold or vertebrate-specific):")
    print(f"  {'Family':30s} {'ratio':>8}  {'type':35s}  "
          f"{'worm':>5} {'fly':>5} {'zfish':>6} {'frog':>5} {'mouse':>6} {'human':>6}")
    for r in results:
        if r["expansion_type"] in ("vertebrate-specific", "vertebrate-expanded (>=3x)"):
            ratio_str = str(r["vert_invert_ratio"])
            print(f"  {r['pfam_family']:30s} {ratio_str:>8}  "
                  f"{r['expansion_type']:35s}  "
                  f"{r['worm']:>5} {r['fly']:>5} {r['zfish']:>6} "
                  f"{r['frog']:>5} {r['mouse']:>6} {r['human']:>6}")

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
    print("Supplementary Analyses")
    print(f"  data_dir   : {data_dir}")
    print(f"  output_dir : {output_dir}")
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────────────
    print("\n[0] Loading data …")
    rbp_rows  = load_csv(os.path.join(data_dir, "new_rbp_db.csv"))
    tf_rows   = load_csv(os.path.join(data_dir, "new_tf_db.csv"))
    pfam_rows = load_tsv(os.path.join(data_dir, "eurbpdb_pfam_mapping.tsv"))

    # ── 1. Non-canonical ──────────────────────────────────────────────────
    print("\n[1] Non-canonical RBP analysis …")
    nc_results = run_noncanonical(rbp_rows)
    write_csv(
        os.path.join(output_dir, "noncanonical_correlations.csv"),
        nc_results,
        ["species", "neurons", "canonical_count",
         "noncanonical_count", "combined_count", "xenopus_excluded"],
    )

    # ── 2. RBP-TF overlap ────────────────────────────────────────────────
    print("\n[2] RBP–TF overlap analysis …")
    overlap_results = run_rbp_tf_overlap(rbp_rows, tf_rows)
    write_csv(
        os.path.join(output_dir, "rbp_tf_overlap.csv"),
        overlap_results,
        ["species", "n_rbps", "n_tfs", "n_overlap", "overlap_pct", "note"],
    )

    # ── 3. Domain expansion ───────────────────────────────────────────────
    print("\n[3] Vertebrate domain expansion analysis …")
    expansion_results = run_domain_expansion(pfam_rows)
    write_csv(
        os.path.join(output_dir, "domain_expansion.csv"),
        expansion_results,
        ["pfam_family", "n_species", "total_genes",
         "mean_invertebrate", "mean_vertebrate", "vert_invert_ratio",
         "expansion_type", "worm", "fly", "zfish", "frog", "mouse", "human"],
    )

    print("\nDone.")
    print("\nExpected results (from manuscript):")
    print("  Non-canonical vs neurons (n=5):  rho=0.900, p=0.037 *")
    print("  Combined vs neurons (n=5):        rho=0.975, p=0.005 **")
    print("  RBP-TF overlap: 2.6–8.3% (excl. zebrafish 16.0%)")
    print("  Human overlap: 194/2961 = 6.6%")
    print("  Top expanded domains: RAP (11x), PARP (8.9x), YTH (8x)")


if __name__ == "__main__":
    main()
