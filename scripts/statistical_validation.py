#!/usr/bin/env python3
"""
Statistical Validation Script
==============================
Reproduces all statistical validation results in Yasuda et al. (iScience, 2026)

Analyses performed:
  1. Bootstrap resampling (n=1,000)  → 95% CI, median rho, % negative
  2. Leave-one-out (LOO) validation  → rho per excluded species
  3. Cohen's d effect sizes          → RBP vs each control protein type
  4. Pfam clan-level robustness      → re-run correlation with clan-level diversity

Usage:
  python3 statistical_validation.py <data_dir> <output_dir>

  <data_dir>   : directory containing input files (default: ../data/processed)
  <output_dir> : directory for output files     (default: ../results)

Input files required:
  main6_all_correlations_final.csv   - primary correlation results
  eurbpdb_pfam_mapping.tsv           - gene-level Pfam assignments (for LOO/bootstrap)
  new_rbp_db.csv                     - RBP database
  new_tf_db.csv                      - TF database
  controls_pfam_mapping.tsv          - Kinase/GPCR Pfam assignments
  Pfam-A_clans_decoded.tsv           - Pfam clan mappings (for clan robustness)

Output files:
  bootstrap_results.csv              - bootstrap CI per protein type
  loo_results.csv                    - leave-one-out rho per excluded species
  cohens_d_results.csv               - effect sizes RBP vs controls
  clan_robustness.csv                - clan-level correlation results

Author: Kyota Yasuda
Date:   March 2026
"""

import sys
import os
import csv
import math
import random
from collections import defaultdict

try:
    import numpy as np
    from scipy import stats
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False
    print("[WARN] scipy/numpy not found – using built-in implementations")


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

# Family diversity values from main6_all_correlations_final.csv
# (worm, fly, zebrafish, frog, mouse, human)
FAMILY_DIVERSITY = {
    "RBP":    [397, 419, 455, 446, 472, 469],
    "TF":     [58,  64,  72,  72,  72,  72],
    "Kinase": [31,  30,  20,  13,  64,  66],
    "GPCR":   [12,  7,   5,   3,   10,  10],
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
# 1. Bootstrap resampling
# ════════════════════════════════════════════════════════════════════════════

def run_bootstrap(n_iter=1000, seed=42):
    """
    Bootstrap resampling (with replacement) of 6 species.
    For each protein type, resample species indices n_iter times
    and compute Spearman rho on family diversity vs neuron count.

    Returns list of result dicts.
    """
    random.seed(seed)
    neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]

    results = []
    for ptype, diversity in FAMILY_DIVERSITY.items():
        rhos = []
        for _ in range(n_iter):
            idx = [random.randint(0, 5) for _ in range(6)]
            x   = [neurons[i]   for i in idx]
            y   = [diversity[i] for i in idx]
            # skip if all x or all y identical (degenerate)
            if len(set(x)) < 2 or len(set(y)) < 2:
                continue
            rho, _ = spearman(x, y)
            rhos.append(rho)

        rhos.sort()
        n     = len(rhos)
        ci_lo = rhos[int(0.025 * n)]
        ci_hi = rhos[int(0.975 * n)]
        med   = rhos[n // 2]
        pct_neg = 100.0 * sum(1 for r in rhos if r < 0) / n

        results.append({
            "protein_type":      ptype,
            "n_iterations":      n,
            "median_rho":        round(med,    3),
            "ci_lower_95":       round(ci_lo,  3),
            "ci_upper_95":       round(ci_hi,  3),
            "pct_negative":      round(pct_neg, 1),
            "observed_rho":      FAMILY_DIVERSITY[ptype],   # stored as list, overwritten below
        })
        # overwrite with scalar observed rho
        obs_rho, _ = spearman(neurons, diversity)
        results[-1]["observed_rho"] = round(obs_rho, 3)

        print(f"  {ptype:6s}  median={med:+.3f}  95%CI=[{ci_lo:+.3f}, {ci_hi:+.3f}]  "
              f"neg={pct_neg:.1f}%")

    return results


# ════════════════════════════════════════════════════════════════════════════
# 2. Leave-one-out (LOO) validation
# ════════════════════════════════════════════════════════════════════════════

def run_loo():
    """
    Remove each species in turn and recompute Spearman rho.
    Returns list of result dicts.
    """
    neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]
    results = []

    for ptype, diversity in FAMILY_DIVERSITY.items():
        for excl_idx, excl_sp in enumerate(SPECIES_ORDER):
            x = [neurons[i]   for i in range(6) if i != excl_idx]
            y = [diversity[i] for i in range(6) if i != excl_idx]
            rho, pval = spearman(x, y)

            def sig(p):
                if p < 0.001: return "***"
                if p < 0.01:  return "**"
                if p < 0.05:  return "*"
                return "ns"

            results.append({
                "protein_type":    ptype,
                "excluded_species": excl_sp,
                "n":               5,
                "rho":             round(rho,  3),
                "p_value":         round(pval, 4),
                "significance":    sig(pval),
            })

    # Print RBP summary
    print("\n  LOO results for RBP family diversity:")
    for r in results:
        if r["protein_type"] == "RBP":
            print(f"    excl. {r['excluded_species']:35s}  "
                  f"rho={r['rho']:+.3f}  p={r['p_value']:.4f}  {r['significance']}")

    return results


# ════════════════════════════════════════════════════════════════════════════
# 3. Cohen's d effect sizes
# ════════════════════════════════════════════════════════════════════════════

def run_cohens_d(n_iter=1000, seed=42):
    """
    Cohen's d = (rho_RBP - rho_control) / pooled_SD
    where SD is estimated from bootstrap distributions.

    Returns list of result dicts.
    """
    random.seed(seed)
    neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]

    # Generate bootstrap rho distributions for all protein types
    boot_rhos = {ptype: [] for ptype in FAMILY_DIVERSITY}
    for _ in range(n_iter):
        idx = [random.randint(0, 5) for _ in range(6)]
        x   = [neurons[i] for i in idx]
        if len(set(x)) < 2:
            continue
        for ptype, diversity in FAMILY_DIVERSITY.items():
            y = [diversity[i] for i in idx]
            if len(set(y)) < 2:
                continue
            rho, _ = spearman(x, y)
            boot_rhos[ptype].append(rho)

    def sd(lst):
        if len(lst) < 2:
            return float("nan")
        m = sum(lst) / len(lst)
        return math.sqrt(sum((v - m) ** 2 for v in lst) / (len(lst) - 1))

    rbp_rhos  = boot_rhos["RBP"]
    rbp_mean  = sum(rbp_rhos) / len(rbp_rhos)
    rbp_sd    = sd(rbp_rhos)

    results = []
    for ctrl in ["TF", "Kinase", "GPCR"]:
        ctrl_rhos = boot_rhos[ctrl]
        ctrl_mean = sum(ctrl_rhos) / len(ctrl_rhos)
        ctrl_sd   = sd(ctrl_rhos)
        pooled_sd = math.sqrt((rbp_sd ** 2 + ctrl_sd ** 2) / 2)
        d         = (rbp_mean - ctrl_mean) / pooled_sd if pooled_sd > 0 else float("nan")

        def interp(d_abs):
            if d_abs >= 0.8: return "large"
            if d_abs >= 0.5: return "medium"
            if d_abs >= 0.2: return "small"
            return "negligible"

        results.append({
            "comparison":       f"RBP vs {ctrl}",
            "rbp_mean_rho":     round(rbp_mean,  3),
            "ctrl_mean_rho":    round(ctrl_mean, 3),
            "rbp_sd":           round(rbp_sd,    3),
            "ctrl_sd":          round(ctrl_sd,   3),
            "pooled_sd":        round(pooled_sd, 3),
            "cohens_d":         round(d,         3),
            "interpretation":   interp(abs(d)),
        })
        print(f"  RBP vs {ctrl:6s}  d={d:+.3f}  ({interp(abs(d))})")

    return results


# ════════════════════════════════════════════════════════════════════════════
# 4. Pfam clan-level robustness
# ════════════════════════════════════════════════════════════════════════════

def run_clan_robustness(data_dir):
    """
    Re-compute RBP family diversity at Pfam clan level using
    Pfam-A_clans_decoded.tsv (clan_id → pfam_id mapping).
    Then run Spearman vs neurons.

    If clan file not found, reports that and skips gracefully.
    """
    clan_file = os.path.join(data_dir, "Pfam-A_clans_decoded.tsv")
    pfam_file = os.path.join(data_dir, "eurbpdb_pfam_mapping.tsv")

    if not os.path.isfile(clan_file):
        print(f"  [SKIP] Pfam-A_clans_decoded.tsv not found in {data_dir}")
        print(f"         Download from: https://ftp.ebi.ac.uk/pub/databases/Pfam/"
              f"current_release/Pfam-A.clans.tsv.gz")
        return []

    if not os.path.isfile(pfam_file):
        raise FileNotFoundError(f"Required: {pfam_file}")

    # Load clan mapping: pfam_id → clan_id (use pfam_id itself if no clan)
    pfam_to_clan = {}
    with open(clan_file, encoding="utf-8") as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                pfam_id = parts[0].strip()
                clan_id = parts[1].strip()
                pfam_to_clan[pfam_id] = clan_id if clan_id else pfam_id

    # Load gene-level Pfam assignments
    pfam_rows = load_tsv(pfam_file)

    # Count distinct clans per species
    species_clans = defaultdict(set)
    for row in pfam_rows:
        sp     = row.get("species", "").strip()
        family = row.get("eurbpdb_family", "").strip()
        pfams  = row.get("pfam_families", "").strip()
        if sp not in SPECIES_ORDER:
            continue
        if family and family != "Non-canonical":
            # Map each pfam to its clan
            for pfam_id in pfams.split(";"):
                pfam_id = pfam_id.strip()
                if pfam_id:
                    clan = pfam_to_clan.get(pfam_id, pfam_id)
                    species_clans[sp].add(clan)

    neurons      = [NEURON_COUNT[sp]      for sp in SPECIES_ORDER]
    clan_counts  = [len(species_clans[sp]) for sp in SPECIES_ORDER]
    rho, pval    = spearman(neurons, clan_counts)

    def sig(p):
        if p < 0.001: return "***"
        if p < 0.01:  return "**"
        if p < 0.05:  return "*"
        return "ns"

    print(f"  Clan-level diversity: {clan_counts}")
    print(f"  Spearman rho={rho:.3f}  p={pval:.4f}  {sig(pval)}")

    return [{
        "level":          "clan",
        "clan_counts":    "→".join(str(c) for c in clan_counts),
        "rho":            round(rho,  3),
        "p_value":        round(pval, 4),
        "significance":   sig(pval),
        "note":           "Should match domain-level rho=0.886 if robust",
    }]


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
    print("Statistical Validation")
    print(f"  data_dir   : {data_dir}")
    print(f"  output_dir : {output_dir}")
    print("=" * 60)

    # ── 1. Bootstrap ────────────────────────────────────────────────────────
    print("\n[1] Bootstrap resampling (n=1,000) …")
    boot = run_bootstrap(n_iter=1000)
    write_csv(
        os.path.join(output_dir, "bootstrap_results.csv"),
        boot,
        ["protein_type", "n_iterations", "observed_rho",
         "median_rho", "ci_lower_95", "ci_upper_95", "pct_negative"],
    )

    # ── 2. LOO ──────────────────────────────────────────────────────────────
    print("\n[2] Leave-one-out validation …")
    loo = run_loo()
    write_csv(
        os.path.join(output_dir, "loo_results.csv"),
        loo,
        ["protein_type", "excluded_species", "n", "rho", "p_value", "significance"],
    )

    # ── 3. Cohen's d ────────────────────────────────────────────────────────
    print("\n[3] Cohen's d effect sizes …")
    cohens = run_cohens_d(n_iter=1000)
    write_csv(
        os.path.join(output_dir, "cohens_d_results.csv"),
        cohens,
        ["comparison", "rbp_mean_rho", "ctrl_mean_rho",
         "rbp_sd", "ctrl_sd", "pooled_sd", "cohens_d", "interpretation"],
    )

    # ── 4. Pfam clan robustness ─────────────────────────────────────────────
    print("\n[4] Pfam clan-level robustness …")
    clan = run_clan_robustness(data_dir)
    if clan:
        write_csv(
            os.path.join(output_dir, "clan_robustness.csv"),
            clan,
            ["level", "clan_counts", "rho", "p_value", "significance", "note"],
        )

    # ── Summary ─────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("Summary (RBP family diversity vs neurons)")
    print("=" * 60)
    rbp_boot = next(r for r in boot if r["protein_type"] == "RBP")
    print(f"  Observed rho      : {rbp_boot['observed_rho']}")
    print(f"  Bootstrap median  : {rbp_boot['median_rho']}")
    print(f"  95% CI            : [{rbp_boot['ci_lower_95']}, {rbp_boot['ci_upper_95']}]")
    print(f"  % negative        : {rbp_boot['pct_negative']}%")
    rbp_loo = [r for r in loo if r["protein_type"] == "RBP"]
    rho_range = f"{min(r['rho'] for r in rbp_loo):.3f}–{max(r['rho'] for r in rbp_loo):.3f}"
    print(f"  LOO rho range     : {rho_range}")
    for c in cohens:
        print(f"  {c['comparison']:20s}  Cohen's d = {c['cohens_d']}")
    print("\nDone.")


if __name__ == "__main__":
    main()
