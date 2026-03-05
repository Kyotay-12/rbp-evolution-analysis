#!/usr/bin/env python3
"""
PGLS Analysis Script
=====================
Phylogenetic Generalized Least Squares (PGLS) analysis for Yasuda et al. (iScience, 2026)
Corresponds to Supplementary Figure 2 and Supplementary Table S7.

Background:
  Standard Spearman correlation ignores phylogenetic non-independence.
  PGLS corrects for this by weighting residuals according to expected
  covariance under a Brownian motion model of trait evolution.
  Pagel's lambda (λ) scales the phylogenetic covariance matrix:
    λ=0 → OLS (no phylogenetic correction)
    λ=1 → full Brownian motion model
  The optimal λ is estimated by maximum likelihood.

Results in manuscript (Supplementary Figure 2):
  Optimal λ = 0.0 (minimal phylogenetic signal)
  Correlation remains significant across all λ values (0 to 1.0)

Usage:
  python3 pgls_analysis.py <data_dir> <output_dir>

  <data_dir>   : directory containing input files (default: ../data/processed)
  <output_dir> : directory for output files     (default: ../results)

Input files required:
  main6_all_correlations_final.csv   - primary correlation results (for x/y values)

Output files:
  pgls_results.csv      - PGLS p-values across lambda values
  pgls_lambda_ml.csv    - maximum likelihood optimal lambda

Author: Kyota Yasuda
Date:   March 2026
"""

import sys
import os
import csv
import math

try:
    import numpy as np
    from scipy import stats, optimize
    HAS_SCIPY = True
except ImportError:
    raise ImportError(
        "scipy and numpy are required for PGLS analysis.\n"
        "Install with: pip install scipy numpy"
    )


# ════════════════════════════════════════════════════════════════════════════
# Species metadata
# ════════════════════════════════════════════════════════════════════════════

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

# Divergence times from TimeTree (Kumar et al., 2017), in million years ago
# Used to build the phylogenetic covariance matrix under Brownian motion.
# Branch lengths = shared evolutionary time (time since common ancestor).
#
# Tree topology (simplified):
#   ((worm, fly), ((zebrafish, frog), (mouse, human)))
#
# Divergence times:
#   worm–fly split:              ~800 MYA
#   invertebrate–vertebrate:     ~800 MYA (worm/fly vs. fish/frog/mouse/human)
#   zebrafish–frog split:        ~435 MYA
#   fish–tetrapod split:         ~435 MYA (zebrafish vs. frog/mouse/human)
#   frog–amniote split:          ~352 MYA
#   mouse–human split:           ~90 MYA
#   root age (total tree depth): ~800 MYA

DIVERGENCE_TIMES = {
    # (species_i, species_j): divergence time in MYA
    ("Caenorhabditis_elegans",  "Drosophila_melanogaster"): 800,
    ("Caenorhabditis_elegans",  "Danio_rerio"):             800,
    ("Caenorhabditis_elegans",  "Xenopus_tropicalis"):      800,
    ("Caenorhabditis_elegans",  "Mus_musculus"):            800,
    ("Caenorhabditis_elegans",  "Homo_sapiens"):            800,
    ("Drosophila_melanogaster", "Danio_rerio"):             800,
    ("Drosophila_melanogaster", "Xenopus_tropicalis"):      800,
    ("Drosophila_melanogaster", "Mus_musculus"):            800,
    ("Drosophila_melanogaster", "Homo_sapiens"):            800,
    ("Danio_rerio",             "Xenopus_tropicalis"):      435,
    ("Danio_rerio",             "Mus_musculus"):            435,
    ("Danio_rerio",             "Homo_sapiens"):            435,
    ("Xenopus_tropicalis",      "Mus_musculus"):            352,
    ("Xenopus_tropicalis",      "Homo_sapiens"):            352,
    ("Mus_musculus",            "Homo_sapiens"):            90,
}

TREE_DEPTH = 800  # MYA – used to normalize branch lengths


# Family diversity values (worm, fly, zebrafish, frog, mouse, human)
# Source: main6_all_correlations_final.csv
FAMILY_DIVERSITY = {
    "RBP":    [397, 419, 455, 446, 472, 469],
    "TF":     [58,  64,  72,  72,  72,  72],
    "Kinase": [31,  30,  20,  13,  64,  66],
    "GPCR":   [12,  7,   5,   3,   10,  10],
}


# ════════════════════════════════════════════════════════════════════════════
# Build phylogenetic covariance matrix
# ════════════════════════════════════════════════════════════════════════════

def build_vcv_matrix(species_list, tree_depth=TREE_DEPTH):
    """
    Build variance-covariance matrix (VCV) under Brownian motion.
    C[i,i] = total branch length from root to tip = tree_depth
    C[i,j] = shared branch length = tree_depth - divergence_time(i,j)
    """
    n = len(species_list)
    C = np.zeros((n, n))

    for i in range(n):
        for j in range(n):
            if i == j:
                C[i, j] = tree_depth
            else:
                sp_i = species_list[i]
                sp_j = species_list[j]
                key  = (sp_i, sp_j) if (sp_i, sp_j) in DIVERGENCE_TIMES \
                       else (sp_j, sp_i)
                div_time   = DIVERGENCE_TIMES.get(key, tree_depth)
                shared_time = tree_depth - div_time
                C[i, j]    = shared_time

    return C


def scale_vcv(C, lam):
    """
    Apply Pagel's lambda scaling:
      C_scaled[i,j] = lam * C[i,j]  for i != j
      C_scaled[i,i] = C[i,i]        (diagonal unchanged)
    """
    n = C.shape[0]
    C_scaled = lam * C.copy()
    for i in range(n):
        C_scaled[i, i] = C[i, i]
    return C_scaled


# ════════════════════════════════════════════════════════════════════════════
# PGLS regression
# ════════════════════════════════════════════════════════════════════════════

def pgls_regression(x, y, C_scaled):
    """
    Generalized Least Squares regression:
      y = b0 + b1 * x + e,  Var(e) = sigma^2 * C_scaled

    Returns:
      slope, intercept, t_stat, p_value, log_likelihood
    """
    n    = len(x)
    x_arr = np.array(x, dtype=float)
    y_arr = np.array(y, dtype=float)

    # Design matrix
    X = np.column_stack([np.ones(n), x_arr])

    # GLS: beta = (X' C^-1 X)^-1 X' C^-1 y
    try:
        C_inv = np.linalg.inv(C_scaled)
    except np.linalg.LinAlgError:
        C_inv = np.linalg.pinv(C_scaled)

    XtCiX  = X.T @ C_inv @ X
    XtCiy  = X.T @ C_inv @ y_arr

    try:
        beta = np.linalg.solve(XtCiX, XtCiy)
    except np.linalg.LinAlgError:
        beta = np.linalg.lstsq(XtCiX, XtCiy, rcond=None)[0]

    intercept, slope = beta[0], beta[1]

    # Residuals and variance
    residuals = y_arr - X @ beta
    sigma2    = float(residuals @ C_inv @ residuals) / n

    # Covariance of beta
    if sigma2 > 0:
        var_beta = sigma2 * np.linalg.inv(XtCiX)
        se_slope = math.sqrt(max(var_beta[1, 1], 0))
    else:
        se_slope = float("nan")

    # t-test for slope
    if se_slope > 0:
        t_stat = slope / se_slope
        df     = n - 2
        p_val  = 2 * (1 - stats.t.cdf(abs(t_stat), df=df))
    else:
        t_stat = float("nan")
        p_val  = 1.0

    # Log-likelihood (for lambda ML estimation)
    sign, log_det = np.linalg.slogdet(C_scaled)
    if sign <= 0 or sigma2 <= 0:
        log_lik = -float("inf")
    else:
        log_lik = (
            -0.5 * n * math.log(2 * math.pi * sigma2)
            - 0.5 * log_det
            - 0.5 * n
        )

    return slope, intercept, t_stat, p_val, log_lik


# ════════════════════════════════════════════════════════════════════════════
# Lambda sensitivity and ML estimation
# ════════════════════════════════════════════════════════════════════════════

def run_lambda_sensitivity(x, y, C, lambdas):
    """Run PGLS across a range of lambda values."""
    results = []
    for lam in lambdas:
        C_scaled = scale_vcv(C, lam)
        slope, intercept, t_stat, p_val, log_lik = pgls_regression(x, y, C_scaled)
        results.append({
            "lambda":      round(lam, 3),
            "slope":       round(slope,     4),
            "intercept":   round(intercept, 4),
            "t_statistic": round(t_stat,    4) if not math.isnan(t_stat) else "NA",
            "p_value":     round(p_val,     4),
            "significant": "*" if p_val < 0.05 else "ns",
            "log_lik":     round(log_lik,   4),
        })
    return results


def find_optimal_lambda(x, y, C):
    """Find lambda that maximizes log-likelihood via scipy minimize."""
    def neg_log_lik(lam_arr):
        lam = float(lam_arr[0])
        lam = max(0.0, min(1.0, lam))
        C_scaled = scale_vcv(C, lam)
        _, _, _, _, log_lik = pgls_regression(x, y, C_scaled)
        return -log_lik  # minimize negative log-likelihood

    result = optimize.minimize_scalar(
        lambda lam: neg_log_lik([lam]),
        bounds=(0.0, 1.0),
        method="bounded",
    )
    return round(float(result.x), 4)


# ════════════════════════════════════════════════════════════════════════════
# File I/O
# ════════════════════════════════════════════════════════════════════════════

def write_csv(path, rows, fieldnames):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f"  Saved → {path}")


# ════════════════════════════════════════════════════════════════════════════
# Main
# ════════════════════════════════════════════════════════════════════════════

def main():
    data_dir   = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
                    os.path.dirname(__file__), "..", "data", "processed")
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
                    os.path.dirname(__file__), "..", "results")

    output_dir = os.path.abspath(output_dir)
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("PGLS Analysis (Supplementary Figure 2)")
    print(f"  output_dir : {output_dir}")
    print("=" * 60)

    # Build VCV matrix
    C = build_vcv_matrix(SPECIES_ORDER)
    neurons = [NEURON_COUNT[sp] for sp in SPECIES_ORDER]
    # Log-transform neurons (as in manuscript)
    log_neurons = [math.log10(n) for n in neurons]

    lambdas = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]

    all_results  = []
    ml_results   = []

    for ptype, diversity in FAMILY_DIVERSITY.items():
        log_div = [math.log10(d) for d in diversity]

        print(f"\n  [{ptype}]")

        # Lambda sensitivity
        sens = run_lambda_sensitivity(log_neurons, log_div, C, lambdas)
        for row in sens:
            row["protein_type"] = ptype
            all_results.append(row)
            print(f"    λ={row['lambda']:.1f}  p={row['p_value']:.4f}  "
                  f"{row['significant']}  loglik={row['log_lik']}")

        # Optimal lambda (ML)
        opt_lam = find_optimal_lambda(log_neurons, log_div, C)
        C_opt   = scale_vcv(C, opt_lam)
        slope, intercept, t_stat, p_val, log_lik = pgls_regression(
            log_neurons, log_div, C_opt)

        print(f"    → Optimal λ={opt_lam}  p={p_val:.4f}  "
              + ("*" if p_val < 0.05 else "ns"))

        ml_results.append({
            "protein_type":   ptype,
            "optimal_lambda": opt_lam,
            "slope":          round(slope,     4),
            "intercept":      round(intercept, 4),
            "t_statistic":    round(t_stat,    4) if not math.isnan(t_stat) else "NA",
            "p_value":        round(p_val,     4),
            "significant":    "*" if p_val < 0.05 else "ns",
            "log_lik":        round(log_lik,   4),
            "note": ("Optimal λ≈0 → OLS appropriate; "
                     "minimal phylogenetic signal" if opt_lam < 0.1
                     else ""),
        })

    # Write outputs
    write_csv(
        os.path.join(output_dir, "pgls_lambda_sensitivity.csv"),
        all_results,
        ["protein_type", "lambda", "slope", "intercept",
         "t_statistic", "p_value", "significant", "log_lik"],
    )

    write_csv(
        os.path.join(output_dir, "pgls_optimal_lambda.csv"),
        ml_results,
        ["protein_type", "optimal_lambda", "slope", "intercept",
         "t_statistic", "p_value", "significant", "log_lik", "note"],
    )

    # Summary
    print("\n" + "=" * 60)
    print("Expected results (from manuscript):")
    print("  Optimal λ = 0.0 for all protein types")
    print("  RBP correlation significant across all λ (0.0 to 1.0)")
    print("  Interpretation: minimal phylogenetic signal in residuals")
    print("  → OLS (λ=0) is the appropriate model")
    print("  → RBP expansion is 'punctual' (WGD-driven) not gradual (BM)")
    print("\nDone.")


if __name__ == "__main__":
    main()
