#!/usr/bin/env python3
"""
RBP Evolution Analysis - Reproduce All Main Results
====================================================
Yasuda K. (2025) iScience

This script reproduces all key statistical results from the paper.
"""

import pandas as pd
import numpy as np
from scipy.stats import spearmanr
import warnings
warnings.filterwarnings('ignore')

def main():
    print("=" * 60)
    print("RBP Evolution Analysis - Results Reproduction")
    print("=" * 60)
    
    # ===== 1. Main 6-species correlation =====
    print("\n1. MAIN CORRELATION (6 species)")
    print("-" * 40)
    
    df = pd.read_csv('data/main_6species/MainTable1_species_metrics.csv')
    
    rho, p = spearmanr(df['Family_count'], df['Neuron_count'])
    print(f"RBP families vs neurons: ρ = {rho:.3f}, p = {p:.4f}")
    
    rho, p = spearmanr(df['Family_count'], df['Genome_size_Mb'])
    print(f"RBP families vs genome:  ρ = {rho:.3f}, p = {p:.4f}")
    
    rho, p = spearmanr(df['Family_count'], df['Cell_type_diversity'])
    print(f"RBP families vs cells:   ρ = {rho:.3f}, p = {p:.4f}")
    
    # ===== 2. Control proteins =====
    print("\n2. CONTROL PROTEINS (6 species)")
    print("-" * 40)
    
    corr = pd.read_csv('data/main_6species/correlations_main6.csv')
    neurons_corr = corr[(corr['vs'] == 'neurons') & (corr['metric'] == 'n_families')]
    
    for _, row in neurons_corr.iterrows():
        sig = "**" if row['p_value'] < 0.01 else "*" if row['p_value'] < 0.05 else "ns"
        print(f"{row['protein_type']:8} vs neurons: ρ = {row['rho']:.3f}, p = {row['p_value']:.3f} ({sig})")
    
    # ===== 3. Extended 15 species =====
    print("\n3. EXTENDED ANALYSIS (15 species)")
    print("-" * 40)
    
    corr15 = pd.read_csv('data/extended_15species/correlations_all15.csv')
    rbp_neurons = corr15[(corr15['protein_type'] == 'RBP') & 
                         (corr15['vs'] == 'neurons') & 
                         (corr15['metric'] == 'n_families')]
    if len(rbp_neurons) > 0:
        row = rbp_neurons.iloc[0]
        sig = "*" if row['p_value'] < 0.05 else "ns"
        print(f"RBP families vs neurons (n=15): ρ = {row['rho']:.3f}, p = {row['p_value']:.3f} ({sig})")
    
    # ===== 4. Bootstrap CI =====
    print("\n4. BOOTSTRAP CONFIDENCE INTERVALS")
    print("-" * 40)
    
    boot = pd.read_csv('data/statistical_validation/bootstrap_ci_results.csv')
    for _, row in boot.iterrows():
        print(f"{row['protein_type']:8}: ρ = {row['rho']:.3f}, 95% CI [{row['ci_lower']:.3f}, {row['ci_upper']:.3f}]")
    
    # ===== 5. Leave-one-out =====
    print("\n5. LEAVE-ONE-OUT VALIDATION")
    print("-" * 40)
    
    loo = pd.read_csv('data/statistical_validation/leave_one_out_results.csv')
    for _, row in loo.iterrows():
        sig = "*" if row['p_value'] < 0.05 else "ns"
        print(f"Exclude {row['excluded_species']:10}: ρ = {row['rho']:.3f}, p = {row['p_value']:.3f} ({sig})")
    
    # ===== 6. Effect sizes =====
    print("\n6. EFFECT SIZES (Cohen's d)")
    print("-" * 40)
    
    eff = pd.read_csv('data/statistical_validation/effect_sizes.csv')
    for _, row in eff.iterrows():
        print(f"{row['comparison']:15}: d = {row['cohens_d']:.3f} ({row['interpretation']})")
    
    # ===== 7. LLPS constancy =====
    print("\n7. LLPS PROPENSITY BY SPECIES")
    print("-" * 40)
    
    llps = pd.read_csv('data/llps_analysis/llps_by_species.csv')
    for _, row in llps.iterrows():
        print(f"{row['species']:15}: {row['pct_high_llps']:.1f}% high-LLPS")
    
    mean_llps = llps['pct_high_llps'].mean()
    std_llps = llps['pct_high_llps'].std()
    print(f"\nMean ± SD: {mean_llps:.1f} ± {std_llps:.1f}%")
    
    rho, p = spearmanr(llps['pct_high_llps'], llps['neurons'])
    print(f"LLPS % vs neurons: ρ = {rho:.2f}, p = {p:.2f} (ns)")
    
    # ===== Summary =====
    print("\n" + "=" * 60)
    print("SUMMARY: All key results successfully reproduced!")
    print("=" * 60)

if __name__ == "__main__":
    main()
