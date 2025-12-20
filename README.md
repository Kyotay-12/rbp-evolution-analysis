# RBP Evolution Analysis

**RNA-binding protein family diversification correlates with neural complexity across metazoan evolution**

Yasuda K. (2025) *iScience* (submitted)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17958862.svg)](https://doi.org/10.5281/zenodo.17958862)

---

## Overview

This repository contains all data and code to reproduce the analyses in Yasuda (2025). We demonstrate that RNA-binding protein (RBP) family diversity shows a strong, specific correlation with neural complexity across metazoan evolution (Spearman ρ = 0.943, *p* = 0.0048).

### Key Findings

- **RBP family diversity** increases ~3.8-fold from invertebrates (11 families) to humans (42 families)
- **Strong correlation** with neuronal count (ρ = 0.943), genome size (ρ = 1.000), and cell type diversity (ρ = 1.000)
- **RBP-specific**: Transcription factors, kinases, and GPCRs show weaker, non-significant correlations
- **LLPS propensity** remains constant (~14%) despite RBP expansion, validated by 6 independent predictors
- **Vertebrate-expanded families** (hnRNP, CELF, CPEB) are enriched for neural/synaptic functions

---

## Repository Structure

```
rbp-evolution-analysis/
├── data/
│   ├── main_6species/           # Primary 6-species analysis
│   │   ├── MainTable1_species_metrics.csv
│   │   ├── MainTable2_correlations.csv
│   │   ├── summary_statistics_main6.csv
│   │   ├── correlations_main6.csv
│   │   ├── family_diversity_6species.csv
│   │   ├── 3utr_rbp_summary_by_species.csv
│   │   └── 3utr_correlation_results.csv
│   │
│   ├── extended_15species/      # Extended 15-species validation
│   │   ├── family_diversity_15species_unified.csv
│   │   ├── correlations_all15.csv
│   │   ├── summary_statistics_all15.csv
│   │   ├── integrated_rbps_all15species_final.csv
│   │   └── integrated_controls_all15species_final.csv
│   │
│   ├── llps_analysis/           # LLPS prediction data
│   │   ├── llps_by_species.csv
│   │   ├── all_rbps_with_llphyscore.csv
│   │   ├── phasepred_high_llps_by_species.csv
│   │   └── [species]_rbp_phasepred_scores.csv (5 files)
│   │
│   ├── statistical_validation/  # Bootstrap, LOO, effect sizes
│   │   ├── bootstrap_ci_results.csv
│   │   ├── leave_one_out_results.csv
│   │   ├── effect_sizes.csv
│   │   └── fishers_z_results.csv
│   │
│   └── family_analysis/         # Family-level expansion
│       ├── vertebrate_expansion.csv
│       ├── family_level_correlations_42family_categorized.csv
│       └── orthogroup_expansion_analysis.csv
│
├── figures/                     # Publication figures (PNG/PDF)
├── scripts/                     # Analysis scripts (Python)
├── docs/                        # Additional documentation
└── README.md
```

---

## Data Description

### Main Tables (Paper Tables 1-2)

| File | Description |
|------|-------------|
| `MainTable1_species_metrics.csv` | Species-level RBP counts, family counts, complexity metrics |
| `MainTable2_correlations.csv` | Spearman correlations for all protein types vs complexity |

### Primary Dataset (6 Species)

| Species | Common Name | RBP Families | Neurons | Genome (Mb) |
|---------|-------------|--------------|---------|-------------|
| C. elegans | Worm | 11 | 302 | 100 |
| D. melanogaster | Fly | 13 | 100,000 | 143 |
| D. rerio | Zebrafish | 25 | 10,000,000 | 1,674 |
| X. tropicalis | Frog | 20 | 20,000,000 | 1,510 |
| M. musculus | Mouse | 40 | 71,000,000 | 2,700 |
| H. sapiens | Human | 42 | 86,000,000,000 | 3,200 |

### Key Statistics

| Analysis | Value |
|----------|-------|
| RBP families vs neurons (6 sp.) | ρ = 0.943, *p* = 0.0048 |
| RBP families vs neurons (15 sp.) | ρ = 0.581, *p* = 0.023 |
| Bootstrap 95% CI | [0.515, 1.000] |
| TF families vs neurons | ρ = 0.696, *p* = 0.125 (ns) |
| Kinase families vs neurons | ρ = 0.618, *p* = 0.191 (ns) |
| GPCR families vs neurons | ρ = 0.802, *p* = 0.055 (ns) |

---

## LLPS Predictors Used

| Predictor | Reference | Threshold |
|-----------|-----------|-----------|
| LLPhyScore | You et al., 2020 | >0.5 |
| catGRANULE | Mitchell et al., 2013 | rank >0.5 |
| PScore | Vernon et al., 2018 | rank >0.5 |
| PLAAC | Lancaster et al., 2014 | NLLR >0 |
| PhaSePred SaPS | Wang et al., 2022 | >0.5 |

All predictors show consistent results: LLPS propensity does not correlate with neural complexity (all *p* > 0.05).

---

## Reproducing Results

### Requirements

```bash
pip install pandas numpy scipy matplotlib seaborn
```

### Quick Start

```python
import pandas as pd
from scipy.stats import spearmanr

# Load main data
df = pd.read_csv('data/main_6species/MainTable1_species_metrics.csv')

# Reproduce main correlation
rho, p = spearmanr(df['Family_count'], df['Neuron_count'])
print(f"RBP families vs neurons: ρ = {rho:.3f}, p = {p:.4f}")
# Output: ρ = 0.943, p = 0.0048
```

### Full Analysis Pipeline

```bash
# Run all analyses
python scripts/run_all_analyses.py

# Generate figures
python scripts/generate_figures.py
```

---

## File Formats

All data files are CSV format with UTF-8 encoding. Key column definitions:

### MainTable1_species_metrics.csv
- `species`: Scientific name abbreviation
- `common_name`: Common name
- `RBP_count`: Total number of RBPs
- `Family_count`: Number of RBP families (42-family classification)
- `Neuron_count`: Estimated neuron number
- `Genome_size_Mb`: Genome size in megabases
- `Cell_type_diversity`: Estimated cell types
- `Mean_LLPS_score`: Mean LLPhyScore
- `High_LLPS_count`: RBPs with LLPhyScore >0.5
- `High_LLPS_percent`: Percentage high-LLPS

### correlations_main6.csv
- `protein_type`: RBP, TF, Kinase, or GPCR
- `metric`: n_families
- `vs`: Complexity metric (neurons, genome_mb, cell_types)
- `rho`: Spearman correlation coefficient
- `p_value`: Two-tailed p-value
- `n`: Sample size

---

## Citation

If you use this data or code, please cite:

```bibtex
@article{yasuda2025rbp,
  title={RNA-binding protein family diversification correlates with neural 
         complexity across metazoan evolution},
  author={Yasuda, Kyota},
  journal={iScience},
  year={2025},
  doi={10.5281/zenodo.17958862}
}
```

---

## Contact

Kyota Yasuda  
Graduate School of Integrated Sciences for Life  
Hiroshima University  
Email: kyotay12@hiroshima-u.ac.jp

---

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Changelog

- **v1.0.0** (2025-12-20): Initial release with complete dataset for iScience submission
