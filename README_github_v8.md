# RBP Family Diversification and Neural Complexity

**Repository for:** Yasuda K. "RNA-Binding Protein Family Diversification Correlates with Neural Complexity Across Metazoan Evolution" (iScience, 2026)

---

## Overview

This repository contains all analysis scripts and processed data for a systematic comparative analysis of RNA-binding protein (RBP) family diversity across six metazoan species. The central finding is that RBP Pfam domain family diversity correlates strongly with neural complexity (Spearman's ρ = 0.886, p = 0.019), a relationship that is RBP-specific and robust to multiple statistical and phylogenetic corrections.

**Species analyzed (primary):** *C. elegans*, *D. melanogaster*, *D. rerio*, *X. tropicalis*, *M. musculus*, *H. sapiens*

**Archived data:** Zenodo DOI: 10.5281/zenodo.18002634

---

## Key Results

| Analysis | Result |
|---|---|
| RBP family diversity vs neurons | ρ = 0.886, p = 0.019 ★ |
| TF family diversity vs neurons | ρ = 0.845 (vertebrate saturation at 72 families) |
| Kinase family diversity vs neurons | ρ = 0.429, ns |
| GPCR family diversity vs neurons | ρ = −0.116, ns |
| IDR longest length vs neurons | ρ = 0.829, p = 0.042 ★ |
| LLPS propensity vs neurons | ρ = 0.714, ns (conserved ~14–17%) |
| 3'UTR length vs neurons | ρ = 0.943, p = 0.0048 ★★ |
| Non-canonical RBPs vs neurons | ρ = 0.900, p = 0.037 ★ |

---

## Data Sources

### Primary RBP Data
- **EuRBPDB** (Li et al., 2020): RBP annotations for all 6 primary species  
  Files: `*_RBP.txt` format  
  Counts: worm 1,442 / fly 1,633 / zebrafish 2,349 / frog 1,758 / mouse 2,995 / human 2,961

### Non-canonical RBP Data
- **RBPWorld** (Li et al., 2025): Canonical/Non-canonical classification for 5 species  
  Files: `*RBPs.csv` format (*X. tropicalis* not available)

### Control Protein Data
- **AnimalTFDB 4.0** (Shen et al., 2023): Transcription factor annotations with tf_family
- **UniProt** KW-0418 (Kinase) and KW-0297 (GPCR): Retrieved via REST API

### Family Classification
- **Pfam domain-based** (unified across all protein types)
- Gene symbols → UniProt accessions (REST API) → Pfam domains → EuRBPDB 686-family reference
- Multiple Pfam hits: prioritize family with highest gene count in EuRBPDB
- No Pfam hit → Non-canonical
- Robustness confirmed at clan level (274–325 clans): identical ρ = 0.886

### IDR Analysis
- **MobiDB** (Piovesan et al., 2023): Consensus disorder predictions via REST API
- 12,332 of 13,138 RBPs (93.9% coverage)
- Metrics: IDR content (%), longest IDR length (aa), total disordered residues (aa)

### LLPS Analysis
- **LLPhyScore** (Cai et al., 2022): Computed from full-length UniProt sequences
- Within-proteome percentile ranks (species-normalized; 13,071 RBP genes)
- **PhaSePred** (Wang et al., 2022): Cross-validation for 5 species
  - catGRANULE, PScore, PLAAC
  - Note: zebrafish coverage only 18.8% (SwissProt restriction)

### Transcript Length Data
- **Ensembl BioMart** (release 115): 5'UTR, CDS, 3'UTR lengths for all protein-coding genes

### Phylogenetic Analysis
- **TimeTree** database: Time-calibrated tree for PGLS
- **R packages**: ape, nlme

---

## Repository Structure

```
rbp-evolution-analysis/
├── data/
│   ├── eurbpdb/              # EuRBPDB species files (*_RBP.txt)
│   ├── rbpworld/             # RBPWorld files (*RBPs.csv)
│   ├── processed/
│   │   ├── new_rbp_db.csv           # Unified RBP database (16,822 entries)
│   │   ├── new_tf_db.csv            # TF database (10,961 entries)
│   │   ├── eurbpdb_pfam_mapping.tsv # Pfam assignments for all RBPs
│   │   ├── controls_pfam_mapping.tsv# Pfam assignments for Kinase/GPCR
│   │   ├── idr_mobidb_results.tsv   # MobiDB IDR predictions
│   │   ├── llphyscore_per_gene.tsv  # LLPhyScore (gene-level)
│   │   └── summary_statistics_main6.csv
│   └── biomart/              # Ensembl BioMart UTR/CDS data
├── scripts/
│   ├── 01_pfam_classification.py    # Pfam domain assignment pipeline
│   ├── 02_correlations.py           # Spearman correlations + bootstrap + LOO
│   ├── 03_idr_analysis.py           # MobiDB API retrieval and analysis
│   ├── 04_llps_analysis.py          # LLPhyScore within-proteome ranking
│   ├── 05_pgls_analysis.R           # PGLS with ape/nlme
│   ├── 06_figures.py                # All main and supplementary figures
│   └── utils/
│       ├── uniprot_api.py           # UniProt REST API helpers
│       └── mobidb_api.py            # MobiDB REST API helpers
├── figures/                  # Generated figures (PDF/PNG)
├── results/
│   ├── main6_all_correlations_final.csv
│   ├── idr_correlations.csv
│   └── pgls_results.csv
└── README.md
```

---

## Reproducing the Analysis

### Requirements
```
Python 3.10+
pandas >= 1.5
scipy >= 1.9
numpy >= 1.23
matplotlib >= 3.6
requests >= 2.28

R >= 4.2
ape >= 5.6
nlme >= 3.1
```

### Step-by-step

**1. RBP family classification (Pfam-based)**
```bash
python scripts/01_pfam_classification.py \
  --eurbpdb_dir data/eurbpdb/ \
  --output data/processed/eurbpdb_pfam_mapping.tsv
```

**2. Correlation analysis**
```bash
python scripts/02_correlations.py \
  --pfam_mapping data/processed/eurbpdb_pfam_mapping.tsv \
  --output results/main6_all_correlations_final.csv
```

**3. IDR analysis (MobiDB API)**
```bash
python scripts/03_idr_analysis.py \
  --rbp_db data/processed/new_rbp_db.csv \
  --output data/processed/idr_mobidb_results.tsv
# Note: Requires internet access; ~2-3 hours for full retrieval
```

**4. LLPS analysis (LLPhyScore)**
```bash
# Requires LLPhyScore installation (fukunagatsu/LLPhyScore on GitHub)
# Run per species FASTA, then aggregate:
python scripts/04_llps_analysis.py \
  --llphyscore_dir /path/to/llphyscore_outputs/ \
  --rbp_db data/processed/new_rbp_db.csv \
  --output data/processed/llphyscore_per_gene.tsv
```

**5. PGLS analysis (R)**
```bash
Rscript scripts/05_pgls_analysis.R \
  --input results/main6_all_correlations_final.csv \
  --tree data/timetree_6species.nwk \
  --output results/pgls_results.csv
```

**6. Generate figures**
```bash
python scripts/06_figures.py
# Outputs to figures/
```

---

## Key Files for Reproducing Main Results

| File | Content | Used for |
|---|---|---|
| `data/processed/new_rbp_db.csv` | 16,822 RBP entries (EuRBPDB unified) | All primary analyses |
| `data/processed/eurbpdb_pfam_mapping.tsv` | Pfam domain assignments per RBP | Family diversity calculation |
| `results/main6_all_correlations_final.csv` | All Spearman correlations | Main Table / Figure 1 |
| `data/processed/idr_mobidb_results.tsv` | MobiDB IDR predictions | Figure 4 |
| `data/processed/llphyscore_per_gene.tsv` | LLPhyScore per gene | Figure 2 |

---

## Citation

Yasuda K. (2026). RNA-Binding Protein Family Diversification Correlates with Neural Complexity Across Metazoan Evolution. *iScience*.

---

## Contact

Kyota Yasuda  
Graduate School of Integrated Sciences for Life, Hiroshima University  
kyotay12@hiroshima-u.ac.jp
