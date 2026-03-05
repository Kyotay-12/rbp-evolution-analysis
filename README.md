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
│   ├── fetch_pfam_from_uniprot.py   # Pfam domain retrieval via UniProt REST API
│   ├── fetch_controls_pfam.py       # Pfam retrieval for Kinase/GPCR controls
│   ├── fetch_idr_mobidb_v2.py       # MobiDB IDR predictions retrieval
│   ├── extract_llphyscore.py        # LLPhyScore extraction from output files
│   ├── aggregate_llphyscore.py      # LLPhyScore aggregation (gene-level)
│   ├── extract_phasepred_v2.py      # PhaSePred multi-method extraction
│   ├── analyze_correlations.py      # Main Spearman correlations (Figure 1)
│   ├── statistical_validation.py    # Bootstrap / LOO / Cohen's d
│   ├── idr_analysis.py              # IDR correlation analysis (Figure 4)
│   ├── supplementary_analyses.py    # Non-canonical / RBP-TF overlap / domain expansion
│   └── pgls_analysis.py             # PGLS with Pagel's lambda (Supp Figure 2)
├── results/
│   ├── main6_all_correlations_final.csv
│   ├── idr_correlations.csv
│   ├── idr_species_summary.csv
│   ├── 3utr_correlation_results.csv
│   ├── 3utr_rbp_summary_by_species.csv
│   └── utr_cds_correlations.csv
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

```

### Step-by-step

**1. Pfam domain retrieval (requires UniProt API access, ~hours)**
```bash
python scripts/fetch_pfam_from_uniprot.py \
  data/eurbpdb/ data/processed/eurbpdb_pfam_mapping.tsv

python scripts/fetch_controls_pfam.py \
  data/processed/ data/processed/controls_pfam_mapping.tsv
```

**2. Main correlation analysis (Figure 1)**
```bash
python scripts/analyze_correlations.py \
  data/processed/ results/
```

**3. Statistical validation (Bootstrap / LOO / Cohen s d / Pfam clan)**
```bash
python scripts/statistical_validation.py \
  data/processed/ results/
```

**4. IDR analysis (Figure 4)**
```bash
# Step 4a: Retrieve MobiDB predictions (~2-3 hours, requires internet)
python scripts/fetch_idr_mobidb_v2.py \
  data/processed/new_rbp_db.csv \
  data/processed/idr_mobidb_results_v2.tsv

# Step 4b: Compute species-level statistics and correlations
python scripts/idr_analysis.py \
  data/processed/ results/
```

**5. LLPS analysis (LLPhyScore)**
```bash
# Requires LLPhyScore installation: github.com/fukunagatsu/LLPhyScore
# Run LLPhyScore per species FASTA (see data/fasta/), then:
python scripts/extract_llphyscore.py
python scripts/aggregate_llphyscore.py
```

**6. Supplementary analyses (Non-canonical / RBP-TF overlap / domain expansion)**
```bash
python scripts/supplementary_analyses.py \
  data/processed/ results/
```

**7. PGLS analysis (Supplementary Figure 2)**
```bash
python scripts/pgls_analysis.py \
  data/processed/ results/
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
