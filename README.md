# Hedgehog pathway activity and ICB signature correlation in ovarian cancer (GSE130000)

Single-cell analysis of Hedgehog (HH) pathway gene expression and its correlation
with published immune-checkpoint-blockade (ICB) response signatures in the
GSE130000 ovarian cancer scRNA-seq dataset (Kan et al., *Oncogene* 2022),
re-processed via TISCH2.

This repository contains all R code, intermediate result tables, and figure
files supporting the manuscript (in preparation).

---

## Data availability

This repository does **not** redistribute the source expression matrix.
Download the harmonised processed data from TISCH2:

- **TISCH2 dataset page:** http://tisch.comp-genomics.org → search `OV_GSE130000`
- **Files required (place in the same directory as the R scripts, or one level
  above if using the layout below):**
  - `OV_GSE130000_expression.h5`  (≈ 135 MB; sparse log-normalised counts)
  - `OV_GSE130000_CellMetainfo_table.tsv` (≈ 1.1 MB; cell-type / patient / tissue annotations)

The original raw data are deposited at GEO under accession
[**GSE130000**](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130000)
(Kan T et al., *Oncogene* 2022; DOI:
[10.1038/s41388-021-02139-z](https://doi.org/10.1038/s41388-021-02139-z)).

---

## Repository contents

```
.
├── Hedgehog_ICB_analysis.R              # initial analysis (composite ICB) — Fig 1–8
├── Hedgehog_ICB_signatures_analysis.R   # 9 published ICB signatures — Fig 9–14
├── README.md                            # this file
├── LICENSE                              # MIT
├── CITATION.cff                         # citation metadata
├── sessionInfo.txt                      # R + package versions used
├── .gitignore                           # excludes data files and OS cruft
│
├── Fig1_*.{pdf,png}      ── Main Fig 1A/B: HH pathway by cell type & tissue (violin)
├── Fig1_alt_*.{pdf,png}  ── Main Fig 1C:   HH pathway dot plot (%active × mean score)
├── Fig14_*.{pdf,png}     ── Main Fig 1D:   Tumor vs Relapsed paired dot plot (Fisher z)
├── Fig10_*.{pdf,png}     ── Main Fig 1E:   HH pathway vs 9 ICB signatures (malignant)
├── Fig13_*.{pdf,png}     ── Main Fig 1F:   GLI2 vs 9 ICB signatures
│
├── Fig4_*.{pdf,png}      ── Supp Fig 1:    HH gene heatmap (tumour-only, all cell types)
├── Fig9_*.{pdf,png}      ── Supp Fig 2:    HH gene × ICB signature correlation heatmap
│
├── Fig2,3,5–8,11,12_*    ── exploratory / superseded panels (kept for transparency)
│
└── *.csv                                 # all numeric result tables (see "Result tables")
```

### Result tables (CSV)

| File | What it contains |
|---|---|
| `HH_pathway_scores_by_celltype_tissue.csv` | Median HH pathway score per cell type × tissue |
| `HH_gene_expression_all_celltypes.csv` | Mean/% expression per HH gene per cell type |
| `HH_gene_expression_malignant.csv` | Per-tissue HH gene expression in malignant cells |
| `HH_genes_vs_ICB_signatures_correlations.csv` | Spearman ρ per HH gene × 9 ICB signatures (malignant) |
| `HH_vs_ICB_signatures_correlations.csv` | Pathway-level HH × 9 ICB signatures (malignant) |
| `HH_ICB_correlation_by_celltype.csv` | HH × composite ICB Spearman ρ per cell type |
| `HH_ICB_correlation_malignant_by_tissue.csv` | HH × ICB ρ split by Tumor/Met/Relapse |
| `HH_ICB_patient_pseudobulk.csv` | Patient-level pseudobulk HH vs ICB |
| `HH_gene_vs_ICB_correlation.csv` | Per-HH-gene Spearman ρ vs composite ICB |
| `HH_vs_ICB_patient_pseudobulk.csv` | Patient-level pseudobulk pairing |

---

## Reproducing the analysis

### 1. System requirements

- R ≥ 4.2 (developed and verified on R 4.5.3; see `sessionInfo.txt`)
- macOS, Linux, or Windows
- ≈ 4 GB free RAM (sparse matrix is small)

### 2. Install R packages

```r
install.packages(c("hdf5r", "Matrix", "dplyr", "tidyr", "ggplot2", "patchwork"))
```

Tested versions (from `sessionInfo.txt`):

| Package    | Version |
|------------|---------|
| hdf5r      | 1.3.12  |
| Matrix     | 1.7-4   |
| dplyr      | 1.2.0   |
| tidyr      | 1.3.2   |
| ggplot2    | 4.0.2   |
| patchwork  | 1.3.2   |

### 3. Get the data

Download from TISCH2 (`OV_GSE130000`) and place both files
(`OV_GSE130000_expression.h5`, `OV_GSE130000_CellMetainfo_table.tsv`) into
**this folder** (the same directory as the R scripts), **or** into the parent
directory and uncomment `setwd("..")` near the top of each script.

### 4. Run

```r
# from R, with the working directory set so that the data files are visible:
source("Hedgehog_ICB_analysis.R")             # generates Fig 1–8 + CSVs
source("Hedgehog_ICB_signatures_analysis.R")  # generates Fig 9–14 + ICB-signature CSVs
```

Each script is self-contained and re-creates every figure / CSV it depends on.
Run time: ≈ 2–4 min on a modern laptop.

---

## Cohort caveat

Per the original publication's Supplementary Table 1, the GSE130000 cohort
comprises 8 samples of mixed histology:

| Sample | Tissue | Histology |
|--------|---------|---|
| P1 | Primary tumour | High-grade serous |
| P2 | Primary tumour | Low-grade serous |
| P3 | Primary tumour | Low-grade serous |
| P4 | Primary tumour | Endometrioid |
| M1 | Peritoneal metastasis | High-grade serous |
| M2 | Peritoneal metastasis | Low-grade serous |
| R1 | Relapse tumour | Serous (grade unspecified) |
| R2 | Relapse tumour | Low-grade serous |

**Patient age, germline BRCA1/2 status, FIGO stage, and prior treatment were
not reported in the source publication and are therefore not available for
stratified analysis.**

---

## Repository

GitHub: **https://github.com/mythreyelab/Hedgehog-ovarian-scRNAseq**

## Citation

If you use this code, please cite:

1. **This repository** — see `CITATION.cff` (GitHub renders a "Cite this" button).
2. **The source dataset** — Kan T, *et al.* Single-cell RNA-seq recognized the
   initiator of epithelial ovarian cancer recurrence. *Oncogene* 41:895–906 (2022).
   DOI: [10.1038/s41388-021-02139-z](https://doi.org/10.1038/s41388-021-02139-z).
3. **TISCH2** — Han Y, *et al.* TISCH2: expanded datasets and new tools for
   single-cell transcriptome analyses of the tumour microenvironment.
   *Nucleic Acids Res* 51 (D1):D1425–D1431 (2023). DOI:
   [10.1093/nar/gkac959](https://doi.org/10.1093/nar/gkac959).

---

## License

MIT — see `LICENSE`.

---

## Contact

Mythreye Karthikeyan — University of Alabama at Birmingham
