# Epigenetic Modifications and Gene Expression in Colorectal Cancer
### Integrative Analysis of GEO Datasets | R | 2024

---

## Overview

In 2024, I was bored. I was also interested to learning some bioinformatics and epigenetics. This project investigates the relationship between epigenetic modifications and gene expression changes in colorectal cancer (CRC), using publicly available GEO microarray datasets. Three datasets were analysed to identify differentially expressed genes (DEGs) associated with chromosomal instability, early-onset CRC, and epigenetic regulation.

The goal was to identify candidate biomarkers relevant to CRC diagnosis, prognosis, and targeted therapy.

---

## Datasets

| GEO Accession | Description | Samples | Key Comparison |
|---|---|---|---|
| **GSE30540** | CRC patient tumour vs normal mucosa | 25 CIN-high / 10 CIN-low | CIN-high vs CIN-low phenotype |
| **GSE4107** | Early-onset CRC vs healthy controls | Mucosa of CRC patients vs healthy | CRC mucosa vs normal |
| **GSE32323** | 5-aza-2'-deoxycytidine treatment in CRC cell lines | Treated vs untreated | Methylation-regulated gene expression |

All datasets were retrieved programmatically from NCBI GEO using the `GEOquery` Bioconductor package. **No raw data files are stored in this repository** — all downloads are handled by the analysis scripts.

---

## Repository Structure

```
CRC-GEO-Analysis/
│
├── R/
│   ├── 00_setup.R             # Package installation and loading
│   ├── 01_GSE30540_CIN.R      # CIN-high vs CIN-low DEA (GSE30540)
│   ├── 02_GSE4107_earlyonset.R  # Early-onset CRC DEA (GSE4107)
│   ├── 03_GSE32323_methylation.R # 5-aza treatment DEA (GSE32323)
│   └── 04_cross_dataset_summary.R # Shared DEG summary across datasets
│
├── results/                   # Output tables (auto-generated, gitignored)
├── figures/                   # Output plots (auto-generated, gitignored)
│
└── README.md
```

---

## Methods

### Data Retrieval & Preprocessing
- Datasets downloaded via `GEOquery::getGEO()`
- Expression matrices extracted from `ExpressionSet` objects
- Phenotype metadata parsed from GEO series matrices
- Low-variance probes filtered prior to DEA
- Data normalised (quantile normalisation where applicable)

### Differential Expression Analysis
- **limma** (Linear Models for Microarray Data) used for all DEA
- Design matrices constructed from phenotype metadata
- Empirical Bayes moderation (`eBayes`) applied
- DEGs defined as: **|log2FC| > 1, adjusted p-value (BH) < 0.05**

### Visualisation
- Volcano plots per dataset (ggplot2)
- Heatmaps of top DEGs (pheatmap)
- PCA of samples coloured by group
- Gene expression distribution (bar/pie charts)

---

## Key Findings

- **GSE30540**: 104 DEGs between CIN-high and CIN-low (34 up, 70 down); nominal p<0.05, |log2FC|>1. Exploratory:underpowered for FDR correction (n=10 CIN-low).
- **GSE4107**: 2,421 DEGs in CRC mucosa vs healthy controls (1,866 up, 555 down); FDR-corrected (BH adj.p<0.05, |log2FC|>1).
- **GSE32323**: 3,433 DEGs in matched tumour vs normal tissue (1,755 up, 1,678 down); paired design blocking by patient, FDR-corrected. 60 epigenetically silenced gene candidates identified from 5-aza-dC cell line sub-experiment.

### Cross-Dataset Candidates
6 genes were significant across all three independent datasets:

| Gene | Full Name | Significance |
|---|---|---|
| **H19** | H19 imprinted maternally expressed transcript | Epigenetically regulated lncRNA, well-validated in CRC |
| **PIGR** | Polymeric immunoglobulin receptor | Mucosal immunity; consistently lost in CRC tumours |
| **GNG4** | G protein subunit gamma 4 | G-protein signalling, altered in multiple cancers |
| **CYP2B6** | Cytochrome P450 family 2 subfamily B member 6 | Drug metabolism; tumour microenvironment relevance |
| **PRAP1** | Proline rich acidic protein 1 | Expressed in normal colon epithelium, lost in CRC |
| **MAB21L2** | Mab-21 like 2 | Developmental transcription factor, emerging CRC role |

---

## How to Run

1. Clone this repository
2. Open `R/00_setup.R` and run it to install all dependencies
3. Run each numbered script in order, or source them individually

```r
source("R/00_setup.R")
source("R/01_GSE30540_CIN.R")
source("R/02_GSE4107_earlyonset.R")
source("R/03_GSE32323_methylation.R")
source("R/04_cross_dataset_summary.R")
```

All figures are saved to `figures/` and result tables to `results/`.

---

## Dependencies

- R >= 4.2.0
- Bioconductor packages: `GEOquery`, `limma`, `Biobase`
- CRAN packages: `ggplot2`, `pheatmap`, `dplyr`, `tidyr`, `ggrepel`

See `R/00_setup.R` for automated installation.

---

## Limitations & Future Work

- Sample sizes are small (n = 10–35 per group); findings should be treated as exploratory
- Microarray platforms differ across datasets; cross-dataset comparisons are indicative only
- **Planned v2**: Integrative re-analysis using TCGA-COAD (matched RNA-seq + DNA methylation 450k array, n ≈ 460) via Google BigQuery / ISB-CGC, with machine learning classification of CRC subtypes

---

## References

- GSE30540: Habermann et al. (2012)
- GSE4107: Sabates-Bellver et al. (2007)
- GSE32323: Khamas et al. (2012)
- Ritchie et al. (2015) limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Research*
- Davis & Meltzer (2007) GEOquery: a bridge between the Gene Expression Omnibus (GEO) and BioConductor. *Bioinformatics*
