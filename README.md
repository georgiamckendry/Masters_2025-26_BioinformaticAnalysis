# Placental Cytokine and Cancer Bioinformatic Analysis

Scripts for the bioinformatic analysis of placental cytokine expression and comparisons with cancer.

## Repository Contents

- **`01_1stv3rd_PLDEG_filtering.R`** — Filtering of first- versus third-trimester placental differentially expressed gene (DEG) data.
- **`02_PLDEG_CytokineOverlap_PathwayAnalysis.R`** — Overlap of the ImmPort Cytokine Registry with filtered placental DEG data and pathway analysis.
- **`03_GdcTransferTool`** — Instructions for using the GDC Data Transfer Tool to download TCGA files.
- **`04_ChemokineCandidates_PL_TCGA.R`** — Analysis of CXCL2, CXCL12, CCL5 and CCL25 expression (TPM) across two placental datasets and five TCGA tumour tissue types.
- **`05_PLDeconvolution.R`** — Rearrangement of the placental DEG dataset and construction of a single-cell signature matrix for deconvolution analysis.

## Data Sources

All datasets are publicly available and were obtained for secondary analysis, which did not require ethical approval.

- **Placental DESeq2 and TPM data:** Gene Expression Omnibus (GEO), accession GSE247382.
- **Cytokine list:** ImmPort Cytokine Registry.
- **Cancer TPM data:** The Cancer Genome Atlas (TCGA).
- **Placental single-cell FPKM data:** Okae et al. (2018), Supplementary Table 1, *Transcriptome Profiling of Primary and Cultured Trophoblast Cells*.

## Software and Packages

Analyses were performed using R version 4.3.2. The R packages used are specified in the individual scripts.

### External Software Platforms

- Reactome
- STRING
- Cytoscape
- CIBERSORT

## Reproducibility

The scripts document the analytical steps used to generate the results reported in the thesis. Access to the original datasets and external software platforms may be required to reproduce the analyses.
