# Masters_2025-26_BioinformaticAnalysis
Scripts for the bioinformatic analysis of placental cytokine expression and comparisons with cancer

Repository contents
01_1stv3rd_PLDEG_filtering.R -  Filtering of first- versus third-trimester placental DEG data.
02_PLDEG_CytokineOverlap_PathwayAnalysis - Overlap of Cytokine Registry with filtered placental DEG data and pathway analysis.
03_GdcTransferTool - instructions for using the GDC transfer tool to download TCGA files.
04_ChemokineCandidates_PL_TCGA - CXCL2, CXCL12, CCL5 and CCL25 expression (TPM) in 2x placental datasets and 5x TCGA tumour tissue types.
05_PLDeconvolution - Placental DEG dataset rearrangement and single cell signature matrix construction for deconvolution analysis.

Data sources
All datasets are publicly available and required no ethical approval for secondary analysis. 
Placental DESeq2 and TPM data - obtained from GEO, accession code GSE247382 
Cytokine list - obtained from ImmPort, Cytokine Registry
Cancer TPM data - obtained from The Cancer Genome Atlas (TCGA)
Placental single cell FPKM data - obtained from Okae et al., (2018) Supplementary Table 1. ‘Transcriptome Profiling of Primary and Cultured Trophoblast Cells’

Software and packages
Analyses were performed in R v4.3.2. The packages used are specified in the individual scripts.

External software platforms required 
Reactome
STRING 
Cytoscape
Cibersort

Reproducibility
The scripts document the analytical steps used to generate the results reported in my thesis. Access to the original datasets, and external software platforms, may be required to reproduce the analyses.
