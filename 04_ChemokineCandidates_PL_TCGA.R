
##PLOTTING CHEMOKINE CANDIDATES - 2xPL datasets and TCGA
#Connect to high capacity storage -TCGA files have been downloaded here using BASH and gdc transfer tool


rm(list = ls())

#setting to general directory with all cancer folders
setwd("C:/Users/mckge834/R studio/Datasets/TCGA_data")

#takes a little while to read in
combined_BRCA_tumour <- read.csv("combined_TPM_BRCA_tumour.csv")
combined_BRCA_tumour<- combined_BRCA_tumour[, -1]
dim(combined_BRCA_tumour) #59428 413

combined_BRCA_normal <- read.csv("combined_TPM_BRCA_normal.csv")
combined_BRCA_normal<- combined_BRCA_normal[, -1]
dim(combined_BRCA_normal) #59428 114

combined_LUSC_tumour <- read.csv("combined_TPM_LUSC_tumour.csv")
combined_LUSC_tumour <- combined_LUSC_tumour[, -1]
dim(combined_LUSC_tumour) #59428 512

combined_LUSC_normal <- read.csv("combined_TPM_LUSC_normal.csv")
combined_LUSC_normal <- combined_LUSC_normal[, -1]
dim(combined_LUSC_normal) #59428 52

combined_STAD_tumour <- read.csv("combined_TPM_STAD_tumour.csv")
combined_STAD_tumour<- combined_STAD_tumour[, -1]
dim(combined_STAD_tumour) #59428 413

combined_STAD_normal <- read.csv("combined_TPM_STAD_normal.csv")
combined_STAD_normal<- combined_STAD_normal[, -1]
dim(combined_STAD_normal) #59428 37

combined_LUAD_tumour <- read.csv("combined_TPM_LUAD_tumour.csv")
combined_LUAD_tumour<- combined_LUAD_tumour[, -1]
dim(combined_LUAD_tumour) #59428 541

combined_LUAD_normal <- read.csv("combined_TPM_LUAD_normal.csv")
combined_LUAD_normal<- combined_LUAD_normal[, -1]
dim(combined_LUAD_normal) #59428 60

combined_COAD_tumour <- read.csv("combined_TPM_COAD_tumour.csv")
combined_COAD_tumour<- combined_COAD_tumour[, -1]
dim(combined_COAD_tumour) #59428 482

combined_COAD_normal <- read.csv("combined_TPM_COAD_normal.csv")
combined_COAD_normal<- combined_COAD_normal[, -1]
dim(combined_COAD_normal) #59428 42


## BRINGING IN PLACENTAL DATA
setwd("C:/Users/mckge834/R studio/Datasets/Okae_TPM")

# Read in 1st trimester TPM file and filter for candidates
CVS_TPM<- read.table("CVS_TPM.txt", header= T)
colnames(CVS_TPM)[1] <- "gene_name"
View(CVS_TPM)

# Repeat for 3rd trimester
PL_TPM <- read.table("PL_TPM.txt", header= T)
colnames(PL_TPM)[1] <- "gene_name"



#CANDIDATE EXTRACTION
# Load required libraries
library(tidyr)
library(dplyr)
library(ggplot2)


#BRINGING IN CHIS DATA
# CHI ALL HUB CHEMOKINES ACROSS 1ST AND 3RD TRIMESTER
library(dplyr)
library(readr)
library(purrr)
library(tidyr)

# Defining target genes
target_genes <- c("CXCL2", "CXCL12","CCL5", "CCL25")

## Defining first trimester samples
#sample A removed as a CXCL12 outlier
samples_1st <- c("C", "E", "G", "I", "J", "K", "L", "M", "O", "P", "Q", "S", "T")

#reading in 1st trimester files, creating table with candidates across all samples
#bind rows, for the function(s) which is referred to samples
Chi_1st <- map_dfr(samples_1st, function(s) {
  
  file <- paste0(
    "C:/Users/mckge834/R studio/Datasets/Chi_data/1st_trim/gene_abundance_P_",
    s,
    "_e.out"
  )
  
  read_tsv(file) %>%
    filter(`Gene Name` %in% target_genes) %>%
    select(`Gene Name`, TPM) %>%
    mutate(sample = s)
})


# Make wide format
Chi_1st_wide <- Chi_1st %>% pivot_wider(
  names_from = sample,
  values_from = TPM
)
View(Chi_1st_wide)
#rename gene column
library(dplyr)
Chi_1st_wide <- Chi_1st_wide %>%
  rename(gene_name = `Gene Name`)


## Defining third trimester samples
samples_3rd <- c("23", "24", "25", "26", "27", "28", "29", "30", "31", "32", "33", "34", "35", "36", "37", "38")

#reading in 3rd trimester files, creating table with candidates across all samples
#bind rows, for the function(s) which is referred to samples
Chi_3rd <- map_dfr(samples_3rd, function(s) {
  
  file <- paste0(
    "C:/Users/mckge834/R studio/Datasets/Chi_data/3rd_trim/gene_abundance_SRR31115",
    s,
    "_e.out"
  )
  
  read_tsv(file) %>%
    filter(`Gene Name` %in% target_genes) %>%
    select(`Gene Name`, TPM) %>%
    mutate(sample = s)
})


# Make wide format
Chi_3rd_wide <- Chi_3rd %>% pivot_wider(
  names_from = sample,
  values_from = TPM
)
View(Chi_3rd_wide)

Chi_3rd_wide <- Chi_3rd_wide %>%
  rename(gene_name = `Gene Name`)



# Define cancer datasets using abbreviations
cancer_datasets <- list(
  `BRCA tumour` = combined_BRCA_tumour,
  `BRCA normal` = combined_BRCA_normal,
  `LUAD tumour` = combined_LUAD_tumour,
  `LUAD normal` = combined_LUAD_normal,
  `LUSC tumour` = combined_LUSC_tumour,
  `LUSC normal` = combined_LUSC_normal,
  `STAD tumour` = combined_STAD_tumour,
  `STAD normal` = combined_STAD_normal,
  `COAD tumour` = combined_COAD_tumour,
  `COAD normal` = combined_COAD_normal,
  `Gonzalez first` = CVS_TPM,
  `Gonzalez third` = PL_TPM, 
  `Chi first` = Chi_1st_wide,
  `Chi third` = Chi_3rd_wide
)


# Initialize an empty list to store IL6 data
CXCL12_list <- list()
CXCL2_list <- list()
CCL5_list <- list()
CCL25_list <- list()


# Extract CXCL12 expression data from each dataset
library(dplyr)
library(tidyr)
for (cancer in names(cancer_datasets)) {
  df <- cancer_datasets[[cancer]]
  
  CXCL12_df <- df %>%
    filter(gene_name == "CXCL12") %>%
    select(-gene_name) %>%
    gather(key = "sample_id", value = "TPM") %>%
    mutate(cancer_type = cancer)
  
  CXCL12_list[[cancer]] <-  CXCL12_df
}

# Combine all CXCL12 data into a single data frame
CXCL12_all <- bind_rows(CXCL12_list)
View(IL1B_all) # run stat test on this file > for each cancer type, then could bind together in one file
#is already in long format


# Check means (optional summary)
CXCL12_all %>%
  group_by(cancer_type) %>%
  summarise(mean_TPM = mean(TPM, na.rm = TRUE))




##barplot

ggplot(CXCL12_all,
       aes(x = cancer_type,
           y = log2(TPM + 1),
           fill = cancer_type)) +
  
  stat_summary(
    fun = mean,
    geom = "bar",
    color = "darkgrey",
    linewidth = 0.7,
    width = 0.9
  ) +
  
  stat_summary(
    fun.data = mean_se,
    geom = "errorbar",
    width = 0.2,
    linewidth = 0.8,
    color = "black"
  ) +
  
  labs(y = "CXCL12  log2(TPM + 1)") +
  
  geom_hline(yintercept = 0,
             color = "darkgrey") +
  
  scale_fill_manual(values = c(
    "BRCA tumour" = "#51806a",
    "BRCA normal" = "#325756",
    "LUAD tumour" = "#51806a",
    "LUAD normal" = "#325756",
    "LUSC tumour" = "#51806a",
    "LUSC normal" = "#325756",
    "STAD tumour" = "#51806a",
    "STAD normal" = "#325756",
    "COAD tumour" = "#51806a",
    "COAD normal" = "#325756",
    "Gonzalez first" = "#E8C6DB",
    "Gonzalez third" = "#C582B2",
    "Chi first" = "#E8C6DB",
    "Chi third" = "#C582B2"
  )) +
  
  scale_x_discrete(limits = c(
    "Gonzalez first",
    "Gonzalez third",
    "Chi first",
    "Chi third",
    "STAD tumour",
    "STAD normal",
    "LUAD tumour",
    "LUAD normal",
    "LUSC tumour",
    "LUSC normal",
    "BRCA tumour",
    "BRCA normal",
    "COAD tumour",
    "COAD normal"
  )) +
  
  coord_cartesian(ylim = c(0, 10.5)) +
  
  scale_y_continuous(
    breaks = seq(0, 12.5, 2.5),
    expand = expansion(mult = c(0.05, 0.10))
  ) +
  
  theme_classic() +
  theme(
    axis.line = element_line(color = "darkgrey", linewidth = 1),
    axis.text.x = element_text(color = "black", size = 14),
    axis.text.y = element_text(color = "black", size = 12),
    axis.title.y = element_text(color = "black", size = 14),
    axis.ticks.x = element_line(color = "black"),
    legend.position = "none"
  )





##STAT TESTS
library(dplyr)
library(tidyr)

CXCL12_break <- CXCL12_all %>%
  separate(
    col = cancer_type,          # the column with "BRCA tumour"
    into = c("cancer", "tissue_type"), # new columns
    sep = " ",                  # split at the space
    remove = FALSE              # keep original column if you want
  )
View(CXCL12_break)

library(rstatix)
CXCL12_between_tissue <- CXCL12_break %>%
  group_by(cancer) %>%                     # one t-test per cancer
  t_test(TPM ~ tissue_type) %>%  # tumour vs normal
  adjust_pvalue(method = "BH") %>%
  add_significance()
View(CXCL12_between_tissue)
