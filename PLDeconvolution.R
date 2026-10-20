
##ARRANGING 1ST AND 3RD TRIMESTER PLACENTAL TISSUE EXPRESSION DATA FOR DECONVOLUTION 
TMP_data <- read.csv("~/Downloads/Masters_2025/Datasets/TPM_Pisarska-mRNA2019_124CVS-43Placenta.csv", header=T)
View(TMP_data)
nrow(TMP_data) #27567
#removing ERCC spike-in control identifiers
TMP_ensembl <- TMP_data[1:27475,]

#retaining gene symbols >> following ensembl ID and _ in first column
colnames(TMP_ensembl)[1] <- "ensembl_id"
library(dplyr)

TMP_gene_symbol <- TMP_ensembl %>%
  mutate(ensembl_id = sub("^[^_]+_", "", ensembl_id))
View(TMP_gene_symbol)

#CVS (1st trimester) samples
TPM_CVS <- TPM_gene_symbol %>%
  select(1, starts_with("CVS"))
View(TPM_CVS) #125
write.table(TPM_CVS, "~/Downloads/Masters_2025/Filtering_files/Deconvolution/CVS_TPM", sep = "\t", row.names = FALSE, col.names = TRUE)
#Note this downloads as a "file" and has to be resaves as a "text file"

#PL (3rd trimester) samples
TPM_PL <- TPM_gene_symbol %>%
  select(1, starts_with("PL"))
View(TPM_PL) #44
write.table(TPM_PL, "~/Downloads/Masters_2025/Filtering_files/Deconvolution/PL_TPM", sep = "\t", row.names = FALSE, col.names = TRUE)
#Same note as above



###OKAE SIGNATURE MATRIX
library(readxl)
Okae_FPKM <- read_excel("~/Downloads/Masters_2025/Datasets/Okae_FPKM.xlsx")

#averaging duplicates in Okae
Okae_subset <- Okae_FPKM[,c(1, 3:13)]
colnames(Okae_subset)
nrow(Okae_subset) #24103
#removing # from column names 
colnames(Okae_FPKM) <- gsub("#", "_", colnames(Okae_FPKM))

sum(duplicated(Okae_subset$`Gene symbol`)) #321

library(dplyr)
Okae_dedup <- Okae_subset %>%
  group_by(`Gene symbol`) %>%
  summarise(across(.cols = where(is.numeric), .fns = ~mean(.x, na.rm = TRUE)), .groups = "drop")
nrow(Okae_dedup) #23782 = 24103-321

#log2(FPKM+1) > unlog and subtract 1 
#all columns except the first containing gene symbols
Okae_dedup[ , -1] <- 2^(Okae_dedup[ , -1]) - 1

#TPMx=(FPKMx/∑(FPKMall))×10^6 for genex
Okae_TPM <- Okae_dedup %>%
  mutate(across(where(is.numeric), ~ .x / sum(.x) * 1e6))

#average replicates
Okae_averaged <- data.frame(
  Gene = Okae_TPM[, 1],  # Keep the first column as Gene names
  CT = rowMeans(Okae_TPM[, c("CT _1", "CT _2", "CT _3")], na.rm = TRUE),
  EVT = rowMeans(Okae_TPM[, c("EVT _1", "EVT _2", "EVT _3")], na.rm = TRUE),
  ST = rowMeans(Okae_TPM[, c("ST _1", "ST _2", "ST _3")], na.rm = TRUE),
  Stroma = rowMeans(Okae_TPM[, c("Stroma _1", "Stroma _2")], na.rm = TRUE)
)

# Remove rows where all values are 0 > these genes are not expressed in placental cell types
Okae_no0 <- Okae_averaged[rowSums(Okae_averaged != 0) > 0, ]

write.table(Okae_no0, "~/Downloads/Masters_2025/Filtering_files/Deconvolution/Okae_TPM", sep = "\t", row.names = FALSE, col.names = TRUE)
#Note this downloads as a "file" and has to be resaves as a "text file"






##PLOTTING DECONVOLUTION RESULTS
rm(list = ls())

getwd()
setwd("C:/Users/mckge834/R studio/Datasets/Cytokine_bioinformatics_files")

CVS_immune_absolute <- read.csv("Deconvolution/CIBERSORTx_Job10_Results_CVS_TPM.csv")
CVS_immune_absolute <- CVS_immune_absolute[,1:23]
View (CVS_immune_absolute)
CVS_Okae_absolute <- read.csv("Deconvolution/CIBERSORTx_Job15_Results_CVS_Okae.csv")
CVS_Okae_absolute <- CVS_Okae_absolute[,1:5]
View (CVS_Okae_absolute)

library("dplyr")
CVS_join <- left_join(x=CVS_immune_absolute, y=CVS_Okae_absolute, by= "Mixture")
CVS_join$Mixture <- "CVS"
View(CVS_join)

CVS_Okae_absolute <- CVS_Okae_absolute %>% mutate(CVS_Okae_absolute, Signature_matrix = "Okae")


#Read in PL immune and arrange
PL_immune_absolute <- read.csv("Deconvolution/CIBERSORTx_Job11_Results_PL_TPM.csv")
PL_immune_absolute <- PL_immune_absolute[,1:23]
PL_Okae_absolute <- read.csv("Deconvolution/CIBERSORTx_Job17_Results_PL_Okae.csv")
PL_Okae_absolute <- PL_Okae_absolute[,1:5]

PL_join <- left_join(x=PL_immune_absolute, y=PL_Okae_absolute, by= "Mixture")
PL_join$Mixture <- "PL"
View(PL_join)

#combine datasets
CVS_PL_combined <- bind_rows(CVS_join, PL_join)
dim(CVS_PL_combined) #167 27
View(CVS_PL_combined)

#Re-naming columns - removing . and replacing with spaces
colnames(CVS_PL_combined) <- gsub("\\.", " ", colnames(CVS_PL_combined))
View(CVS_PL_combined)

#Reshape both to long format
library(tidyr)
CVS_PL_long <- CVS_PL_combined %>%
  pivot_longer(cols = -c(Mixture), names_to = "CellType", values_to = "Absolute score")
View(CVS_PL_long)

#log transform
CVS_PL_long_log2 <- CVS_PL_long %>%
  mutate(`Absolute score` = log2(`Absolute score` + 1))
View(CVS_PL_long_log2)


library(tidyverse)

#Define the special cell types
Placenta <- c( "CT",  "ST", "EVT", "Stroma")


#Assign facet group based on cell type
CVS_PL_long <- CVS_PL_long %>%
  mutate(FacetGroup = if_else(CellType %in% Placenta, "Placenta", "Immune"),
         FacetGroup = factor(FacetGroup, levels = c("Immune", "Placenta")))
levels(CVS_PL_long_log2$FacetGroup)
View(CVS_PL_long_log2)


CVS_PL_long <- CVS_PL_long %>%
  mutate(New_Joined_Column = paste0(.[[1]], "_", .[[4]]))
View(CVS_PL_long_log2)


#Keep CellType order the same as in the original dataframe
CVS_PL_long$CellType <- factor(
  CVS_PL_long$CellType,
  levels = colnames(CVS_PL_combined)[colnames(CVS_PL_combined) != "Mixture"]
)


custom_colors <- c(
  "CVS_Immune" = "#93b1cf",
  "PL_Immune" = "#647588",
  "CVS_Placenta" = "#E8C6DB",
  "PL_Placenta" = "#C582B2")

# Plot
ggplot(CVS_PL_long, aes(x = CellType, y = `Absolute score`, fill = New_Joined_Column)) +
  geom_bar(
    stat = "summary",
    fun = "mean",
    position = position_dodge(0.8),
    colour = "darkgrey",
    linewidth = 0.3
  ) +
  geom_errorbar(
    stat = "summary",
    fun.data = mean_se,
    position = position_dodge(0.8),
    width = 0.3,
    size = 0.3
  ) +
  scale_fill_manual(
    values = custom_colors,
    name = "Group",
    labels = c(
      "CVS_Immune" = "1st trimester",
      "PL_Immune" = "3rd trimester",
      "CVS_Placenta" = "1st trimester",
      "PL_Placenta" = "3rd trimester"
    )
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgray"),
    axis.ticks = element_line(color = "darkgray"),
    legend.title = element_blank(),
    legend.position = "bottom",
    plot.margin = margin(20, 20, 20, 20)
  ) +
  facet_grid(. ~ FacetGroup, scales = "free_x", space = "free_x", switch = "x") +
  labs(x = "Cell type", y = "Absolute score")






