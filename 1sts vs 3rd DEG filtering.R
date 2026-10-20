##FILTERING THE GONZALEZ ET AL FIRST VERSUS THIRD TRIMESTER PLACENTAL DESEQ2 DATASET

getwd()

#loading 1st vs 3rd DEGs
DEG_data <- read.csv("GSE247382_DESeq2_First-versus-Third-Trimester_sexAdjusted_124CVS-43PL.csv", skip = 3, header = F)
#excluding 3 rows of headers
nrow(DEG_data) #25312 = 25315-3

#naming important columns
colnames(DEG_data)[2] <- "gene_symbol"
colnames(DEG_data)[9] <- "log2FC"
colnames(DEG_data)[14] <- "pFDR"

#creating a subset dataframe with useful info
DEG_subset <- DEG_data[, c(as.character("gene_symbol"), "log2FC", "pFDR")]
  nrow(DEG_subset) #25312
  
#removing NA values
sum(is.na(DEG_subset$pFDR)) #2624
DEG_subset_NAfilter <- DEG_subset[!is.na(DEG_subset$pFDR), ]
sum(is.na(DEG_subset_NAfilter$pFDR)) #0
nrow(DEG_subset_NAfilter)#22688 = 25312-2624
  
#averaging log2FC and pFDR across duplicated genes
sum(duplicated(DEG_subset_NAfilter$gene_symbol)) #74 duplicates
library(dplyr)
DEG_dedup <- DEG_subset_NAfilter %>%
  group_by(gene_symbol) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    pFDR = mean(pFDR, na.rm = TRUE),
    .groups = "drop"  # removes grouping
  )
nrow(DEG_dedup) #22614 = 22688-74

#pvalue cut off 
sum(DEG_dedup$pFDR > 0.05) #6046 #How many gene symbols have p>0.05 and need to be removed
DEG_p0.05 <- DEG_dedup[DEG_dedup$pFDR <= 0.05, ] #retaining only gene symbols with p<=0.05
nrow(DEG_p0.05) #16568 = 22614-6046
max(DEG_p0.05$pFDR) #0.04998037 - only p<0.05
min(DEG_p0.05$pFDR) #0

#log2FC 1/-1 threshold
sum(DEG_p0.05$log2FC > 1) #2,515
sum (DEG_p0.05$log2FC < -1) #4,703
# 2,515 + 4,703 = 7,218 gene symbols to retain
sum(DEG_p0.05$log2FC <= 1 & DEG_p0.05$log2FC >= -1) #9,350 gene symbols to remove
#7,218 + 9,350 = 16,568
# Create a new dataframe with genes where log2FC is >1 or <-1
DEG_log2FC1 <- DEG_p0.05[DEG_p0.05$log2FC > 1 | DEG_p0.05$log2FC < -1, ]
nrow(DEG_log2FC1) #7,218

write.csv (DEG_log2FC1, "~/Downloads/Masters_2025/Filtering_files/DEG_filtered_log2FC1(2).csv", row.names = F)



##VOLCANO PLOT 
#working from NA, duplicate average data - to display all genes on volcano

#adjusting 0 values before log
DEG_dedup$pFDR[DEG_dedup$pFDR == 0] <- 0.0001

#Significance column
DEG_dedup$significance <- "Not significant"
DEG_dedup$significance[DEG_dedup$pFDR < 0.05 & DEG_dedup$log2FC > 1] <- "Significant upregulation"
DEG_dedup$significance[DEG_dedup$pFDR < 0.05 & DEG_dedup$log2FC < -1] <- "Significant downregulation"

count_not_sig <- sum(DEG_dedup$significance == "Not significant")
count_up <- sum(DEG_dedup$significance == "Significant upregulation")
count_down <- sum(DEG_dedup$significance == "Significant downregulation")


#Plot
ggplot() +
  geom_point(data = DEG_dedup, aes(x = log2FC, y = -log10(pFDR), color = significance), 
             alpha = 0.6, size = 1.5) +
  scale_color_manual(
    values = c(
      "Significant upregulation" = "#6B1F3A", 
      "Significant downregulation" = "#334264", 
      "Not significant" = "grey"
    ),
    labels = c(
      "Significant upregulation" = "Significant\nupregulation", 
      "Significant downregulation" = "Significant\ndownregulation", 
      "Not significant" = "Not\nsignificant"
    )
  ) +
  theme_minimal() +
  labs(x = "Log2 fold change", y = "-Log10 FDR adjusted p value") +
  theme(
    axis.text = element_text(size = 13),
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    legend.position = "bottom",  # Position the legend at the bottom
    legend.title = element_blank(), # Remove the legend title (the "Significance" label)
    legend.text = element_text(size = 13),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "darkgray"),
    axis.ticks = element_line(color = "darkgray")
  ) +
  geom_hline(yintercept = 1.3010, linetype = "dashed", color = "black") +  # p-value threshold of 0.05
  geom_vline(xintercept = c(1, -1), linetype = "dashed", color = "black") +  # log2FC thresholds
  annotate("text", x = -7.5, y = 300, 
           label = count_down, 
           hjust = 0.5, color = "#334264", size = 4.5) +
  annotate("text", x = 0, y = 300, 
           label = count_not_sig, 
           hjust = 0.5, color = "grey30", size = 4.5) +
  annotate("text", x = 6, y = 300, 
           label = count_up, 
           hjust = 0.5, color = "#6B1F3A", size = 4.5)
