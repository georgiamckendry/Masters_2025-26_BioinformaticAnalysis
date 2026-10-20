rm(list = ls())
getwd()
setwd("C:/Users/mckge834/R studio/Datasets/Cytokine_bioinformatics_files")

#INNER JOIN OF filtered lst vs 3rd DEGS WITH CYTOKINE REGISTRY
#reading in filtered DEGs

DEG_filtered <- read.csv("DEG_filtered_log2FC1(2).csv")
colnames(DEG_filtered) #"gene_symbol"
nrow(DEG_filtered) #7218

#reading in Cytokine Registry
library(readxl)
Cytokine_registry <- read_excel ("CytokineRegistry.November_2015__gene symbols edited.xls")
colnames(Cytokine_registry)  #EntrezGene Symbol (Human)
#ommiting column headers, 1x duplicate
library(dplyr)
Cytokine_registry_unique <- Cytokine_registry$`EntrezGene Symbol (Human)` %>%
  na.omit() %>% #267  
  unique() #as is only a list of gene names, no averaging across duplicates is required
length(Cytokine_registry_unique) #266

#vector of cytokines from the registry > name for join
Cytokine_registry_unique <- data.frame(
  gene_symbol = Cytokine_registry_unique
)

#inner join
DEG_cytokines_innerjoin <- inner_join(x=DEG_filtered, y=Cytokine_registry_unique, by = "gene_symbol")
nrow(DEG_cytokines_innerjoin)#101 

write.csv(DEG_cytokines_innerjoin, 
          "~/Downloads/Masters_2025/Filtering_files/DEG_filtered_cytokines_innerjoin.csv", 
          row.names = FALSE)


##RUN INNER JOIN THROUGH REACTOME WEB BROWSER
#...

#VENNDIAGRAM PLOT OF DEG-CYTOKINE REGISTRY OVERLAP
install.packages("VennDiagram")
library(VennDiagram)
library(grid)
grid.newpage() +
draw.pairwise.venn(
  area1 = 7218,
  area2 = 266,
  cross.area = 101,
  category = c("Filtered DEGs", "Cytokine Registry"),
  fill = c("grey", "#93b1cf"),
  cat.cex = 0.8,
  cat.dist = c(0.1, 0.1),     
  fontfamily = "Arial",
  cat.fontfamily = "Arial",
  cat.fontface = "bold",
  margin = 0.1,
  print.mode = c("raw", "percent"),
  position.method = "centroid",
  alpha = 0.6
)


#VOLCANOPLOT HIGHLIGHTING CANDIDATE CYTOKINES
#loading 1st vs 3rd DEGs
DEG_data <- read.csv("GSE247382_DESeq2_First-versus-Third-Trimester_sexAdjusted_124CVS-43PL.csv", skip = 3, header = F)
#excluding 3 rows of headers
nrow(DEG_data) #25312 (25315-3)
View(DEG_data)

#naming important columns
colnames(DEG_data)[2] <- "gene_symbol"
colnames(DEG_data)[9] <- "log2FC"
colnames(DEG_data)[14] <- "padj_FDR"

#removing NA values
sum(is.na(DEG_data$padj_FDR)) #2624
DEG_data_NAfilter <- DEG_data[!is.na(DEG_data$padj_FDR), ]
sum(is.na(DEG_data_NAfilter$padj_FDR)) #0
nrow(DEG_data_NAfilter)#22688 = 25312-2624

#creating a subset dataframe with useful info
DEG_subset <- DEG_data_NAfilter[, c(as.character("gene_symbol"), "log2FC", "padj_FDR")]
nrow(DEG_subset) #22688

#averaging log2FC and padj_FDR across duplicated genes
sum(duplicated(DEG_subset$gene_symbol)) #74 duplicates

DEG_dedup <- DEG_subset %>%
  group_by(gene_symbol) %>%
  summarise(
    log2FC = mean(log2FC, na.rm = TRUE),
    padj_FDR = mean(padj_FDR, na.rm = TRUE),
    .groups = "drop"  # removes grouping
  )
nrow(DEG_dedup) #22614 = 22688-74

#removing 0s for log transformation
DEG_dedup$padj_FDR[DEG_dedup$padj_FDR == 0] <- 0.0001

# Create a significance column
DEG_dedup$significance <- "Not significant"
DEG_dedup$significance[DEG_dedup$padj_FDR < 0.05 & DEG_dedup$log2FC > 1] <- "Significant upregulation"
DEG_dedup$significance[DEG_dedup$padj_FDR < 0.05 & DEG_dedup$log2FC < -1] <- "Significant downregulation"
# Overwrite label if gene is cytokine candidate
DEG_dedup$significance[DEG_dedup$gene_symbol %in% DEG_cytokines_innerjoin$gene_symbol] <- "Cytokine"

count_not_sig <- sum(DEG_dedup$significance == "Not significant")
count_up <- sum(DEG_dedup$significance == "Significant upregulation")
count_down <- sum(DEG_dedup$significance == "Significant downregulation")
count_cytokine <- sum(DEG_dedup$significance == "Cytokine")
count_cytokine #101

# Split the data into two subsets
DEG_dedup_cytokine <- DEG_dedup[DEG_dedup$significance == "Cytokine", ]
DEG_dedup_other <- DEG_dedup[DEG_dedup$significance != "Cytokine", ]

# Plot
ggplot() +
  # Plot non-cytokine genes first
  geom_point(data = DEG_dedup_other, aes(x = log2FC, y = -log10(padj_FDR), color = significance), 
             alpha = 0.6, size = 2) +
  # Plot cytokine genes on top
  geom_point(data = DEG_dedup_cytokine, aes(x = log2FC, y = -log10(padj_FDR), color = significance), 
             alpha = 0.6, size = 2) +
  scale_color_manual(
    values = c(
      "Significant upregulation" = "#6B1F3A", 
      "Significant downregulation" = "#334264", 
      "Not significant" = "grey", 
      "Cytokine" = "#93b1cf"
    ),
    labels = c(
      "Significant upregulation" = "Significant\nup-regulation", 
      "Significant downregulation" = "Significant\ndown-regulation", 
      "Not significant" = "Not\nsignificant",
      "Cytokine" = "Cytokine\n candidates"
    )
  ) +
  theme_minimal() +
  labs(x = "Log2 fold change", y = "-Log10 adjusted p-value") +
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
           hjust = 0.5, color = "#334264", size = 5) +
  annotate("text", x = 0, y = 300, 
           label = count_not_sig, 
           hjust = 0.5, color = "grey30", size = 5) +
  annotate("text", x = 6, y = 300, 
           label = count_up, 
           hjust = 0.5, color = "#6B1F3A", size = 5) +
  annotate("text", x = -7.5, y = 275, 
           label = count_cytokine, 
           hjust = 0.5, color = "#93b1cf", size = 5)


#REACTOME PATHWAY ANALYSIS BAR GRAPH
#Read in Pathway analysis results csv 
DEG_cytokines_reactome <- read.csv("DEG_cytokines_reactome.csv", header = T)

#Create a new column and add counts of candidates in the Submitted.entites column - number associated with each pathway
DEG_cytokines_reactome$Mapped_entities <- sapply(DEG_cytokines_reactome$Submitted.entities.found, function(x) length(strsplit(x, ";")[[1]]))

write.csv(DEG_cytokines_reactome, "~/Downloads/Masters_2025/Filtering_files/DEG_cytokines_reactome.csv")

#Data frame of top 20 pathways - ranked by p value
DEG_cytokines_reactome_10 <- DEG_cytokines_reactome[1:10,]

#descending order for plotting
DEG_cytokines_reactome_10$Pathway.name <- factor(
  DEG_cytokines_reactome_10$Pathway.name,
  levels = rev(DEG_cytokines_reactome_10$Pathway.name)
)

#horizontal bar pot
library(ggplot2)
library(stringr)
ggplot(DEG_cytokines_reactome_10, aes(x = Mapped.entities, y = Pathway.name)) +
  geom_bar(stat = "identity", aes(fill = Entities.pValue)) +
  scale_y_discrete(labels = function(x) str_wrap(x, width = 45)) +
  scale_fill_gradient(low = "#93b1cf", high = "#647588", 
                      name = "FDR\nadjusted\np value",
                      guide = guide_colorbar(reverse = TRUE))+
  labs(x = "Enriched gene count",
       y = "Reactome pathway") +
  theme(legend.position = "bottom") +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 10, lineheight = 0.5, margin = margin(r = 5)),  # Adjust axis text size
        axis.title = element_text(size = 11),  # Adjust axis title size
        plot.margin = margin(20, 70, 20, 20),  # Add space around the plot
        legend.text = element_text(size = 9), 
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),  # Smaller keys
        legend.spacing.x = unit(0.5, "cm"),  # Space between legend items
        legend.margin = margin(t = 5, b = 5),  # Add a bit of top/bottom padding
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(color = "darkgray"),          # Adds both x and y axis lines
        axis.ticks = element_line(color = "darkgray"))





                