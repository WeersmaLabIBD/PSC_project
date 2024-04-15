######################################################
# Plot heritability analysis results (specificity gene)
# Date: 06-04-2024
######################################################

library(ggplot2)
library(dplyr)

# Read prioritization output file
PSC_prioritization <- read.csv("/Users/s.qs/Documents/Chapters/4 PSC_Werna/20240405_heritability analysis_specificity_merged plasma/3_Plot heritability results/Input/prioritization.csv")

# Read cell type mapping
cell_type_mapping <-  read.csv("/Users/s.qs/Documents/Chapters/4 PSC_Werna/20240405_heritability analysis_specificity_merged plasma/1_CELLEX analysis/3_Conversion of Gene Name to Ensembl ID/Output/column_mapping.csv")

# Merge two files
PSC_prioritization_full <- merge(PSC_prioritization, cell_type_mapping, by.x = "annotation", by.y = "NewColumnName")

# Change column name
colnames(PSC_prioritization_full)[which(names(PSC_prioritization_full) == "OriginalColumnName")] <- "cell_type"

# benjamini-hochberg correction
PSC_prioritization_full$BH <- p.adjust(PSC_prioritization_full$pvalue, method = 'BH')

# Calculate -log10(p-value) and add to DataFrame
PSC_prioritization_full$neg_log10_BH <- -log10(PSC_prioritization_full$BH)

# Create a new column that contains an asterisk (*) if p-value < 0.05
PSC_prioritization_full$star <- ifelse(PSC_prioritization_full$BH < 0.05, "*", "")

# Set cell type order
unique(PSC_prioritization_full$cell_type)
cell_type_order <- c("Stem", "TA 1", "TA 2", "Cycling TA", "Immature Enterocytes 1", 
                     "Immature Enterocytes 2", "Enterocyte Progenitors", "Enterocytes", 
                     "M cells", "Best4+ Enterocytes",  "DUOX2 enterocytes", "Secretory TA",
                     "Immature Goblet", "Goblet", "Tuft", "Enteroendocrine", "WNT2B+ Fos-hi",
                     "WNT2B+ Fos-lo 1", "WNT2B+ Fos-lo 2", "RSPO3+", "WNT5B+ 1", "WNT5B+ 2", 
                     "Inflammatory Fibroblasts", "Myofibroblasts", "Endothelial", "Microvascular",
                     "Post-capillary Venules", "Pericytes", "Glia", "Plasma",
                     "Follicular", "GC", "Cycling B", "Macrophages", "DC1", "DC2", "CD69+ Mast",
                     "CD69- Mast", "Cycling Monocytes", "Inflammatory Monocytes", "Cycling T", 
                     "NKs", "ILCs", "CD8+ IELs", "CD8+ LP", "CD8+ IL17+", "CD4+ Activated Fos-hi", 
                     "CD4+ Activated Fos-lo", "CD4+ Memory", "CD4+ PD1+", "Tregs")
PSC_prioritization_full$cell_type <- factor(PSC_prioritization_full$cell_type, levels = cell_type_order)


# Create a heatmap and add asterisks
p <- ggplot(PSC_prioritization_full, aes(x = cell_type, y = gwas)) +
  geom_tile(aes(fill = neg_log10_BH)) + 
  scale_fill_gradient(low = "#deebfa", high = "#7fb6f5") + 
  geom_text(aes(label = star), color = "white", size = 4) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  labs(fill = "-log10(P CELLECT-LDSC)", x = "Cell Type", y = "GWAS") 

# Set the file path and name for saving the graph
file_path <- "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20240405_heritability analysis_specificity_merged plasma/3_Plot heritability results/Output/benjamini-hochberg correction/"

# save the graph
ggsave("PSC_prioritization_heatmap.pdf", plot = p, width = 18, height = 2, path = file_path)
