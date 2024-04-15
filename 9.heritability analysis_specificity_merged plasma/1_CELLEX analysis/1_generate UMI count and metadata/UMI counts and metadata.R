########################################################################
# Generate UMI counts and metadata from PSC-UC to CELLEX analysis
# Date: 05-04-2024
# Author: Shiqiang Sun
########################################################################

#############################################
# libraries                                 #
#############################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(RColorBrewer)
library(pheatmap)
library(reshape2)

#############################################
# main codes                                #
#############################################

# Load the PSC bile duct Seurat dataset
PSC_data <- readRDS('/groups/umcg-weersma/tmp01/projects/PSC/ongoing/Seurat_object/PSC_colitis_allcomps_mergedcelltypes.rds')

# Merge IgA, IgG, IgM, and Ig_negative to Plasma
PSC_data$mergedcelltype2 <- PSC_data$mergedcelltype
PSC_data$mergedcelltype2[PSC_data$mergedcelltype2 %in% c("IgA", "IgG", "IgM", "Ig_negative")] <- "Plasma"

# Get UMI count from Seurat object
umi_counts <- GetAssayData(PSC_data, assay = "RNA", slot = "counts")
umi_counts_df <- as.data.frame(umi_counts)

# Save UMI count
write.csv(umi_counts_df, "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/input/PSC_UC_UMI_count_merged_Plasma.csv", row.names = TRUE)

# Get meta data from Seurat object
metadata <- as.data.frame(PSC_data@meta.data)

# Clean metadata, only keep 'mergedcelltype2' column
metadata <- metadata[, "mergedcelltype2", drop = FALSE]

# Save UMI count
write.csv(metadata, "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/input/PSC_UC_metadata_merged_Plasma.csv")

