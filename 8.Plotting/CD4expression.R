##################################################
# Check CD4 expression in CD8+IL17 for fig 4c,d  #
##################################################

library(Seurat)
library(patchwork)
library(ggplot2)
library(dplyr) 
library(tidyr) 

imm <- readRDS("~/Documents/R/PSC/Nieuw/imm_azimuth_with_plasma.rds")
DefaultAssay(imm) = "RNA"
Idents(imm) <- "predicted.cell_type.pred"

# define the list of genes of interest
genes <- c("CD4")
list.celltypes <- c("CD8+ IL17+", "Cycling T", "Tregs", "CD4+ PD1+")

#subset to have only PSC-I condition or UC-I condition
psc_inflamed <- imm[ ,imm@meta.data$state == 'PSC-I' & imm@meta.data$predicted.cell_type.pred %in% list.celltypes]
uc_inflamed <- imm[ ,imm@meta.data$state == 'UC-I'& imm@meta.data$predicted.cell_type.pred %in% list.celltypes]

#define the levels for the celltypes
psc_inflamed@meta.data$predicted.cell_type.pred <- factor(psc_inflamed@meta.data$predicted.cell_type.pred, levels = list.celltypes)
uc_inflamed@meta.data$predicted.cell_type.pred <- factor(uc_inflamed@meta.data$predicted.cell_type.pred, levels = list.celltypes)

#make plots
VlnPlot(psc_inflamed, features = genes, cols = c('#B3CDE3', '#CCEBC5', '#FFFFCC', '#F1B6DA'), idents = list.celltypes,
        group.by = "predicted.cell_type.pred") +NoLegend()
ggsave("Results/Figures_new/PSC_CD4_v2.pdf", width = 4, height = 4)
VlnPlot(uc_inflamed, features = genes, cols = c('#B3CDE3', '#CCEBC5', '#FFFFCC', '#F1B6DA'), idents = list.celltypes,
        group.by = "predicted.cell_type.pred") +NoLegend()
ggsave("Results/Figures_new/UC_CD4_v2.pdf", width = 4, height = 4)

