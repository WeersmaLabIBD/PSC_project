
###################################################
# Title: input for figure 3b
# Date: 28-08-2023
###################################################


######################
# Library
######################

library(ggplot2) 
library(dplyr) 
library(tidyr) 
library(stringr)
library(patchwork)
library(grid)
library(gridExtra) 
library(ggpubr)
library(RColorBrewer) 
library(Seurat)
library(readr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)
require(tidyverse)
library(packcircles)

######################
# Main codes
######################

# read in dataframes
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")

DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"

# Create a table with rows celltype and column disease/analysis, with number of DE genes per analysis for each celltype for UC and PSC
files <- list.files(path="Results/DE_2023/DE_PSCINI_all/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
total_PSC <- data.frame()
lapply(files, function(x) {
  t <- read.csv(x, header=TRUE) # load file
  celltype <- tools::file_path_sans_ext(basename(x))
  t$celltype <- celltype
  out <- t[t$p_val_adj < 0.05,]
    # write to file
  total_PSC <<- rbind(total_PSC, out)
})
total_PSC <- total_PSC[,c(6,8,9)]
DEpercelltype_PSC <- total_PSC %>%
  group_by(celltype) %>%
  summarize(num_genes = n())
#write_csv(DEpercelltype_PSC, "Results/DE_2023/DEpercelltype_PSC.csv")

files <- list.files(path="Results/DE_2023/DE_UCINI_all/", pattern="*.csv", full.names=TRUE, recursive=FALSE)
total_UC <- data.frame()
lapply(files, function(x) {
  t <- read.csv(x, header=TRUE) # load file
  celltype <- tools::file_path_sans_ext(basename(x))
  t$celltype <- celltype
  out <- t[t$p_val_adj < 0.05,]
  # write to file
  total_UC <<- rbind(total_UC, out)
})
total_UC <- total_UC[,c(6,8,9)]
DEpercelltype_UC <- total_UC %>%
  group_by(celltype) %>%
  summarize(num_genes = n())
#write_csv(DEpercelltype_UC, "Results/DE_2023/DEpercelltype_UC.csv")
