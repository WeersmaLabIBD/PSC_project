#####################################################
# Generate DE genes and GO terms list per cell type #
# Input for basic_DE_2023_v2 and supp 3,4,5,6       #
#####################################################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(MAST)
library(readxl)
library(openxlsx)
library(enrichR)

#load in datasets
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")

#set idents to final celltyping and assay to RNA
DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"

# create full DE genes lists (no p-value cutoff) for PSC-I v PSC-NI and UC-I v UC-NI
  # PSC
data = str # Repeat for all three datasets (epi, stromal, immune)
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  path <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_PSCINI_all/",x,".csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "PSC-NI")
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"PSC_I","PSC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # save dataframes
  write.csv(x, path)
}
  
  # UC
data = str # Repeat for all three datasets (epi, stromal, immune)
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  path <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_UCINI_all/",x,".csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "UC-NI")
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"UC_I","UC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # save dataframes
  write.csv(x, path)
}
      ## Repeat for all three datasets (epi, stromal, immune)

#create GO terms and significant DE genes lists analysis PSC-I vs PSC-NI
data = epi
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  pathPSCI <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/",x,"_PSC_I_up.csv",sep="")
  pathPSCNI <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/",x,"_PSC_NI_up.csv",sep="")
  pathPSCI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/",x,"_PSC_I_up_GO.csv",sep="")
  pathPSCNI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/",x,"_PSC_NI_up_GO.csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "PSC-I", ident.2 = "PSC-NI")
  x <- filter(x, x$p_val_adj < 0.05)
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"PSC_I","PSC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # separate PSC_I up gene list and PSC_NI up gene list
  PSC_I_up_x = filter(x, x$state == "PSC_I")  
  PSC_NI_up_x = filter(x, x$state == "PSC_NI")
  # GO terms for each gene list 
  GO_PSCIup_x <- enrichr(PSC_I_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  GO_PSCNIup_x <- enrichr(PSC_NI_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  # save dataframes
  write.csv(PSC_I_up_x, pathPSCI)
  write.csv(PSC_NI_up_x, pathPSCNI)
  write.csv(GO_PSCIup_x, pathPSCI_GO)
  write.csv(GO_PSCNIup_x, pathPSCNI_GO)
}

    ## Repeat for all three datasets (epi, stromal, immune)

#DE analysis UC-I vs UC-NI
data = epi
celltypes <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # safe pathnames
  pathUCI <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/",x,"_UC_I_up.csv",sep="")
  pathUCNI <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/",x,"_UC_NI_up.csv",sep="")
  pathUCI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/",x,"_UC_I_up_GO.csv",sep="")
  pathUCNI_GO <- paste("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/",x,"_UC_NI_up_GO.csv",sep="")
  # find markers
  x <- FindMarkers(data, subset.ident = x, group.by = "state", test.use = "MAST", ident.1 = "UC-I", ident.2 = "UC-NI")
  x <- filter(x, x$p_val_adj < 0.05)
  if(nrow(x) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'state' column
  x$state = NA
  x = as.data.frame(x)
  x$state = ifelse(x$avg_log2FC > 0,"UC_I","UC_NI")
  # create a 'gene' column
  x$gene = NA
  x$gene = rownames(x)
  # separate UC_I up gene list and UC_NI up gene list
  UC_I_up_x = filter(x, x$state == "UC_I")  
  UC_NI_up_x = filter(x, x$state == "UC_NI")
  # GO terms for each gene list 
  GO_UCIup_x <- enrichr(UC_I_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  GO_UCNIup_x <- enrichr(UC_NI_up_x$gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
  # save dataframes
  write.csv(UC_I_up_x, pathUCI)
  write.csv(UC_NI_up_x, pathUCNI)
  write.csv(GO_UCIup_x, pathUCI_GO)
  write.csv(GO_UCNIup_x, pathUCNI_GO)
}

      ## Repeat for all three datasets (epi, stromal, immune)