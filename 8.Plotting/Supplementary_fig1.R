#####################################################
# Script to combine supplementary tables into files #
#####################################################

library(openxlsx)
library(dplyr)
library(Seurat)
library(stringr)
library(readr)

# Create supp table 1 with overview of celltypes and markers
  # Generate list with marker genes for each celltype
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

data <- epi
celltypes <- as.character(unique(Idents(data)))
celltypes_epi <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # find markers
  y <- FindMarkers(data, group.by = "celltype.final", test.use = "MAST", ident.1 = x, only.pos = T)
  y <- filter(y, y$p_val_adj < 0.05)
  if(nrow(y) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'gene' column
  y$gene = NA
  y$gene = rownames(y)
  # save dataframe
  assign( paste(x) , y) 
}

data <- imm
celltypes <- as.character(unique(Idents(data)))
celltypes_imm <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # find markers
  y <- FindMarkers(data, group.by = "predicted.cell_type.pred", test.use = "MAST", ident.1 = x, only.pos = T)
  y <- filter(y, y$p_val_adj < 0.05)
  if(nrow(y) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'gene' column
  y$gene = NA
  y$gene = rownames(y)
  # save dataframe
  assign( paste(x) , y)
}

data <- str
celltypes <- as.character(unique(Idents(data)))
celltypes_str <- as.character(unique(Idents(data)))
for (x in celltypes) {
  # find markers
  y <- FindMarkers(data, group.by = "predicted.cell_type.pred", test.use = "MAST", ident.1 = x, only.pos = T)
  y <- filter(y, y$p_val_adj < 0.05)
  if(nrow(y) == 0)
  {
    #skip iteration because of zero significant genes
    next
  }
  #rest of iteration for case of no error
  # create a 'gene' column
  y$gene = NA
  y$gene = rownames(y)
  # save dataframe
  assign( paste(x) , y) 
}

#export each data frame to separate sheets in same Excel file
name.vec <- c(celltypes_epi, celltypes_imm, celltypes_str)
df.list <- c()
for (i in 1:length(name.vec)){
  df.list <- c(df.list, as.name(name.vec[i]))
}

library("openxlsx")
wb <- createWorkbook()
for(i in 1:length(df.list)){
  print(i)
  print(df.list[[i]])
  addWorksheet(wb, name.vec[i])
  writeData(wb, name.vec[i], eval(df.list[[i]]), startRow = 1, startCol = 1)
}
saveWorkbook(wb, file = "Results/DE_2023/markers.xlsx", overwrite = TRUE)

# Create supp table 2 with pibble results
PSCIvPSCNI <- read_tsv("Results/Pibble_results/PSC_I_vs_NI/psc_i_vs_ni_table_p95.tsv") 
PSCIvPSCNI$idx <- NULL
PSCIvPSCNI <- subset(PSCIvPSCNI, PSCIvPSCNI$.width == "0.95")

UCIvUCNI <- read_tsv("Results/Pibble_results/UC_I_vs_NI/uc_i_vs_ni_table_p95.tsv") 
UCIvUCNI$idx <- NULL
UCIvUCNI <- subset(UCIvUCNI, UCIvUCNI$.width == "0.95")

UCNIvHCNI <- read_tsv("Results/Pibble_results/UC_NI_vs_HC/uc_ni_vs_hc_table_p95.tsv") 
UCNIvHCNI$idx <- NULL
UCNIvHCNI <- subset(UCNIvHCNI, UCNIvHCNI$.width == "0.95")

PSCNIvHCNI <- read_tsv("Results/Pibble_results/PSC_NI_vs_HC/psc_ni_vs_hc_table_p95.tsv") 
PSCNIvHCNI$idx <- NULL
PSCNIvHCNI <- subset(PSCNIvHCNI, PSCNIvHCNI$.width == "0.95")

pibble <- createWorkbook()
addWorksheet(pibble, "PSCIvPSCNI")
writeData(pibble, "PSCIvPSCNI", PSCIvPSCNI, startRow = 1, startCol = 1)
addWorksheet(pibble, "UCIvUCNI")
writeData(pibble, "UCIvUCNI", UCIvUCNI, startRow = 1, startCol = 1)
addWorksheet(pibble, "UCNIvHCNI")
writeData(pibble, "UCNIvHCNI", UCNIvHCNI, startRow = 1, startCol = 1)
addWorksheet(pibble, "PSCNIvHCNI")
writeData(pibble, "PSCNIvHCNI", PSCNIvHCNI, startRow = 1, startCol = 1)

saveWorkbook(pibble, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults2_pibble.xlsx", overwrite = TRUE)


# Create supp table 3 with DE analyses for PSC I vs NI
setwd("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_PSCINI_all/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_PSCINI_all/")
name_list <- str_extract(file_list, '.*(?=\\.csv)')
dataset <- data.frame()
wb <- createWorkbook()
for (i in 1:length(file_list)){
  print(name_list[i])
  print(file_list[i])
  temp_data <- read_csv(file_list[i]) 
  temp_data$...1 <- NULL
  addWorksheet(wb, name_list[i])
  writeData(wb, name_list[i], temp_data, startRow = 1, startCol = 1)
}

saveWorkbook(wb, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults3_DEgenesPSCINI.xlsx", overwrite = TRUE)

# Create supp table 4 with DE analyses for UC I vs NI
setwd("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_UCINI_all/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DE_UCINI_all/")
name_list <- str_extract(file_list, '.*(?=\\.csv)')
dataset <- data.frame()
wb <- createWorkbook()
for (i in 1:length(file_list)){
  print(name_list[i])
  print(file_list[i])
  temp_data <- read_csv(file_list[i]) 
  temp_data$...1 <- NULL
  addWorksheet(wb, name_list[i])
  writeData(wb, name_list[i], temp_data, startRow = 1, startCol = 1)
}

saveWorkbook(wb, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults4_DEgenesUCINI.xlsx", overwrite = TRUE)

# Create supp table 5 with GO terms for PSC I vs NI
  # some file names were abbreviations because of character max
setwd("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/GO_only/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/PSCIvPSCNI/GO_only/")
name_list <- str_extract(file_list, '.*(?=\\.csv)')
dataset <- data.frame()
wb <- createWorkbook()
for (i in 1:length(file_list)){
  print(name_list[i])
  print(file_list[i])
  temp_data <- read_csv(file_list[i]) 
  temp_data$...1 <- NULL
  addWorksheet(wb, name_list[i])
  writeData(wb, name_list[i], temp_data, startRow = 1, startCol = 1)
}

saveWorkbook(wb, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults5_GOtermsPSCINI.xlsx", overwrite = TRUE)

# Create supp table 6 with GO terms for UC I vs NI
# some file names were abbreviations because of character max
setwd("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/GO_only/")
file_list <- list.files(path="/Users/amberbangma/Documents/R/PSC/Results/DE_2023/UCIvUCNI/GO_only/")
name_list <- str_extract(file_list, '.*(?=\\.csv)')
dataset <- data.frame()
wb <- createWorkbook()
for (i in 1:length(file_list)){
  print(name_list[i])
  print(file_list[i])
  temp_data <- read_csv(file_list[i]) 
  temp_data$...1 <- NULL
  addWorksheet(wb, name_list[i])
  writeData(wb, name_list[i], temp_data, startRow = 1, startCol = 1)
}

saveWorkbook(wb, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults6_GOtermsUCINI.xlsx", overwrite = TRUE)

# Supp table 7 was manually created with  PSC riskgenes from Jiang et al. (2017)

# Create supp table 8 with cell-cell interaction results
setwd("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/Cellchat/")
PSCI_DUOX2sender <- read.xlsx("PSC_I_cellchat_DUOX2sender.xlsx") 
PSCI_Inflammonosender <- read.xlsx("PSC_I_cellchat_inflamMonosender.xlsx") 
UCI_Inflammonosender <- read.xlsx("UC_I_cellchat_inflamMonosender.xlsx") 
UCI_DUOX2sender <- read.xlsx("UC_I_cellchat_DUOX2sender.xlsx") 

wb <- createWorkbook()
addWorksheet(wb, "PSCI_DUOX2sender")
writeData(wb, "PSCI_DUOX2sender", PSCI_DUOX2sender, startRow = 1, startCol = 1)
addWorksheet(wb, "PSCI_Inflammonosender")
writeData(wb, "PSCI_Inflammonosender", PSCI_Inflammonosender, startRow = 1, startCol = 1)
addWorksheet(wb, "UCI_Inflammonosender")
writeData(wb, "UCI_Inflammonosender", UCI_Inflammonosender, startRow = 1, startCol = 1)
addWorksheet(wb, "UCI_DUOX2sender")
writeData(wb, "UCI_DUOX2sender", UCI_DUOX2sender, startRow = 1, startCol = 1)

saveWorkbook(wb, file = "/Users/amberbangma/Documents/R/PSC/Results/DE_2023/SuppResults8_Cellcell.xlsx", overwrite = TRUE)


