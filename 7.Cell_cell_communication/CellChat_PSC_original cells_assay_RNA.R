
#################################################################################################################
# Cellchat analysis in PSC-IBD data 
# Date: 02-11-2023
# Assay = RNA
#################################################################################################################

####################
# libraries
####################

library(Seurat)
library(Matrix)
library(patchwork)
library(ggpubr)
library(NMF)
library(ggalluvial)
library(circlize)
library(ComplexHeatmap)
library(reshape2)
library(CellChat)

####################
# Functions        
####################

# set some options
options(stringsAsFactors = FALSE)

# functions
init_cellchat_object <- function(seurat_object, assay = 'RNA', slot = 'data', ident = 'final_celltype'){
  # set the default assay
  DefaultAssay(seurat_object) <- assay
  # extract the data
  data.input <- GetAssayData(seurat_object, assay = assay, slot = slot)
  meta <- seurat_object@meta.data
  # create the object
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = ident)
  return(cellchat)
}

preprocess_cellchat_object <- function(chat_object, nthreads=8){
  # set the database
  CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
  CellChatDB.use <- CellChatDB # simply use the default CellChatDB
  chat_object@DB <- CellChatDB.use
  chat_object <- subsetData(chat_object) # This step is necessary even if using the whole database
  # set multithreading options
  future::plan("multisession", workers = nthreads) # do parallel
  # get genes
  chat_object <- identifyOverExpressedGenes(chat_object)
  chat_object <- identifyOverExpressedInteractions(chat_object)
  # project gene expression data onto PPI network (optional)
  chat_object <- projectData(chat_object, PPI.human)
  return(chat_object)
}

inference_communication_network <- function(chat_object, min.cells=10, thresh=1){
  # Compute the communication probability and infer cellular communication network
  chat_object <- computeCommunProb(chat_object)
  # Filter out the cell-cell communication if there are only few number of cells in certain cell groups
  chat_object <- filterCommunication(chat_object, min.cells = min.cells)
  # Infer the cell-cell communication at a signaling pathway level
  chat_object <- computeCommunProbPathway(chat_object, thresh = thresh)
  # Calculate the aggregated cell-cell communication network
  chat_object <- aggregateNet(chat_object, thresh = thresh)
  return(chat_object)
}


####################
# Main codes        
####################

#read the seurat object
PSC_data <- readRDS('/groups/umcg-weersma/tmp01/projects/PSC/ongoing/Seurat_object/PSC_colitis_allcomps_mergedcelltypes.rds')

# select the final celltypes
PSC_data@meta.data$final_celltype <- PSC_data@meta.data$mergedcelltype
PSC_data@meta.data$final_celltype <- gsub('\\+', 'pos', PSC_data@meta.data$final_celltype)
PSC_data@meta.data$final_celltype <- gsub('- ', 'neg_', PSC_data@meta.data$final_celltype)
PSC_data@meta.data$final_celltype <- gsub(' |-', '_', PSC_data@meta.data$final_celltype)

############################################# Cellchat in PSC_I sample ########################################################################

# Subset PSC_I sample
PSC_data_PSC_I <- PSC_data[, PSC_data@meta.data$state == "PSC-I"]

# Initialization of CellChat object
PSC_data_PSC_I_cellchat <- init_cellchat_object(PSC_data_PSC_I)

# Preprocess the expression data for cell-cell communication analysis
PSC_data_PSC_I_cellchat <- preprocess_cellchat_object(PSC_data_PSC_I_cellchat)

# Inference of cell-cell communication network
PSC_data_PSC_I_cellchat <- inference_communication_network(PSC_data_PSC_I_cellchat)

# Extract significant L-R pairs (DUXO2 as sender)
PSC_data_PSC_I_cellchat_DUOX2 <- subsetCommunication(PSC_data_PSC_I_cellchat, sources.use = "DUOX2_enterocytes")
# Sort DataFrame
PSC_data_PSC_I_cellchat_DUOX2 <- PSC_data_PSC_I_cellchat_DUOX2[order(PSC_data_PSC_I_cellchat_DUOX2$target),]
PSC_data_PSC_I_cellchat_DUOX2 <- PSC_data_PSC_I_cellchat_DUOX2[order(PSC_data_PSC_I_cellchat_DUOX2$source),]
write.csv(PSC_data_PSC_I_cellchat_DUOX2, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_PSC_I_cellchat/PSC_data_PSC_I_cellchat_DUOX2_sender.csv")

# Extract significant L-R pairs (inflammatory monocytes as sender)
PSC_data_PSC_I_cellchat_inflammatory_monocytes <- subsetCommunication(PSC_data_PSC_I_cellchat, sources.use = "Inflammatory_Monocytes")
# Sort DataFrame
PSC_data_PSC_I_cellchat_inflammatory_monocytes <- PSC_data_PSC_I_cellchat_inflammatory_monocytes[order(PSC_data_PSC_I_cellchat_inflammatory_monocytes$target),]
PSC_data_PSC_I_cellchat_inflammatory_monocytes <- PSC_data_PSC_I_cellchat_inflammatory_monocytes[order(PSC_data_PSC_I_cellchat_inflammatory_monocytes$source),]
write.csv(PSC_data_PSC_I_cellchat_inflammatory_monocytes, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_PSC_I_cellchat/PSC_data_PSC_I_cellchat_inflammatory_monocytes_sender.csv")

# Extract significant L-R pairs (CD8+IL17+ as sender)
PSC_data_PSC_I_cellchat_CD8pos_IL17pos <- subsetCommunication(PSC_data_PSC_I_cellchat, sources.use = "CD8pos_IL17pos")
# Sort DataFrame
PSC_data_PSC_I_cellchat_CD8pos_IL17pos <- PSC_data_PSC_I_cellchat_CD8pos_IL17pos[order(PSC_data_PSC_I_cellchat_CD8pos_IL17pos$target),]
PSC_data_PSC_I_cellchat_CD8pos_IL17pos <- PSC_data_PSC_I_cellchat_CD8pos_IL17pos[order(PSC_data_PSC_I_cellchat_CD8pos_IL17pos$source),]
write.csv(PSC_data_PSC_I_cellchat_CD8pos_IL17pos, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_PSC_I_cellchat/PSC_data_PSC_I_cellchat_CD8pos_IL17pos_sender.csv")

# Extract significant L-R pairs (CD4+PD1+ as sender)
PSC_data_PSC_I_cellchat_CD4pos_PD1pos <- subsetCommunication(PSC_data_PSC_I_cellchat, sources.use = "CD4pos_PD1pos")
# Sort DataFrame
PSC_data_PSC_I_cellchat_CD4pos_PD1pos <- PSC_data_PSC_I_cellchat_CD4pos_PD1pos[order(PSC_data_PSC_I_cellchat_CD4pos_PD1pos$target),]
PSC_data_PSC_I_cellchat_CD4pos_PD1pos <- PSC_data_PSC_I_cellchat_CD4pos_PD1pos[order(PSC_data_PSC_I_cellchat_CD4pos_PD1pos$source),]
write.csv(PSC_data_PSC_I_cellchat_CD4pos_PD1pos, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_PSC_I_cellchat/PSC_data_PSC_I_cellchat_CD4pos_PD1pos_sender.csv")

# Save CellChat object
saveRDS(PSC_data_PSC_I_cellchat, file = "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_PSC_I_cellchat/PSC_data_PSC_I_cellchat.rds")


############################################# Cellchat in UC-I sample ########################################################################

# Subset UC-I sample
PSC_data_UC_I <- PSC_data[, PSC_data@meta.data$state == 'UC-I']

# Initialization of CellChat object
PSC_data_UC_I_cellchat <- init_cellchat_object(PSC_data_UC_I)

# Preprocess the expression data for cell-cell communication analysis
PSC_data_UC_I_cellchat <- preprocess_cellchat_object(PSC_data_UC_I_cellchat)

# Inference of cell-cell communication network
PSC_data_UC_I_cellchat <- inference_communication_network(PSC_data_UC_I_cellchat)


# Extract significant L-R pairs (DUXO2 as sender)
PSC_data_UC_I_cellchat_DUOX2 <- subsetCommunication(PSC_data_UC_I_cellchat, sources.use = "DUOX2_enterocytes")
# Sort DataFrame
PSC_data_UC_I_cellchat_DUOX2 <- PSC_data_UC_I_cellchat_DUOX2[order(PSC_data_UC_I_cellchat_DUOX2$target),]
PSC_data_UC_I_cellchat_DUOX2 <- PSC_data_UC_I_cellchat_DUOX2[order(PSC_data_UC_I_cellchat_DUOX2$source),]
write.csv(PSC_data_UC_I_cellchat_DUOX2, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_UC_I_cellchat/PSC_data_UC_I_cellchat_DUOX2_sender.csv")

# Extract significant L-R pairs (inflammatory monocytes as sender)
PSC_data_UC_I_cellchat_inflammatory_monocytes <- subsetCommunication(PSC_data_UC_I_cellchat, sources.use = "Inflammatory_Monocytes")
# Sort DataFrame
PSC_data_UC_I_cellchat_inflammatory_monocytes <- PSC_data_UC_I_cellchat_inflammatory_monocytes[order(PSC_data_UC_I_cellchat_inflammatory_monocytes$target),]
PSC_data_UC_I_cellchat_inflammatory_monocytes <- PSC_data_UC_I_cellchat_inflammatory_monocytes[order(PSC_data_UC_I_cellchat_inflammatory_monocytes$source),]
write.csv(PSC_data_UC_I_cellchat_inflammatory_monocytes, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_UC_I_cellchat/PSC_data_UC_I_cellchat_inflammatory_monocytes_sender.csv")

# Extract significant L-R pairs (CD8+IL17+ as sender)
PSC_data_UC_I_cellchat_CD8pos_IL17pos <- subsetCommunication(PSC_data_UC_I_cellchat, sources.use = "CD8pos_IL17pos")
# Sort DataFrame
PSC_data_UC_I_cellchat_CD8pos_IL17pos <- PSC_data_UC_I_cellchat_CD8pos_IL17pos[order(PSC_data_UC_I_cellchat_CD8pos_IL17pos$target),]
PSC_data_UC_I_cellchat_CD8pos_IL17pos <- PSC_data_UC_I_cellchat_CD8pos_IL17pos[order(PSC_data_UC_I_cellchat_CD8pos_IL17pos$source),]
write.csv(PSC_data_UC_I_cellchat_CD8pos_IL17pos, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_UC_I_cellchat/PSC_data_UC_I_cellchat_CD8pos_IL17pos_sender.csv")

# Extract significant L-R pairs (CD4+PD1+ as sender)
PSC_data_UC_I_cellchat_CD4pos_PD1pos <- subsetCommunication(PSC_data_UC_I_cellchat, sources.use = "CD4pos_PD1pos")
# Sort DataFrame
PSC_data_UC_I_cellchat_CD4pos_PD1pos <- PSC_data_UC_I_cellchat_CD4pos_PD1pos[order(PSC_data_UC_I_cellchat_CD4pos_PD1pos$target),]
PSC_data_UC_I_cellchat_CD4pos_PD1pos <- PSC_data_UC_I_cellchat_CD4pos_PD1pos[order(PSC_data_UC_I_cellchat_CD4pos_PD1pos$source),]
write.csv(PSC_data_UC_I_cellchat_CD4pos_PD1pos, "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_UC_I_cellchat/PSC_data_UC_I_cellchat_CD4pos_PD1pos_sender.csv")

# Save CellChat object
saveRDS(PSC_data_UC_I_cellchat, file = "/Users/s.qs/Documents/Chapters/4 PSC_Werna/20231102_CellChat_original cells_assay RNA/Output/PSC_data_UC_I_cellchat/PSC_data_UC_I_cellchat.rds")

