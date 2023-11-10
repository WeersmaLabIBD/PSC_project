#######################################
# Generate DUOX2 specific markers     #
# Script to create figure 2b and 2c   #
#######################################

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
library(ggalluvial)
library(reshape2)

#load in datasets
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
DefaultAssay(epi) = "RNA"
Idents(epi) <- "celltype.final"

#Create DUOX2 marker pathways
DUOX2markers <- FindMarkers(epi, group.by = "celltype.final", test.use = "MAST", ident.1 = "DUOX2 enterocytes",only.pos = T)
DUOX2markers$Gene <- rownames(DUOX2markers)
DUOX2markers <- filter(DUOX2markers, DUOX2markers$p_val_adj < 0.05)
DUOX2pathways <- enrichr(DUOX2markers$Gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
DUOX2pathways <- filter(DUOX2pathways, DUOX2pathways$Adjusted.P.value < 0.05)

#Create enterocyte marker pathways
Enteromarkers <- FindMarkers(epi, group.by = "celltype.final", test.use = "MAST", ident.1 = "Enterocytes", only.pos = T)
Enteromarkers$Gene <- rownames(Enteromarkers)
Enteromarkers <- filter(Enteromarkers, Enteromarkers$p_val_adj < 0.05)
Enteropathways <- enrichr(Enteromarkers$Gene, databases = "GO_Biological_Process_2018")$GO_Biological_Process
Enteropathways <- filter(Enteropathways, Enteropathways$Adjusted.P.value < 0.05)

#Create specific marker pathways list and write df's
DUOX2specific <- subset(DUOX2pathways,!(Term%in%Enteropathways$Term))
DUOX2nonspecific <- subset(DUOX2pathways,(Term%in%Enteropathways$Term))
DUOX2genesspecific <- subset(DUOX2markers,!(Gene%in%Enteromarkers$Gene))
write.csv(DUOX2specific, "Results/DE_2023/DUOX2/DUOX2specific.csv")
write.csv(DUOX2nonspecific, "Results/DE_2023/DUOX2/DUOX2nonspecific.csv")
write.csv(DUOX2pathways, "Results/DE_2023/DUOX2/DUOX2pathways.csv")
write.csv(Enteropathways, "Results/DE_2023/DUOX2/Enteropathways.csv")
write.csv(DUOX2genesspecific, "Results/DE_2023/DUOX2/DUOX2genesspecific.csv")
write.csv(DUOX2markers, "Results/DE_2023/DUOX2/DUOX2markers.csv")
write.csv(Enteromarkers, "Results/DE_2023/DUOX2/Enteromarkers.csv")

# create figure 2c
  # MHC-II in DUOX2: HLA-DPA1, HLA-DQB1, HLA-DRA, HLA-DRB1
  # MHC-II in entero: none
genes = c("HLA-DPA1", "HLA-DQB1", "HLA-DRA", "HLA-DRB1")
VlnPlot(epi,features = genes, cols = c("#E8B094","#59402E99"), idents = c("DUOX2 enterocytes", "Enterocytes"), 
        split.by = "celltype.final", ncol = 4)

sel.genes = c("HLA-DRA", "HLA-DRB1")
VlnPlot(epi,features = sel.genes, cols = c("#E8B094","#59402E99"), idents = c("DUOX2 enterocytes", "Enterocytes"), 
        split.by = "celltype.final", pt.size = 0.001)

#ggsave("Results/Figures_new/vlnDUOX_HLA.pdf", width = 8, height = 6)

# create figure 2b
# get top 10 DUOX2 specific pathways and get genes per pathway
Top10Pathways <- DUOX2specific[order(DUOX2specific$Adjusted.P.value), ][1:10, ]

Pathway1 <- unlist(strsplit(Top10Pathways$Genes[1], split = ";"))
Pathway2 <- unlist(strsplit(Top10Pathways$Genes[2], split = ";"))
Pathway3 <- unlist(strsplit(Top10Pathways$Genes[3], split = ";"))
Pathway4 <- unlist(strsplit(Top10Pathways$Genes[4], split = ";"))
Pathway5 <- unlist(strsplit(Top10Pathways$Genes[5], split = ";"))
Pathway6 <- unlist(strsplit(Top10Pathways$Genes[6], split = ";"))
Pathway7 <- unlist(strsplit(Top10Pathways$Genes[7], split = ";"))
Pathway8 <- unlist(strsplit(Top10Pathways$Genes[8], split = ";"))
Pathway9 <- unlist(strsplit(Top10Pathways$Genes[9], split = ";"))
Pathway10 <- unlist(strsplit(Top10Pathways$Genes[10], split = ";"))

Pathwaylist <- c(Pathway1,Pathway2,Pathway3,Pathway4,Pathway5,Pathway6,Pathway7,Pathway8,Pathway9,Pathway10)

#check overlap and log2FC per gene
genesinpathways <- data.frame(table(c(Pathway1,Pathway2,Pathway3,Pathway4,Pathway5,Pathway6,Pathway7,Pathway8,Pathway9,Pathway10)))
genesinpathways <- merge(genesinpathways, DUOX2markers[,c(2,6)], by.x="Var1",by.y="Gene")

#select genes of interest (genes with 3 or more pathways and HLA-E)
interest <- c("STAT1","TMSB4X","ZC3H12A","CASP8","IFNGR1","IFNGR2","JAK1","TIMP1","VEGFA","HLA-E")
genesinpathways <- genesinpathways[genesinpathways$Var1 %in% interest,]

#create df with genes, log2FC and pathway ther are in
'HLA-E' %in% Pathway6
row1 <- c("STAT1",Top10Pathways[1,7],"Pathway 1")
row2 <- c("STAT1",Top10Pathways[2,7],"Pathway 2")
row3 <- c("STAT1",Top10Pathways[8,7],"Pathway 8")
row4 <- c("STAT1",Top10Pathways[9,7],"Pathway 9")
row5 <- c("TMSB4X",Top10Pathways[3,7],"Pathway 3")
row6 <- c("TMSB4X",Top10Pathways[4,7],"Pathway 4")
row7 <- c("TMSB4X",Top10Pathways[7,7],"Pathway 7")
row8 <- c("TMSB4X",Top10Pathways[8,7],"Pathway 8")
row9 <- c("ZC3H12A",Top10Pathways[1,7],"Pathway 1")
row10 <- c("ZC3H12A",Top10Pathways[3,7],"Pathway 3")
row11 <- c("ZC3H12A",Top10Pathways[8,7],"Pathway 8")
row12 <- c("ZC3H12A",Top10Pathways[10,7],"Pathway 10")
row13 <- c("CASP8",Top10Pathways[2,7],"Pathway 2")
row14 <- c("CASP8",Top10Pathways[8,7],"Pathway 8")
row15 <- c("CASP8",Top10Pathways[5,7],"Pathway 5")
row16 <- c("IFNGR1",Top10Pathways[1,7],"Pathway 1")
row17 <- c("IFNGR1",Top10Pathways[2,7],"Pathway 2")
row18 <- c("IFNGR1",Top10Pathways[9,7],"Pathway 9")
row19 <- c("IFNGR2",Top10Pathways[1,7],"Pathway 1")
row20 <- c("IFNGR2",Top10Pathways[2,7],"Pathway 2")
row21 <- c("IFNGR2",Top10Pathways[9,7],"Pathway 9")
row22 <- c("JAK1",Top10Pathways[1,7],"Pathway 1")
row23 <- c("JAK1",Top10Pathways[2,7],"Pathway 2")
row24 <- c("JAK1",Top10Pathways[9,7],"Pathway 9")
row25 <- c("TIMP1",Top10Pathways[1,7],"Pathway 1")
row26 <- c("TIMP1",Top10Pathways[4,7],"Pathway 4")
row27 <- c("TIMP1",Top10Pathways[7,7],"Pathway 7")
row28 <- c("VEGFA",Top10Pathways[1,7],"Pathway 1")
row29 <- c("VEGFA",Top10Pathways[4,7],"Pathway 4")
row30 <- c("VEGFA",Top10Pathways[7,7],"Pathway 7")
row31 <- c("HLA-E",Top10Pathways[6,7],"Pathway 6")

final.df <-  as_data_frame(t(data.frame(row1,row2,row3,row4,row5,row6,row7,row8,row9,row10,
                          row11,row12,row13,row14,row15,row16,row17,row18,row19,
                          row20,row21,row22,row23,row24,row25,row26,row27,row28,row29,row30,row31))
)
colnames(final.df) <- c("Gene", "Odds Ratio", "Pathway")
rownames(final.df) <- NULL

final.df$Pathway  <- sub("Pathway 10", Top10Pathways$Term[10], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 1", Top10Pathways$Term[1], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 2", Top10Pathways$Term[2], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 3", Top10Pathways$Term[3], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 4", Top10Pathways$Term[4], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 5", Top10Pathways$Term[5], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 6", Top10Pathways$Term[6], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 7", Top10Pathways$Term[7], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 8", Top10Pathways$Term[8], final.df$Pathway)
final.df$Pathway  <- sub("Pathway 9", Top10Pathways$Term[9], final.df$Pathway)

#order on most significant pathway
final.df$Pathway <- factor(final.df$Pathway, levels = c("cellular response to cytokine stimulus (GO:0071345)", 
                                                        "regulation of cytokine-mediated signaling pathway (GO:0001959)",
                                                        "negative regulation of cytokine secretion (GO:0050710)",
                                                        "regulated exocytosis (GO:0045055)",
                                                        "activation of cysteine-type endopeptidase activity involved in apoptotic process (GO:0006919)",
                                                        "antigen processing and presentation of endogenous peptide antigen (GO:0002483)", 
                                                        "platelet degranulation (GO:0002576)",
                                                        "negative regulation of I-kappaB kinase/NF-kappaB signaling (GO:0043124)",
                                                        "regulation of response to interferon-gamma (GO:0060330)",
                                                        "cellular response to organic cyclic compound (GO:0071407)"))
#order on largest log2FC
final.df$Gene <- factor(final.df$Gene, levels = c("HLA-E", "VEGFA", "IFNGR1","CASP8","ZC3H12A","IFNGR2","TMSB4X", "TIMP1","STAT1","JAK1"))


final.df$`Odds Ratio` <- as.numeric(final.df$`Odds Ratio`)
final.df$ORgroup <- cut(final.df$`Odds Ratio`, breaks = c(0,5,10,15,20,Inf),
                         labels = c("<5","5-10","10-15","15-20",">20"))
final.df$freq <- 1



# create alluvial plot
ggplot(data = final.df,
       aes(axis1 = Pathway, axis2 = Gene,
           y = freq)) +
  scale_x_discrete(limits = c("Pathway", "Gene"), expand = c(.2, .05)) +
  xlab("Demographic") +
  geom_alluvium(aes(fill = ORgroup)) +
  geom_stratum() +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  theme_minimal()+NoLegend()


ggsave("/Users/amberbangma/Documents/R/PSC/Results/DE_2023/DUOX2/alluvial.pdf",  width=6, height=6)
