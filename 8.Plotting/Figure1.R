####################################################
# Script to generate figure 1                      #
####################################################
# labels were slightly edited after, label specification can be found in supplementary file 1

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(readr)
library(tidyr)
library(readxl)
library(openxlsx)
library(reshape2)

# load in dataset
epi <- readRDS("Nieuw/epi_azimuth_duox2.rds")
imm <- readRDS("Nieuw/imm_azimuth_with_plasma.rds")
str <- readRDS("Nieuw/stro_azimuth.rds")
all <- readRDS("Nieuw/PSC_colitis_allcomps_mergedcelltypes.rds")

DefaultAssay(epi) = "RNA"
DefaultAssay(imm) = "RNA"
DefaultAssay(str) = "RNA"
DefaultAssay(all) = "RNA"

Idents(epi) <- "celltype.final"
Idents(imm) <- "predicted.cell_type.pred"
Idents(str) <- "predicted.cell_type.pred"
Idents(all) <- "mergedcelltype"

# make umaps per subtype 

DimPlot(str, pt.size = 0,   cols = c("#FF5733", "#D95F0E", "#B38CC6",
                                    "#8c564b", "#FADADD", "#7f7f7f", "#d62728", "#FF6F61",
                                    "#aec7e8", "#ffbb78","#C7A995", "#98df8a", "#ff7f0e")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Stromal cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

#ggsave("Results/Figures_new/Dimplot_str.pdf", width = 8, height = 6)

DimPlot(str, pt.size = 0,   cols = c("#FF5733", "#D95F0E", "#B38CC6",
                                     "#8c564b", "#FADADD", "#7f7f7f", "#d62728", "#FF6F61",
                                     "#aec7e8", "#ffbb78","#C7A995", "#98df8a", "#ff7f0e")) + 
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.text.x = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.ticks = element_blank())
#ggsave("Results/Figures_new/Dimplot_str_notext.pdf", width = 6, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#1f77b4", "#aec7e8", "#800080","#6baed6", "#08519c", "#555599",
                                    "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                                    "#8c564b", "#c49c94", "#e377c2","#5E00EB", "#f7b6d2", "#7f7f7f",
                                    "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
                                    "#FDE1E2", "#ffbb78", "#ff7f0e", "#2ca02c")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Immune cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

#ggsave("Results/Figures_new/Dimplot_imm.pdf", width = 8, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#1f77b4", "#aec7e8", "#800080","#6baed6", "#08519c", "#555599",
                                     "#d62728", "#ff9896", "#9467bd", "#c5b0d5",
                                     "#8c564b", "#c49c94", "#e377c2","#5E00EB", "#f7b6d2", "#7f7f7f",
                                     "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5",
                                     "#FDE1E2", "#ffbb78", "#ff7f0e", "#2ca02c")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                         axis.title.y = element_blank(),
                                                                                         axis.text.x = element_blank(),
                                                                                         axis.text.y = element_blank(),
                                                                                         axis.ticks = element_blank())

#ggsave("Results/Figures_new/Dimplot_imm_notext.pdf", width = 6, height = 6)

DimPlot(epi, pt.size = 0,   cols = c("#2D6100", "#d9ef8b", "#ffffbf",
                                     "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#E8E23A", "#66bd63", "#1f9915",
                                     "#1f77b4", "#bcbd22", "#b2df8a", "#8c564b", "#5E00EB", "#9e0142")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) + custom_theme +
  ggtitle("Immune cells") +   xlab("UMAP 1") +
  ylab("UMAP 2")

#ggsave("Results/Figures_new/Dimplot_epi.pdf", width = 8, height = 6)



DimPlot(epi, pt.size = 0,   cols = c("#2D6100", "#d9ef8b", "#ffffbf",
                                     "#fee08b", "#fdae61", "#f46d43", "#d53e4f", "#E8E23A", "#66bd63", "#1f9915",
                                     "#1f77b4", "#bcbd22", "#b2df8a", "#8c564b", "#5E00EB", "#9e0142")) + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                                                             axis.title.y = element_blank(),
                                                                                                                             axis.text.x = element_blank(),
                                                                                                                             axis.text.y = element_blank(),
                                                                                                                             axis.ticks = element_blank())

#ggsave("Results/Figures_new/Dimplot_epi_notext.pdf", width = 6, height = 6)

# make umaps colored by inflammation
Idents(imm) <- "inflammation"
DimPlot(imm, pt.size = 0,   cols = c("#1f77b4","grey")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) 
#ggsave("Results/Figures_new/Dimplot_infl_imm.pdf", width = 8, height = 6)

DimPlot(imm, pt.size = 0,   cols = c("#1f77b4","grey")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))  +NoLegend()+ theme(axis.title.x = element_blank(),
                                                                             axis.title.y = element_blank(),
                                                                             axis.text.x = element_blank(),
                                                                             axis.text.y = element_blank(),
                                                                             axis.ticks = element_blank()) 
#ggsave("Results/Figures_new/Dimplot_infl_imm_notext.pdf", width = 6, height = 6)

Idents(epi) <- "inflammation"
DimPlot(epi, pt.size = 0,   cols = c("grey","#1f9915")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) 
#ggsave("Results/Figures_new/Dimplot_infl_epi.pdf", width = 8, height = 6)
DimPlot(epi, pt.size = 0,   cols = c("grey","#1f9915")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))  + NoLegend() + theme(axis.title.x = element_blank(),
                                                                                          axis.title.y = element_blank(),
                                                                                          axis.text.x = element_blank(),
                                                                                          axis.text.y = element_blank(),
                                                                                          axis.ticks = element_blank()) 
                                                                                           
#ggsave("Results/Figures_new/Dimplot_infl_epi_notext.pdf", width = 6, height = 6)

Idents(str) <- "inflammation"
DimPlot(str, pt.size = 0,   cols = c("#d62728","grey")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1)) 
#ggsave("Results/Figures_new/Dimplot_infl_str.pdf", width = 8, height = 6)

DimPlot(str, pt.size = 0,   cols = c("#d62728","grey")) +
  guides(color = guide_legend(override.aes = list(size=3), ncol=1))+ NoLegend() + theme(axis.title.x = element_blank(),
                                                                           axis.title.y = element_blank(),
                                                                           axis.text.x = element_blank(),
                                                                           axis.text.y = element_blank(),
                                                                           axis.ticks = element_blank()) 
#ggsave("Results/Figures_new/Dimplot_infl_str_notext.pdf", width = 6, height = 6)

# stacked barplots
# create groups
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("IgA", "IgG", "IgM", "Ig_negative")] <- "Plasma cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Follicular", "GC", "Cycling B")] <- "B cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("CD4+ Activated Fos-hi", "CD4+ Activated Fos-lo", "CD4+ Memory", "CD4+PD1+", "Tregs" )] <- "CD4+ T cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("CD8+ IELs", "CD8+ IL17+", "CD8+ LP" )] <- "CD8+ T cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Cycling T")] <- "Cycling T"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("NKs")] <- "NK cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("ILCs")] <- "ILCs"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("CD69- Mast", "CD69+ Mast","Cycling Monocytes","DC1","DC2","Inflammatory Monocytes","Macrophages" )] <- "Myeloid cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Endothelial", "Microvascular","Pericytes","Post-capillary Venules")] <- "Endothelial cells"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Myofibroblasts","RSPO3+","WNT2B+ Fos-hi","WNT2B Fos-lo 1", "WNT2B Fos-lo 2","WNT5B+ 1","WNT5B+ 2")] <- "Fibroblasts"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Inflammatory Fibroblasts")] <- "Inflammatory Fibroblasts"
all@meta.data$groups[all@meta.data$mergedcelltype %in% c("Glia")] <- "Glial cells"
all@meta.data$groups[all@meta.data$celltype.final %in% c("Best4+ Enterocytes", "Enterocytes")] <- "Enterocytes"
all@meta.data$groups[all@meta.data$celltype.final %in% c("DUOX2 enterocytes")] <- "DUOX2+ enterocytes"
all@meta.data$groups[all@meta.data$celltype.final %in% c("Cycling TA", "Secretory TA", "Stem", "TA 1", "TA 2")] <- "Undifferentiated"
all@meta.data$groups[all@meta.data$celltype.final %in% c("Enterocyte Progenitors", "Immature Enterocytes 1", "Immature Enterocytes 2")] <- "Immature absorptive"
all@meta.data$groups[all@meta.data$celltype.final %in% c("Enteroendocrine", "Goblet", "Immature Goblet", "M cells", "Tuft")] <- "Secretory"
table(all$groups)

# make functions
get_nr_cells_per_group <- function(seurat_object, groups, remove_empty=F){
  # get the numbers table
  numbers_table <- data.frame(table(seurat_object@meta.data[, groups]))
  # remove entries that have no numbers
  if (remove_empty) {
    numbers_table <- numbers_table[numbers_table$Freq > 0, ]
  }
  return(numbers_table)
}

create_barplot_percentages <- function(numbers_table, number_column, stack_column, group_columns, stacks_to_plot=NULL, group_order=NULL){
  # add a new column that has the groups
  numbers_table[['group']] <- apply(numbers_table, 1, function(x){
    # get the value of each column
    group_parts <- rep(NA, times = length(group_columns))
    for(i in 1:length(group_columns)){
      group_parts[i] <- x[group_columns[i]]
    }
    # and paste it together
    return(paste(group_parts, collapse = '_'))
  })
  # add a percentage
  numbers_table[['group_pct']] <- apply(numbers_table, 1, function(x){
    # get the group for this row
    group <- x['group']
    # and the number
    number <- as.numeric(x[number_column])
    # get the total of the group
    group_total <- sum(as.numeric(numbers_table[numbers_table[['group']] == group, number_column]))
    # calculate what percentage this entry makes up of the total
    pct <- number / group_total
    return(pct)
  })
  # subset if we don't want to plot all groups
  if (!is.null(stacks_to_plot)) {
    numbers_table <- numbers_table[numbers_table[[stack_column]] %in% stacks_to_plot, ]
  }
  # order the stack
  if (!is.null(group_order)) {
    numbers_table[['group']] <- factor(numbers_table[['group']], levels = group_order)
  }
  # grab the colours
  cc <- get_color_coding_dict()
  # subset to what we have
  colour_map <- cc[intersect(names(cc), numbers_table[[stack_column]])]
  # make the fill
  fillScale <- scale_fill_manual(name = "cell type",values = colour_map)
  # create the plot
  p <- ggplot(data = numbers_table, mapping = aes(x = numbers_table[['group']], y = numbers_table[['group_pct']], fill = numbers_table[[stack_column]])) +
    geom_bar(stat = 'identity', position = 'stack') +
    fillScale +
    xlab('group') +
    ylab('fraction') +
    theme(panel.border = element_rect(color="black", fill=NA, size=1.1), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), strip.background = element_rect(colour="white", fill="white"))
  return(p)
}

get_color_coding_dict <- function(){
  color_coding <- list()
  #ecco box
  color_coding[["Enterocytes"]] <- "#2D6100"
  color_coding[["DUOX2+ enterocytes"]] <- "#E8B094"
  color_coding[["Immature absorptive"]] <- "#d9ef8b"
  color_coding[["Secretory"]] <- "#fdae61"
  color_coding[["Undifferentiated"]] <- "#ffec80"
  color_coding[["Myeloid cells"]] <- "#FDE1E2"
  color_coding[["Plasma cells"]] <- "#08519c"
  color_coding[["B.cells"]] <- "#17becf"
  color_coding[["B cells"]] <- "#17becf"
  color_coding[["CD8+ T cells"]] <- "#1f77b4"
  color_coding[["CD4+ T cells"]] <- "#9edae5"
  color_coding[["Cycling T"]] <- "#800080"
  color_coding[["NK cells"]] <- "#ff9896"
  color_coding[["ILCs"]] <- "#f46d43"
  color_coding[["Inflammatory Fibroblasts"]] <- "#c5b0d5"
  color_coding[["Fibroblasts"]] <- "#d62728"
  color_coding[["Glial cells"]] <- "#ff7f0e"
  color_coding[["Endothelial cells"]] <- "#ffbb78"
  return(color_coding)}


# get nr of cells per celltype
cell_numbers <- get_nr_cells_per_group(all, c('groups', 'state', 'compartment'), remove_empty = T)

cell_numbers$groups <- factor(cell_numbers$groups, levels = c("Enterocytes","DUOX2+ enterocytes", "Immature absorptive", "Secretory","Undifferentiated",
                                                              "Myeloid cells", "Plasma cells", "B cells", "CD8+ T cells", "CD4+ T cells", "Cycling T", "NK cells", "ILCs",
                                                              "Inflammatory Fibroblasts","Fibroblasts","Glial cells","Endothelial cells"))
cell_numbers$state <- factor(cell_numbers$state)
cell_numbers$state <- factor(cell_numbers$state, levels = c("HC-NI", "UC-NI", "UC-I", "PSC-NI", "PSC-I"))


immune_cells <- as.character(unique(cell_numbers[cell_numbers$compartment == 'immune', 'groups']))
epithelial_cells <- as.character(unique(cell_numbers[cell_numbers$compartment == 'epithelial', 'groups']))
stromal_cells <- as.character(unique(cell_numbers[cell_numbers$compartment == 'stromal', 'groups']))

# plot
create_barplot_percentages(cell_numbers, number_column = 'Freq', stack_column = 'groups', group_columns = c('state'), group_order = c('HC-NI', 'UC-NI', 'UC-I', 'PSC-NI', 'PSC-I'), stacks_to_plot = immune_cells) + NoLegend()
#ggsave("Results/Figures_new/Bar_imm.pdf", width = 4, height = 6)

create_barplot_percentages(cell_numbers, number_column = 'Freq', stack_column = 'groups', group_columns = c('state'), group_order = c('HC-NI', 'UC-NI', 'UC-I', 'PSC-NI', 'PSC-I'), stacks_to_plot = epithelial_cells) +NoLegend()
#ggsave("Results/Figures_new/Bar_epi.pdf", width = 4, height = 6)

create_barplot_percentages(cell_numbers, number_column = 'Freq', stack_column = 'groups', group_columns = c('state'), group_order = c('HC-NI', 'UC-NI', 'UC-I', 'PSC-NI', 'PSC-I'), stacks_to_plot = stromal_cells) + NoLegend()
#ggsave("Results/Figures_new/Bar_str.pdf", width = 4, height = 6)

# boxplots
  #create counts per celltype per sample
data=all
counts=as.data.frame.matrix(table(data$Final_HTO, data$mergedcelltype))

  #calculate proportions
total_row = apply(counts, 1, sum)
pcts = lapply(counts, function(x) {
  x / total_row
})
prop = as.data.frame(pcts)
rownames(counts) = rownames(prop)

prop$disease =  sapply(strsplit(rownames(prop),"-"), `[`, 1)
prop$inflammation =  sapply(strsplit(rownames(prop),"-"), `[`, 2)
prop$state =  apply(prop[ ,c(55,56)] , 1 , paste , collapse = "-" )
prop = prop[,-c(55,56)]
prop$state = as.factor(prop$state)
prop$state <- factor(prop$state, levels = c("HC-NI","UC-NI","UC-I","PSC-NI","PSC-I"))
levels(prop$state)

  #boxplots
ggplot(prop, aes(x = prop$state, y = prop$Stem, fill = prop$state)) +
  geom_boxplot() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("") +
  scale_fill_manual(values = c("grey","#d9ef8b","#E8E23A","#1f9915","#2D6100")) +
  labs(fill = "Variables") +  # Legend title
  guides(fill = guide_legend(override.aes = list(colour = NULL)))+NoLegend()
#ggsave("Results/Figures_new/Box_stem.pdf", width = 4, height = 6)

ggplot(prop, aes(x = prop$state, y = prop$CD4..PD1.)) +
  geom_boxplot(fill = c("grey","#9edae5","#17becf","#1f77b4","#08519c")) +
  theme_bw() +
  theme(legend.position = "none", panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("") +
  labs(fill = "Variables") +  # Legend title
  guides(fill = guide_legend(override.aes = list(colour = NULL)))
#ggsave("Results/Figures_new/Box_cd4pd1.pdf", width = 4, height = 6)

ggplot(prop, aes(x = prop$state, y = prop$Inflammatory.Fibroblasts,fill = state)) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey", "#fee08b", "#f46d43", "#ff9896","#d62728")) +
  theme_bw() +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")) +
  xlab("") +
  ylab("") +
  labs(fill = "Variables") +  # Legend title
  guides(fill = guide_legend(override.aes = list(colour = NULL)))+NoLegend()
#ggsave("Results/Figures_new/Box_inflamfibro.pdf", width = 4, height = 6)

#DUOX2 boxplot for figure 2
levels(prop$state)
ggplot(prop, aes(x = prop$state, y = prop$DUOX2.enterocytes, fill=prop$state )) +
  geom_boxplot() +
  scale_fill_manual(values = c("grey","#FDE1E2","#c49c94","#C7A995","#59402E99")) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),legend.position = "none") +
  xlab("") + ylab("")

#ggsave("Results/Figures_new/Box_DUOX2.pdf", width = 6, height = 6)
