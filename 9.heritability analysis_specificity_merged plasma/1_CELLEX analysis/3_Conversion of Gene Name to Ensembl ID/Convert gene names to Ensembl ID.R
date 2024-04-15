
#####################################################################
# Conversion of Gene Name to Ensembl ID for PSC_UC CELLEX output
# Change column names to uppercase letters 
# Date: 06-04-2024
# Author: Shiqing Sun
#####################################################################


####################
# libraries        #
####################

library(patchwork)
library(biomaRt)
library(ggplot2)
library(readr)
library(tidyr)
library(readxl)
library(dplyr)
library(tibble)
options(stringsAsFactors = FALSE)

####################
# Functions        #
####################

# Define a function to generate column names in R
generate_column_names <- function(n) {
  letters <- LETTERS
  if (n <= 26) {
    return(letters[1:n])
  } else {
    # Initialize column names list with single letters
    column_names <- letters
    # Generate combinations of letters for more than 26 columns
    for (i in 1:length(letters)) {
      for (j in 1:length(letters)) {
        column_names <- c(column_names, paste0(letters[i], letters[j]))
        if (length(column_names) >= n) {
          return(column_names[1:n])
        }
      }
    }
  }
}


####################
# Main Code        #
####################

###################################################################### Conversion of Gene Name to Ensembl ID for PSC CELLEX output ######################################################################  

# we need some more memory
options(future.globals.maxSize = 2000 * 1000 * 1024^2)

# set seed
set.seed(1234)

# Read the PSC CELLEX output
PSC_output <- read.csv("/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/output/PSC_UC_gene_score_combined_Plasma.csv", check.names = FALSE)

# Select ensembl dataset
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")

# grab gene list
genes_list <- unique(PSC_output$gene)

# get Ensembl ID from ensembl dataset
genes_to_ensembl <- getBM(attributes = c('hgnc_symbol', 'ensembl_gene_id'),
                          filters = 'hgnc_symbol',
                          values = genes_list,
                          mart = ensembl)


# Create a mapping that from hgnc_symbol to ensembl_gene_id
map <- setNames(genes_to_ensembl$ensembl_gene_id, genes_to_ensembl$hgnc_symbol)

# Convert gene name of PSC_output to Ensembl ID
PSC_output$gene <- sapply(PSC_output$gene, function(x) ifelse(x %in% names(map), map[x], NA))

# Remove rows including NA
PSC_output2 <- PSC_output[!is.na(PSC_output$gene), ]


###################################################################### Change column names to uppercase letters ######################################################################  

# Calculate the number of column names to generate (excluding the first column)
n_cols <- ncol(PSC_output2) - 1

# Generate combinations of uppercase letters
# Use direct LETTERS for up to 26 columns; for more, generate combinations
new_colnames <- generate_column_names(n_cols)

# Replace column names from the second to the last with the new letter combinations
original_colnames <- colnames(PSC_output2) 
colnames(PSC_output2)[2:length(colnames(PSC_output2))] <- new_colnames[1:n_cols]

# Check the updated column names
colnames(PSC_output2)

# Save PSC_output2
write.csv(PSC_output2, "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLECT/example/PSC_UC_gene_score_merged_Plasma_modified.csv", row.names = FALSE, quote = FALSE)
write.csv(PSC_output2, "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/output/PSC_UC_gene_score_merged_Plasma_modified.csv", row.names = FALSE, quote = FALSE)


# Create a dataframe with the mapping of original to new column names
# Exclude the first column from the mapping, as it remains unchanged
column_mapping <- data.frame(
  OriginalColumnName = original_colnames[-1], # Exclude the first original column name
  NewColumnName =  new_colnames # Include the unchanged first column manually
)

# Save the mapping dataframe
write.csv(column_mapping, "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/output/column_mapping.csv", row.names = F)






