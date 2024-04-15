# CELLEX.py

#Import packages
import loompy # needed for importing data for this tutorial
import numpy as np # needed for formatting data for this tutorial
import pandas as pd # needed for formatting data for this tutorial
import dask.dataframe as dd
import cellex

# Set constants
dirOut = "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/output/" # output directory for results and plots
prefixData = "PSC_UC" # prefix to prepend to files

#Load input data and metadata
path = "/groups/umcg-weersma/tmp01/projects/PSC/ongoing/heritability_analysis_specificity_merged_Plasma/CELLEX/input/"
metadata = pd.read_csv(path + "PSC_UC_metadata_merged_Plasma.csv", index_col=0)
data = pd.read_csv(path + "PSC_UC_UMI_count.csv", index_col=0)

#Create ESObject and compute Expression Specificity
eso = cellex.ESObject(data=data, annotation=metadata, verbose=True)
eso.compute(verbose=True)

# Save results
eso.save_as_csv(path=dirOut, file_prefix=prefixData, verbose=True)

# Save gene score table
eso.results["esmu"].to_csv(dirOut + "PSC_UC_gene_score_combined_Plasma.csv", index=True, header=True)


