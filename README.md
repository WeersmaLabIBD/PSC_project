# PSC_project

The repository for the code used in the PSC project and manuscript called "Single-cell RNA sequencing in inflammatory bowel disease in patients with primary sclerosing cholangitis: a distinct form of colitis".

The steps for replicating the results of this study are grouped into categories, and each in its own sub-directory including a README file.

# Overview 

For each different step of the scRNAseq analysis, separate codes were generated. Languages and packages are listed below: 

- R >= 3.6.1
- Seurat >= 3.1


External tools used were:
- Souporcell v.1: https://github.com/wheaton5/souporcell
- Cellranger v.3.0.2: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/installation
- fido v.1.0.2: https://jsilve24.github.io/fido/articles/introduction-to-fido.html
- Cellchat v.1.4.0: http://www.cellchat.org/


# Workflow
the content for each step includes these steps in this order:

-   Alignment of the sequence data to the GRCh38 reference Human Genome: 
-   Lane integration:
-   Demultiplexing, doublet_detection, and quality control (QC): 
-   General quality control and filtering: 
-   Normalization and data integration: 
-   Dimensional reduction and clustering: 
-   Cell type classification: 
-   Cell abundances analysis:
-   Differential gene expression analysis and pathway analysis: 
-   Cell-cell interaction analysis: 

# Hardware
Analyses were performed on the Gearshift cluster http://docs.gcc.rug.nl/gearshift/cluster/ 
