# scRNAseqTools
Tools for scRNAseq exam project (gene filtering, PCA, clustering, annotation)
This R package contains the  analysis for single-cell RNA-seq:

- Preprocessing and gene filtering
- Gene expression quantification (UMI ≥ 3)
- PCA and variance analysis
- UMAP visualization and clustering
- Cell type annotation via SingleR
- Tissue origin inference

scRNAseqTools is an R package for a complete analysis of single-cell RNA-seq (scRNA-seq) data from 10X Genomics. It integrates preprocessing, gene filtering, visualization, PCA, UMAP, clustering, and automatic cell type annotation using Seurat, SingleR, and other standard tools.

## Installation
clone the GitHub repository:
git clone https://github.com/chiaravit/scRNAseqTools.git

open Rstudio and install:
# If devtools is not installed
install.packages("devtools")

# Install the local package
devtools::install("/Users/chiara/Desktop/scRNAseqTools")


## Features
Data Loading: Read 10X data and GTF annotation.
Gene Filtering: Maintain protein-coding genes and then exclude mitochondrial and ribosomal genes and ribosomal pseudogenes.
Gene Expression Summary: Visualize genes with expression ≥3 UMIs as a violin plot
Dimensionality Reduction: PCA + UMAP with customizable PCs.
Clustering: Community detection via SNN.
Cell type annotation: Automatic cell type labeling using SingleR.
Vignette: Complete walkthrough of usage.

## Vignettes
The package includes a detailed vignette:
browseVignettes("scRNAseqTools")
or view it manually from doc/ folder
