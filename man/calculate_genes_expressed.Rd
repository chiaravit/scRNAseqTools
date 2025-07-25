% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/how_many_expressed_genes_umi.R
\name{calculate_genes_expressed}
\alias{calculate_genes_expressed}
\title{Calculate number of genes expressed per cell (≥3 UMIs)}
\usage{
calculate_genes_expressed(seurat_obj)
}
\arguments{
\item{seurat_obj}{A Seurat object.}
}
\value{
A vector with the number of expressed genes per cell.
}
\description{
For each cell in a Seurat object, counts how many genes have ≥3 UMIs.
create a matrix with TRUE/FALSE for genes with ≥3 UMIs and count how many are them
}
