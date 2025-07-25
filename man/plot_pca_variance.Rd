% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/histogram_variance_first_20PCs.R
\name{plot_pca_variance}
\alias{plot_pca_variance}
\title{Plot PCA histogram of first 20 PCs to assess variance}
\usage{
plot_pca_variance(seurat_obj, num_pcs = 20)
}
\arguments{
\item{seurat_obj}{A Seurat object containing PCA results}

\item{num_pcs}{Integer. The number of top principal components to plot. Default is 20.}
}
\value{
A ggplot object showing the variance explained per PC
}
\description{
This function creates a bar plot showing the proportion of variance explained
by the first 20 principal components (PCs) from a Seurat object.
Extract standard deviations of PCs, calculate variance explained.
Create a data frame of variance explained. (df) keeps only first 20 PCs.
}
\details{
ggplot2 library creates the histogram.
}
\examples{
to use the function: load ggplot2 library and plot_pca_variance(seurat_obj)
}
